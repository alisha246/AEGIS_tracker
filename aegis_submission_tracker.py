#!/usr/bin/env python3

import argparse
import csv
import io
import json
import os
import re
import smtplib
import time
import threading
import xml.etree.ElementTree as ET
import shutil
import subprocess
from collections import Counter, defaultdict
from datetime import date, datetime
from email.message import EmailMessage
from urllib.parse import urljoin

import numpy as np
import requests


# ------------------------------------------------------------
# Configuration
# ------------------------------------------------------------

ENA_API = "https://www.ebi.ac.uk/ena/portal/api/search"
ENA_FILEREPORT = "https://www.ebi.ac.uk/ena/portal/api/filereport"
BIOSAMPLES_API = "https://www.ebi.ac.uk/biosamples/samples"

# AEGIS umbrella
DEFAULT_UMBRELLA = "PRJEB80366"

# Default weekly recipient
DEFAULT_EMAIL = "ahamed@ebi.ac.uk"

# NCBI
NCBI_BIOPROJECT_URL = "https://www.ncbi.nlm.nih.gov/bioproject/{}"
NCBI_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
NCBI_DATASETS_BASE = "https://api.ncbi.nlm.nih.gov/datasets/v2"
NCBI_API_KEY = os.getenv("NCBI_API_KEY")
NCBI_TOOL = "aegis_submission_tracker"
NCBI_EMAIL = DEFAULT_EMAIL
NCBI_MIN_INTERVAL = 1.0 if not NCBI_API_KEY else 0.20
NCBI_MAX_RETRIES = 6
NCBI_ESUMMARY = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"

# Keep ENA OR queries reasonably small
QUERY_CHUNK_SIZE = 40

SAMPLE_FIELDS = [
    "accession",
    "secondary_sample_accession",
    "sample_alias",
    "project_name",
    "study_accession",
    "scientific_name",
    "tax_id",
    "first_public",
    "last_updated",
]

# ------------------------------------------------------------
# Command-line arguments
# ------------------------------------------------------------

parser = argparse.ArgumentParser(
    description=(
        "AEGIS ENA submission tracker. "
        "Includes records where ENA sample project_name is AEGIS, "
        "BioSamples project/project name is AEGIS, OR "
        "the project belongs to umbrella PRJEB80366."
    )
)

parser.add_argument(
    "-u",
    "--umbrella",
    default=DEFAULT_UMBRELLA,
    help=f"Umbrella project accession (default: {DEFAULT_UMBRELLA})"
)

parser.add_argument(
    "-o",
    "--output-prefix",
    default="aegis",
    help="Prefix for output files (default: aegis)"
)

parser.add_argument(
    "--username",
    help="Optional ENA/Webin username"
)

parser.add_argument(
    "--password",
    help="Optional ENA/Webin password"
)

parser.add_argument(
    "--send-email",
    action="store_true",
    help="Send the summary email"
)

parser.add_argument(
    "--email-to",
    action="append",
    help=(
        "Recipient email. Repeat for multiple recipients. "
        f"Default: {DEFAULT_EMAIL}"
    )
)

parser.add_argument(
    "--smtp-host",
    default=os.getenv("SMTP_HOST"),
    help="SMTP host, or set SMTP_HOST"
)

parser.add_argument(
    "--smtp-port",
    type=int,
    default=int(os.getenv("SMTP_PORT", "587")),
    help="SMTP port (default: 587)"
)

parser.add_argument(
    "--smtp-user",
    default=os.getenv("SMTP_USER"),
    help="SMTP username, or set SMTP_USER"
)

parser.add_argument(
    "--smtp-password",
    default=os.getenv("SMTP_PASSWORD"),
    help="SMTP password, or set SMTP_PASSWORD"
)

parser.add_argument(
    "--smtp-from",
    default=os.getenv("SMTP_FROM", DEFAULT_EMAIL),
    help=f"From address (default: {DEFAULT_EMAIL})"
)

parser.add_argument(
    "--smtp-no-tls",
    action="store_true",
    help="Do not use STARTTLS"
)

parser.add_argument(
    "--skip-ncbi",
    action="store_true",
    help="Skip all NCBI public-status checking"
)

parser.add_argument(
    "--skip-ncbi-runs",
    action="store_true",
    help=(
        "Skip NCBI SRA run checks but still check assemblies in NCBI. "
        "Useful on shared HPC systems to avoid NCBI E-utilities rate limits."
    )
)

parser.add_argument(
    "--stale-days",
    type=int,
    default=30,
    help=(
        "Flag a sample as a potential publication delay when its ENA/BioSamples "
        "public date is at least this many days old but no public run or assembly "
        "is detected (default: 30)."
    )
)

args = parser.parse_args()


# ------------------------------------------------------------
# General helpers
# ------------------------------------------------------------

def chunks(values, size):
    values = list(values)

    for i in range(0, len(values), size):
        yield values[i:i + size]


def make_auth():
    if args.username and args.password:
        return (args.username, args.password)

    return None


def parse_date(value):
    if not value:
        return None

    value = value.strip()

    try:
        return datetime.strptime(value[:10], "%Y-%m-%d").date()
    except ValueError:
        return None


def days_since(value):
    parsed = parse_date(value)

    if not parsed:
        return ""

    return (date.today() - parsed).days


# ------------------------------------------------------------
# ENA Portal API
# ------------------------------------------------------------


def get_ena_return_fields(result):
    """
    Return the set of fields currently supported by ENA for a result type.
    """
    try:
        response = requests.get(
            "https://www.ebi.ac.uk/ena/portal/api/returnFields",
            params={"result": result},
            auth=make_auth(),
            timeout=60
        )
        response.raise_for_status()

        reader = csv.DictReader(
            io.StringIO(response.text),
            delimiter="\t"
        )

        fields = set()

        for row in reader:
            field = (
                row.get("columnId")
                or row.get("field")
                or row.get("id")
                or ""
            )
            if field:
                fields.add(field)

        return fields

    except requests.RequestException:
        return set()


def query_ena(result, query, fields):
    params = {
        "result": result,
        "query": query,
        "fields": ",".join(fields),
        "format": "tsv",
        "limit": "0"
    }

    last_error = None

    for attempt in range(1, 6):
        try:
            response = requests.get(
                ENA_API,
                params=params,
                auth=make_auth(),
                timeout=180
            )

            if response.status_code in {500, 502, 503, 504}:
                last_error = requests.HTTPError(
                    f"{response.status_code} Server Error"
                )
                if attempt < 5:
                    wait_seconds = min(30, 2 ** attempt)
                    print(
                        f"\nENA temporary error {response.status_code}; "
                        f"retrying in {wait_seconds}s "
                        f"[attempt {attempt}/5] ..."
                    )
                    time.sleep(wait_seconds)
                    continue

            response.raise_for_status()

            if not response.text.strip():
                return []

            reader = csv.DictReader(
                io.StringIO(response.text),
                delimiter="\t"
            )
            return list(reader)

        except requests.RequestException as error:
            last_error = error
            if attempt < 5:
                wait_seconds = min(30, 2 ** attempt)
                print(
                    f"\nENA request failed; retrying in "
                    f"{wait_seconds}s [attempt {attempt}/5] ..."
                )
                time.sleep(wait_seconds)
                continue

    raise last_error


def get_biosamples_aegis():
    """
    Fetch BioSamples records where the exact characteristic "project name"
    has the value AEGIS, matching common capitalization variants.

    BioSamples server-side attribute filtering can be case-sensitive, so
    query AEGIS, Aegis and aegis separately, merge by accession, and then
    verify the returned characteristic locally using casefold().
    """

    print()
    print(
        'Finding BioSamples samples where "project name" = AEGIS '
        "(case-insensitive) ..."
    )

    headers = {
        "Accept": "application/hal+json"
    }

    samples = {}

    filters = [
        "attr:project name:AEGIS",
        "attr:project name:Aegis",
        "attr:project name:aegis",
    ]

    for biosamples_filter in filters:

        print(
            f"  BioSamples filter: {biosamples_filter}"
        )

        url = BIOSAMPLES_API

        params = {
            "filter": biosamples_filter,
            "size": 200
        }

        while url:

            response = requests.get(
                url,
                headers=headers,
                params=params,
                timeout=180
            )

            if not response.ok:
                print()
                print("BioSamples query failed")
                print("-----------------------")
                print("URL:", response.url)
                print("HTTP status:", response.status_code)
                print("Response:")
                print(response.text[:2000])
                print()

            response.raise_for_status()

            data = response.json()

            for sample in (
                data.get("_embedded", {})
                .get("samples", [])
            ):

                characteristics = sample.get(
                    "characteristics",
                    {}
                )

                values = characteristics.get(
                    "project name",
                    []
                )

                is_aegis = False

                for value in values:

                    if isinstance(value, dict):
                        value_text = value.get(
                            "text",
                            ""
                        )
                    else:
                        value_text = str(value)

                    if (
                        value_text
                        .strip()
                        .casefold()
                        == "aegis"
                    ):
                        is_aegis = True
                        break

                if not is_aegis:
                    continue

                accession = sample.get(
                    "accession",
                    ""
                ).strip()

                if not accession:
                    continue

                samples[accession] = {
                    "accession":
                        accession,

                    "biosamples_name":
                        sample.get(
                            "name",
                            ""
                        ),

                    "biosamples_release":
                        sample.get(
                            "release",
                            ""
                        ),

                    "biosamples_update":
                        sample.get(
                            "update",
                            ""
                        ),

                    "project_name":
                        "AEGIS"
                }

            next_link = (
                data.get("_links", {})
                .get("next", {})
                .get("href")
            )

            if next_link:
                url = urljoin(
                    "https://www.ebi.ac.uk",
                    next_link
                )
                params = None
            else:
                url = None

    print(
        f"Found {len(samples)} unique BioSamples samples "
        'with "project name" = AEGIS/Aegis/aegis'
    )

    return samples
# ------------------------------------------------------------
# Step 1
# Find projects underneath PRJEB80366
# ------------------------------------------------------------

def get_umbrella_projects(umbrella):
    """
    Resolve child PRJEB projects directly from the ENA umbrella project XML.

    ENA umbrella projects store their children as:
        <RELATED_PROJECTS>
            <RELATED_PROJECT>
                <CHILD_PROJECT accession="PRJEB..."/>
            </RELATED_PROJECT>
        </RELATED_PROJECTS>

    The live source is the ENA Browser API XML record for the umbrella.
    The last successful ENA-derived child-project list is cached so a
    temporary ENA connection failure does not abort the weekly tracker.
    """

    print()
    print(f"Finding projects under ENA umbrella {umbrella} ...")

    cache_file = f"{args.output_prefix}_umbrella_projects.txt"
    url = f"https://www.ebi.ac.uk/ena/browser/api/xml/{umbrella}"

    headers = {
        "Accept": "application/xml,text/xml,*/*;q=0.8",
        "User-Agent": (
            "AEGIS-ENA-submission-tracker/7.0 "
            f"({DEFAULT_EMAIL})"
        )
    }

    projects = set()
    last_error = None

    for attempt in range(1, 6):
        try:
            response = requests.get(
                url,
                params={"includeLinks": "true"},
                auth=make_auth(),
                timeout=180,
                headers=headers
            )

            if response.status_code in {500, 502, 503, 504}:
                raise requests.HTTPError(
                    f"ENA Browser API returned {response.status_code}"
                )

            response.raise_for_status()

            root = ET.fromstring(response.text)

            for element in root.iter():
                local_name = element.tag.split("}")[-1]

                if local_name != "CHILD_PROJECT":
                    continue

                accession = str(
                    element.attrib.get("accession", "")
                ).strip()

                if re.fullmatch(r"PRJEB\d+", accession):
                    projects.add(accession)

            # Fallback for any XML serialization/schema variation.
            if not projects:
                projects.update(
                    re.findall(
                        r'<CHILD_PROJECT[^>]+accession=["\'](PRJEB\d+)["\']',
                        response.text,
                        flags=re.IGNORECASE
                    )
                )

            projects.discard(umbrella)

            if projects:
                break

            last_error = RuntimeError(
                "ENA umbrella XML was retrieved, but no CHILD_PROJECT "
                "PRJEB accessions were found."
            )

        except (
            requests.RequestException,
            ET.ParseError
        ) as error:
            last_error = error

        if attempt < 5:
            wait_seconds = min(30, 2 ** attempt)
            print(
                f"  ENA umbrella lookup attempt {attempt}/5 failed; "
                f"retrying in {wait_seconds}s ..."
            )
            time.sleep(wait_seconds)

    if not projects and os.path.exists(cache_file):
        print(
            f"  live ENA umbrella lookup failed; using cached "
            f"ENA child-project list from {cache_file}"
        )

        with open(
            cache_file,
            "r",
            encoding="utf-8"
        ) as handle:
            projects = {
                line.strip()
                for line in handle
                if re.fullmatch(
                    r"PRJEB\d+",
                    line.strip()
                )
            }

        projects.discard(umbrella)

    if not projects:
        raise RuntimeError(
            f"No child projects could be resolved from the ENA umbrella "
            f"XML for {umbrella}. Last lookup error: {last_error}"
        )

    with open(
        cache_file,
        "w",
        encoding="utf-8"
    ) as handle:
        for project in sorted(projects):
            handle.write(project + "\n")

    print(
        f"Found {len(projects)} child projects under {umbrella} "
        "from ENA umbrella XML"
    )

    return projects


def get_project_name_aegis_samples():
    print()
    print(
        'Finding ENA samples where project_name = "AEGIS" ...'
    )

    rows = query_ena(
        result="sample",
        query='project_name="AEGIS"',
        fields=SAMPLE_FIELDS
    )

    print(
        f"Found {len(rows)} ENA samples with "
        'project_name = "AEGIS"'
    )

    return rows


# ------------------------------------------------------------
# Step 3
# Find samples belonging to umbrella child projects
# ------------------------------------------------------------

def get_umbrella_samples(projects):
    all_rows = []

    projects = sorted(projects)

    print()
    print(
        "Finding samples belonging to umbrella projects ..."
    )

    for project_chunk in chunks(
        projects,
        QUERY_CHUNK_SIZE
    ):
        query = " OR ".join(
            f'study_accession="{project}"'
            for project in project_chunk
        )

        rows = query_ena(
            result="sample",
            query=query,
            fields=SAMPLE_FIELDS
        )

        all_rows.extend(rows)

    print(
        f"Found {len(all_rows)} sample/project relationships "
        "through the umbrella"
    )

    return all_rows


# ------------------------------------------------------------
# Fetch ENA metadata for BioSamples AEGIS samples
# ------------------------------------------------------------

def get_ena_metadata_for_biosamples(accessions):
    accessions = sorted(
        set(accessions)
    )

    if not accessions:
        return []

    print()
    print(
        "Fetching ENA sample metadata for BioSamples AEGIS accessions ..."
    )

    all_rows = []

    for accession_chunk in chunks(
        accessions,
        QUERY_CHUNK_SIZE
    ):

        query = " OR ".join(
            f'accession="{accession}"'
            for accession in accession_chunk
        )

        rows = query_ena(
            result="sample",
            query=query,
            fields=SAMPLE_FIELDS
        )

        all_rows.extend(rows)

    print(
        f"Found ENA metadata for {len(all_rows)} BioSamples-derived samples"
    )

    return all_rows


# ------------------------------------------------------------
# Step 4
# Merge samples from ENA project_name + umbrella + BioSamples
# ------------------------------------------------------------

def merge_samples(
    aegis_samples,
    umbrella_samples,
    biosample_ena_samples,
    biosamples_raw,
    umbrella_projects
):
    merged = {}

    # ENA project_name = AEGIS
    for row in aegis_samples:
        accession = row.get("accession")

        if not accession:
            continue

        row = dict(row)
        row["_project_name_match"] = True
        row["_umbrella_match"] = (
            row.get("study_accession", "")
            in umbrella_projects
        )
        row["_biosamples_match"] = (
            accession in biosamples_raw
        )

        merged[accession] = row

    # Umbrella samples
    for row in umbrella_samples:
        accession = row.get("accession")

        if not accession:
            continue

        if accession not in merged:
            row = dict(row)

            row["_project_name_match"] = (
                row.get(
                    "project_name",
                    ""
                ).strip().upper()
                == "AEGIS"
            )

            row["_umbrella_match"] = True

            row["_biosamples_match"] = (
                accession in biosamples_raw
            )

            merged[accession] = row

        else:
            merged[accession][
                "_umbrella_match"
            ] = True

            if accession in biosamples_raw:
                merged[accession][
                    "_biosamples_match"
                ] = True

    # BioSamples-derived ENA records
    for row in biosample_ena_samples:
        accession = row.get("accession")

        if not accession:
            continue

        if accession not in merged:
            row = dict(row)
            row["_project_name_match"] = (
                row.get(
                    "project_name",
                    ""
                ).strip().upper()
                == "AEGIS"
            )
            row["_umbrella_match"] = (
                row.get(
                    "study_accession",
                    ""
                )
                in umbrella_projects
            )
            row["_biosamples_match"] = True

            merged[accession] = row

        else:
            merged[accession][
                "_biosamples_match"
            ] = True

    # BioSamples records with no ENA sample metadata yet
    for accession, biosample in biosamples_raw.items():

        if accession not in merged:
            merged[accession] = {
                "accession": accession,
                "project_name": "AEGIS",
                "scientific_name": biosample.get(
                    "biosamples_name",
                    ""
                ),
                "biosamples_release": biosample.get(
                    "biosamples_release",
                    ""
                ),
                "biosamples_update": biosample.get(
                    "biosamples_update",
                    ""
                ),
                "_project_name_match": False,
                "_umbrella_match": False,
                "_biosamples_match": True
            }

    # Add BioSamples metadata to all matching samples
    for accession, row in merged.items():

        biosample = biosamples_raw.get(
            accession,
            {}
        )

        if biosample:
            row["biosamples_release"] = biosample.get(
                "biosamples_release",
                ""
            )
            row["biosamples_update"] = biosample.get(
                "biosamples_update",
                ""
            )
            row["biosamples_name"] = biosample.get(
                "biosamples_name",
                ""
            )

        reasons = []

        if row.get(
            "_project_name_match",
            False
        ):
            reasons.append(
                "ENA project_name"
            )

        if row.get(
            "_biosamples_match",
            False
        ):
            reasons.append(
                "BioSamples"
            )

        if row.get(
            "_umbrella_match",
            False
        ):
            reasons.append(
                "umbrella"
            )

        row["included_by"] = (
            " + ".join(reasons)
            if reasons
            else "unknown"
        )

    return merged


# ------------------------------------------------------------
# Step 5
# Retrieve experiments and runs
# ------------------------------------------------------------

RUN_FIELDS = [
    "accession",
    "experiment_accession",
    "sample_accession",
    "study_accession",
    "submission_accession",
    "scientific_name",
    "tax_id",
    "instrument_platform",
    "instrument_model",
    "library_strategy",
    "library_source",
    "library_selection",
    "first_created",
    "first_public",
    "last_updated"
]


def get_runs_for_samples(sample_accessions):
    sample_accessions = sorted(set(sample_accessions))
    all_rows = []

    print()
    print("Finding experiments and runs for all AEGIS samples ...")

    initial_batches = list(
        chunks(sample_accessions, QUERY_CHUNK_SIZE)
    )

    def fetch_batch(batch):
        query = " OR ".join(
            f'sample_accession="{sample}"'
            for sample in batch
        )

        try:
            return query_ena(
                result="read_run",
                query=query,
                fields=RUN_FIELDS
            )

        except requests.RequestException as error:
            if len(batch) == 1:
                print(
                    f"\nWARNING: ENA run lookup failed for "
                    f"{batch[0]} after retries: {error}"
                )
                return []

            midpoint = len(batch) // 2
            left = batch[:midpoint]
            right = batch[midpoint:]

            print(
                f"\nENA batch of {len(batch)} samples still failed; "
                f"splitting into {len(left)} + {len(right)} ..."
            )

            return fetch_batch(left) + fetch_batch(right)

    for number, batch in enumerate(initial_batches, start=1):
        print(
            f"\rQuerying batch {number}/{len(initial_batches)}",
            end="",
            flush=True
        )
        all_rows.extend(fetch_batch(batch))

    print()

    unique = {}
    for row in all_rows:
        accession = row.get("accession", "")
        if accession:
            unique[accession] = row

    rows = list(unique.values())

    print(f"Retrieved {len(rows)} run records")
    return rows


def _first_nonempty(mapping, *keys):
    for key in keys:
        value = mapping.get(key)
        if value not in (None, "", []):
            return value
    return ""


def _extract_biosample_accession(item):
    """
    NCBI Assembly ESummary has used slightly different key spellings over time.
    Accept the common variants and also scan nested values as a fallback.
    """
    direct = _first_nonempty(
        item,
        "biosampleaccn",
        "biosample",
        "BioSampleAccn",
        "BioSample"
    )

    if isinstance(direct, str):
        match = re.search(r"\bSAM(?:E|N|D)[A-Z]?\d+\b", direct)
        if match:
            return match.group(0)

    blob = json.dumps(item)

    match = re.search(
        r"\bSAM(?:E|N|D)[A-Z]?\d+\b",
        blob
    )

    return match.group(0) if match else ""


def _extract_bioproject_accession(item, fallback_project):
    direct = _first_nonempty(
        item,
        "bioprojectaccn",
        "bioproject",
        "BioProjectAccn",
        "BioProject"
    )

    if isinstance(direct, str):
        match = re.search(r"\bPRJ[A-Z]{2}\d+\b", direct)
        if match:
            return match.group(0)

    blob = json.dumps(item)

    match = re.search(
        r"\bPRJ[A-Z]{2}\d+\b",
        blob
    )

    return (
        match.group(0)
        if match
        else fallback_project
    )


def _extract_assembly_accession(item):
    direct = _first_nonempty(
        item,
        "assemblyaccession",
        "assemblyAccession",
        "AssemblyAccession"
    )

    if isinstance(direct, str):
        match = re.search(
            r"\bGC[AF]_\d+(?:\.\d+)?\b",
            direct
        )
        if match:
            return match.group(0)

    blob = json.dumps(item)

    match = re.search(
        r"\bGC[AF]_\d+(?:\.\d+)?\b",
        blob
    )

    return match.group(0) if match else ""


def _extract_assembly_release_date(item):
    value = _first_nonempty(
        item,
        "seqreleasedate",
        "releaseDate",
        "releasedate",
        "SeqReleaseDate"
    )

    return str(value) if value else ""


def _extract_assembly_submission_date(item):
    value = _first_nonempty(
        item,
        "submissiondate",
        "submissionDate",
        "SubmissionDate"
    )

    return str(value) if value else ""



_ncbi_lock = threading.Lock()
_ncbi_last_request_time = 0.0


def ncbi_get(url, params, timeout=60):
    """
    Make a rate-limited NCBI E-utilities GET request.

    Without an API key, NCBI allows up to 3 requests/second.
    This helper spaces EVERY request, includes tool/email identification,
    honors Retry-After when present, and retries HTTP 429/5xx errors with
    exponential backoff.

    Optional:
        export NCBI_API_KEY="..."
    """

    global _ncbi_last_request_time

    params = dict(params)
    params["tool"] = NCBI_TOOL
    params["email"] = NCBI_EMAIL

    if NCBI_API_KEY:
        params["api_key"] = NCBI_API_KEY

    for attempt in range(NCBI_MAX_RETRIES):

        with _ncbi_lock:
            now = time.monotonic()
            elapsed = now - _ncbi_last_request_time

            if elapsed < NCBI_MIN_INTERVAL:
                time.sleep(
                    NCBI_MIN_INTERVAL - elapsed
                )

            _ncbi_last_request_time = time.monotonic()

        response = requests.get(
            url,
            params=params,
            timeout=timeout,
            headers={
                "User-Agent": (
                    "AEGIS-submission-tracker/4.1 "
                    f"({NCBI_EMAIL})"
                )
            }
        )

        if response.status_code == 429:
            retry_after = response.headers.get(
                "Retry-After"
            )

            try:
                wait_seconds = float(
                    retry_after
                )
            except (TypeError, ValueError):
                wait_seconds = min(
                    60,
                    2 ** (attempt + 1)
                )

            print(
                f"\nNCBI rate limit hit (429). "
                f"Retrying in {wait_seconds:.1f}s "
                f"[attempt {attempt + 1}/{NCBI_MAX_RETRIES}] ..."
            )

            time.sleep(
                wait_seconds
            )
            continue

        if response.status_code in {
            500, 502, 503, 504
        }:
            wait_seconds = min(
                60,
                2 ** (attempt + 1)
            )

            print(
                f"\nNCBI temporary error "
                f"{response.status_code}. "
                f"Retrying in {wait_seconds:.1f}s "
                f"[attempt {attempt + 1}/{NCBI_MAX_RETRIES}] ..."
            )

            time.sleep(
                wait_seconds
            )
            continue

        response.raise_for_status()
        return response

    raise requests.HTTPError(
        f"NCBI request failed after {NCBI_MAX_RETRIES} retries: {url}"
    )


ASSEMBLY_FIELDS_DESIRED = [
    "accession",
    "sample_accession",
    "secondary_sample_accession",
    "study_accession",
    "secondary_study_accession",
    "scientific_name",
    "tax_id",
    "first_public",
    "last_updated",
    "description"
]


def get_assemblies_for_samples(
    sample_accessions,
    samples=None,
    umbrella_projects=None
):
    """
    Retrieve ENA public assemblies using the exact query form confirmed
    to work manually:

        result=assembly
        query=sample_accession="SAMEA..."
        fields=accession,description,sample_accession

    Only BioSamples-style accessions (SAMEA/SAMN/SAMD/SAMS) are used here.
    """

    if samples is None:
        samples = {}

    selected_samples = set()

    # Use only BioSamples-style primary sample accessions.
    for accession in sample_accessions:
        if accession:
            value = str(accession).strip()
            if value.startswith(("SAMEA", "SAMN", "SAMD", "SAMS")):
                selected_samples.add(value)

    for key, row in samples.items():
        for value in (
            key,
            row.get("accession", ""),
            row.get("sample_accession", ""),
            row.get("biosample_accession", ""),
        ):
            if value:
                value = str(value).strip()

                if value.startswith(
                    ("SAMEA", "SAMN", "SAMD", "SAMS")
                ):
                    selected_samples.add(value)

    selected_samples = sorted(selected_samples)

    if not selected_samples:
        print()
        print(
            "No BioSamples-style AEGIS sample accessions "
            "available for ENA assembly lookup."
        )
        return []

    print()
    print(
        "Finding ENA assemblies for ALL detected AEGIS BioSamples ..."
    )
    print(
        f"ENA assembly lookup BioSamples: {len(selected_samples)}"
    )

    # These are the exact fields returned by the working manual query.
    fields = [
        "accession",
        "description",
        "sample_accession"
    ]

    all_rows = []

    batches = list(
        chunks(
            selected_samples,
            QUERY_CHUNK_SIZE
        )
    )

    for number, sample_chunk in enumerate(
        batches,
        start=1
    ):

        query = " OR ".join(
            f'sample_accession="{sample}"'
            for sample in sample_chunk
        )

        try:
            rows = query_ena(
                result="assembly",
                query=query,
                fields=fields
            )

            all_rows.extend(rows)

        except requests.HTTPError as error:
            print()
            print(
                f"WARNING: ENA assembly batch "
                f"{number} failed: {error}"
            )

        print(
            f"\rChecked ENA assembly batches "
            f"{number}/{len(batches)}",
            end="",
            flush=True
        )

    print()

    # Deduplicate by public assembly accession (GCA_...).
    unique = {}

    for row in all_rows:
        accession = (
            row.get("accession", "")
            .strip()
        )

        if accession:
            unique[accession] = row

    rows = sorted(
        unique.values(),
        key=lambda row: row.get(
            "accession",
            ""
        )
    )

    print(
        f"Retrieved {len(rows)} ENA assembly records "
        "linked to detected AEGIS BioSamples"
    )

    return rows


def ncbi_datasets_get(url, params=None, timeout=120):
    """
    GET helper for NCBI Datasets with retries for 429/temporary server errors.
    """
    params = dict(params or {})
    api_key = os.getenv("NCBI_API_KEY")

    headers = {
        "Accept": "application/json",
        "User-Agent": f"AEGIS-submission-tracker ({DEFAULT_EMAIL})"
    }

    if api_key:
        headers["api-key"] = api_key

    for attempt in range(6):
        response = requests.get(
            url,
            params=params,
            headers=headers,
            timeout=timeout
        )

        if response.status_code == 404:
            return response

        if response.status_code == 429:
            retry_after = response.headers.get("Retry-After")
            try:
                wait = float(retry_after)
            except (TypeError, ValueError):
                wait = min(60, 2 ** (attempt + 1))

            print(
                f"\nNCBI Datasets rate limit hit. "
                f"Retrying in {wait:.1f}s ..."
            )
            time.sleep(wait)
            continue

        if response.status_code in {500, 502, 503, 504}:
            wait = min(60, 2 ** (attempt + 1))
            print(
                f"\nNCBI Datasets temporary error "
                f"{response.status_code}. Retrying in {wait:.1f}s ..."
            )
            time.sleep(wait)
            continue

        response.raise_for_status()
        return response

    raise requests.HTTPError(
        f"NCBI Datasets request failed after retries: {url}"
    )


def walk_dicts(value):
    if isinstance(value, dict):
        yield value
        for child in value.values():
            yield from walk_dicts(child)
    elif isinstance(value, list):
        for child in value:
            yield from walk_dicts(child)


def first_recursive_value(value, keys):
    keys = set(keys)
    for item in walk_dicts(value):
        for key in keys:
            found = item.get(key)
            if found not in (None, "", []):
                return found
    return ""


def find_recursive_records(value, accession_keys):
    """
    Return nested dictionaries that look like assembly/sequence report rows.
    """
    accession_keys = set(accession_keys)
    found = []

    for item in walk_dicts(value):
        if any(
            item.get(key) not in (None, "", [])
            for key in accession_keys
        ):
            found.append(item)

    return found



def canonical_assembly_accession(accession):
    """
    Normalize GCA/GCF accessions for ENA↔NCBI matching by removing
    the optional version suffix, e.g.:
        GCA_977066645.1 -> GCA_977066645
    """
    if not accession:
        return ""

    return str(accession).strip().split(".", 1)[0]


def get_ncbi_project_statuses(project_accessions):
    """
    Check whether ENA/INSDC BioProject accessions (for example PRJEB...) are
    present/public in NCBI BioProject.

    ENA remains the source of truth for which projects belong to the AEGIS
    tracker. NCBI is used only as an independent public-visibility check.

    The lookup is batched using the BioProject accession search field [PRJA].
    Matching ESummary records are mapped back to their PRJ accession.
    """

    project_accessions = sorted({
        str(accession).strip()
        for accession in project_accessions
        if accession
        and str(accession).strip() != "NO_PROJECT_ACCESSION"
        and re.fullmatch(r"PRJ[A-Z]{2}\d+", str(accession).strip())
    })

    results = {
        accession: {
            "ena_project_status": "Public",
            "ncbi_project_status": "Not found / not public",
            "ncbi_project_accession": "",
            "ncbi_project_uid": ""
        }
        for accession in project_accessions
    }

    if not project_accessions:
        return results

    if args.skip_ncbi:
        for accession in results:
            results[accession]["ncbi_project_status"] = "Not checked"
        return results

    print()
    print("Checking ENA projects in NCBI BioProject ...")

    groups = list(chunks(project_accessions, 20))

    for number, group in enumerate(groups, start=1):
        term = " OR ".join(
            f'"{accession}"[PRJA]'
            for accession in group
        )

        try:
            search_response = ncbi_get(
                NCBI_ESEARCH,
                params={
                    "db": "bioproject",
                    "term": term,
                    "retmode": "json",
                    "retmax": "1000"
                },
                timeout=120
            )

            ids = (
                search_response.json()
                .get("esearchresult", {})
                .get("idlist", [])
            )

            if ids:
                summary_response = ncbi_get(
                    NCBI_ESUMMARY,
                    params={
                        "db": "bioproject",
                        "id": ",".join(str(uid) for uid in ids),
                        "retmode": "json"
                    },
                    timeout=120
                )

                payload = summary_response.json().get("result", {})

                for uid in ids:
                    item = payload.get(str(uid), {})
                    if not isinstance(item, dict):
                        continue

                    blob = json.dumps(item)
                    accessions = re.findall(
                        r"\bPRJ[A-Z]{2}\d+\b",
                        blob
                    )

                    # Prefer an accession that was actually requested in this batch.
                    matched_source = next(
                        (acc for acc in group if acc in accessions),
                        ""
                    )

                    if not matched_source:
                        continue

                    # NCBI may display the same INSDC accession (PRJEB...) rather
                    # than creating a PRJNA alias. Preserve the accession NCBI
                    # actually reports, preferring the source accession when present.
                    reported = (
                        matched_source
                        if matched_source in accessions
                        else (accessions[0] if accessions else "")
                    )

                    results[matched_source] = {
                        "ena_project_status": "Public",
                        "ncbi_project_status": "Public",
                        "ncbi_project_accession": reported,
                        "ncbi_project_uid": str(uid)
                    }

        except requests.RequestException as error:
            print(
                f"\nWARNING: NCBI BioProject status batch "
                f"{number} failed: {error}"
            )

            for accession in group:
                results[accession]["ncbi_project_status"] = "Check failed"

        print(
            f"\rChecked NCBI BioProject batches "
            f"{number}/{len(groups)}",
            end="",
            flush=True
        )

        time.sleep(0.3)

    print()
    return results


def build_project_status_table(project_accessions, ncbi_project_status):
    """Build one row per ENA project for the project public-status report."""
    rows = []

    for accession in sorted(set(project_accessions)):
        if not accession or accession == "NO_PROJECT_ACCESSION":
            continue

        status = ncbi_project_status.get(accession, {})

        rows.append({
            "project_accession": accession,
            "ena_project_status": status.get(
                "ena_project_status", "Public"
            ),
            "ncbi_project_status": status.get(
                "ncbi_project_status", "Not checked"
            ),
            "ncbi_project_accession": status.get(
                "ncbi_project_accession", ""
            ),
            "ncbi_project_uid": status.get(
                "ncbi_project_uid", ""
            )
        })

    return rows


def get_ncbi_assembly_statuses(assembly_accessions):
    """
    Check whether each ENA GCA/GCF assembly is public in NCBI Datasets.

    ENA's assembly result may return an unversioned accession such as
    GCA_977066645, while NCBI reports GCA_977066645.1. Matching is therefore
    done on the accession with the version suffix removed.
    """

    input_accessions = sorted(set(
        accession for accession in assembly_accessions if accession
    ))

    results = {
        accession: {
            "ncbi_status": "Not found / not public",
            "ncbi_accession": "",
            "ncbi_release_date": "",
            "ncbi_assembly_status": "",
            "ncbi_assembly_level": ""
        }
        for accession in input_accessions
    }

    if not input_accessions or args.skip_ncbi:
        if args.skip_ncbi:
            for accession in results:
                results[accession]["ncbi_status"] = "Not checked"
        return results

    # base accession -> ENA accession(s)
    canonical_to_input = defaultdict(list)

    for accession in input_accessions:
        canonical_to_input[
            canonical_assembly_accession(accession)
        ].append(accession)

    print()
    print("Checking ENA assemblies in NCBI Datasets ...")

    groups = list(chunks(input_accessions, 40))

    for number, group in enumerate(groups, start=1):

        joined = ",".join(group)

        url = (
            f"{NCBI_DATASETS_BASE}/genome/accession/"
            f"{joined}/dataset_report"
        )

        try:
            response = ncbi_datasets_get(
                url,
                params={"page_size": 1000},
                timeout=120
            )

            if response.status_code != 404:
                payload = response.json()

                # Current NCBI Datasets response uses payload["reports"].
                report_records = payload.get("reports", [])

                # Fallback for any schema/wrapper variation.
                if not isinstance(report_records, list):
                    report_records = []

                if not report_records:
                    report_records = find_recursive_records(
                        payload,
                        {
                            "accession",
                            "assembly_accession",
                            "assemblyAccession"
                        }
                    )

                for record in report_records:
                    reported_accession = str(
                        record.get("accession")
                        or record.get("assembly_accession")
                        or record.get("assemblyAccession")
                        or first_recursive_value(
                            record,
                            {
                                "accession",
                                "assembly_accession",
                                "assemblyAccession"
                            }
                        )
                        or ""
                    ).strip()

                    if not reported_accession:
                        continue

                    base = canonical_assembly_accession(
                        reported_accession
                    )

                    matching_inputs = canonical_to_input.get(
                        base,
                        []
                    )

                    if not matching_inputs:
                        continue

                    # Current NCBI assembly report is nested under assembly_info.
                    assembly_info = (
                        record.get("assembly_info")
                        or record.get("assemblyInfo")
                        or {}
                    )

                    release_date = str(
                        first_recursive_value(
                            record,
                            {
                                "release_date",
                                "releaseDate",
                                "seq_release_date",
                                "seqReleaseDate"
                            }
                        )
                        or ""
                    )

                    assembly_status = str(
                        record.get("assembly_status")
                        or record.get("assemblyStatus")
                        or first_recursive_value(
                            record,
                            {
                                "assembly_status",
                                "assemblyStatus"
                            }
                        )
                        or ""
                    )

                    assembly_level = str(
                        assembly_info.get("assembly_level")
                        or assembly_info.get("assemblyLevel")
                        or first_recursive_value(
                            record,
                            {
                                "assembly_level",
                                "assemblyLevel"
                            }
                        )
                        or ""
                    )

                    for original_accession in matching_inputs:
                        results[original_accession] = {
                            "ncbi_status": "Public",
                            "ncbi_accession": reported_accession,
                            "ncbi_release_date": release_date,
                            "ncbi_assembly_status": assembly_status,
                            "ncbi_assembly_level": assembly_level
                        }

        except requests.RequestException as error:
            print(
                f"\nWARNING: NCBI assembly status batch "
                f"{number} failed: {error}"
            )

            for accession in group:
                results[accession]["ncbi_status"] = "Check failed"

        print(
            f"\rChecked NCBI assembly batches "
            f"{number}/{len(groups)}",
            end="",
            flush=True
        )

        time.sleep(0.2)

    print()

    return results

def classify_assembly_component(record):
    location = str(
        record.get("assignedMoleculeLocationType", "")
        or record.get("assigned_molecule_location_type", "")
    ).strip()

    chromosome = str(
        record.get("chrName", "")
        or record.get("chr_name", "")
    ).strip()

    role = str(
        record.get("role", "")
    ).strip()

    if (
        location.casefold() == "chromosome"
        or (
            chromosome
            and chromosome.casefold() not in {"un", "unknown"}
            and role.casefold() == "assembled-molecule"
        )
    ):
        return "chromosome"

    if "scaffold" in role.casefold():
        return "scaffold"

    if "contig" in role.casefold():
        return "contig"

    if role.casefold() in {
        "unplaced-scaffold",
        "unlocalized-scaffold"
    }:
        return "scaffold"

    return "sequence/contig"


def get_ena_assembly_components_metadata(
    assembly_rows,
    samples
):
    """
    Fast ENA component retrieval using metadata only.

    No assembly FASTA files are downloaded.

    If a BioSample has exactly one detected assembly, component rows are
    mapped directly to that GCA accession. If it has multiple assemblies,
    all candidate GCA accessions are retained and mapping_status is marked
    ambiguous rather than attempting slow FASTA-based disambiguation.
    """

    if not assembly_rows:
        return []

    print()
    print(
        "Fetching assembly component metadata from ENA "
        "(fast metadata-only mode; no FASTA downloads) ..."
    )

    assemblies_by_sample = defaultdict(list)

    for assembly in assembly_rows:
        sample_accession = str(
            assembly.get("sample_accession", "")
        ).strip()

        if sample_accession:
            assemblies_by_sample[sample_accession].append(assembly)

    sample_accessions = sorted(assemblies_by_sample)

    if not sample_accessions:
        return []

    supported = get_ena_return_fields("sequence")

    desired_fields = [
        "accession",
        "description",
        "sample_accession",
        "study_accession",
        "scientific_name",
        "tax_id",
        "base_count",
        "sequence_length",
        "first_public",
        "last_updated"
    ]

    if supported:
        fields = [
            field
            for field in desired_fields
            if field in supported
        ]
    else:
        fields = [
            "accession",
            "description",
            "sample_accession"
        ]

    if "accession" not in fields:
        fields.insert(0, "accession")

    print(
        "ENA sequence metadata fields: "
        + ", ".join(fields)
    )

    components = []
    seen = set()

    batches = list(
        chunks(
            sample_accessions,
            QUERY_CHUNK_SIZE
        )
    )

    for number, sample_chunk in enumerate(
        batches,
        start=1
    ):
        query = " OR ".join(
            f'sample_accession="{sample}"'
            for sample in sample_chunk
        )

        try:
            rows = query_ena(
                result="sequence",
                query=query,
                fields=fields
            )
        except requests.RequestException as error:
            print(
                f"\nWARNING: ENA sequence metadata batch "
                f"{number} failed: {error}"
            )
            rows = []

        for record in rows:
            component_accession = str(
                record.get("accession", "")
            ).strip()

            sample_accession = str(
                record.get("sample_accession", "")
            ).strip()

            if not component_accession or not sample_accession:
                continue

            candidate_assemblies = assemblies_by_sample.get(
                sample_accession,
                []
            )

            candidate_accessions = sorted({
                row.get("assembly_accession", "")
                for row in candidate_assemblies
                if row.get("assembly_accession")
            })

            if len(candidate_accessions) == 1:
                assembly_accession = candidate_accessions[0]
                mapping_status = "exact_single_assembly_for_sample"

            elif len(candidate_accessions) > 1:
                assembly_accession = ";".join(candidate_accessions)
                mapping_status = "ambiguous_multiple_assemblies_for_sample"

            else:
                assembly_accession = ""
                mapping_status = "sample_link_only"

            key = (
                component_accession,
                sample_accession,
                assembly_accession
            )

            if key in seen:
                continue

            seen.add(key)

            description = str(
                record.get("description", "")
            )

            lower = description.casefold()

            if "mitochond" in lower:
                component_type = "mitochondrion"

            elif "chloroplast" in lower or "plastid" in lower:
                component_type = "chloroplast"

            elif "chromosome" in lower:
                component_type = "chromosome"

            elif "scaffold" in lower:
                component_type = "scaffold"

            elif "contig" in lower:
                component_type = "contig"

            else:
                component_type = "assembled sequence"

            chromosome_name = ""

            chromosome_match = re.search(
                r"\bchromosome\s+([^\s,;]+)",
                description,
                flags=re.IGNORECASE
            )

            if chromosome_match:
                chromosome_name = chromosome_match.group(1)

            components.append({
                "assembly_accession":
                    assembly_accession,

                "sample_accession":
                    sample_accession,

                "project_accession":
                    (
                        candidate_assemblies[0].get(
                            "project_accession",
                            ""
                        )
                        if candidate_assemblies
                        else record.get(
                            "study_accession",
                            ""
                        )
                    ),

                "project_name":
                    (
                        candidate_assemblies[0].get(
                            "project_name",
                            ""
                        )
                        if candidate_assemblies
                        else samples.get(
                            sample_accession,
                            {}
                        ).get(
                            "project_name",
                            ""
                        )
                    ),

                "component_accession":
                    component_accession,

                "genbank_accession":
                    component_accession,

                "refseq_accession":
                    "",

                "component_type":
                    component_type,

                "chromosome_name":
                    chromosome_name,

                "sequence_name":
                    "",

                "role":
                    "",

                "molecule_location_type":
                    "",

                "length":
                    (
                        record.get(
                            "sequence_length",
                            ""
                        )
                        or record.get(
                            "base_count",
                            ""
                        )
                    ),

                "description":
                    description,

                "first_public":
                    record.get(
                        "first_public",
                        ""
                    ),

                "last_updated":
                    record.get(
                        "last_updated",
                        ""
                    ),

                "mapping_status":
                    mapping_status,

                "ena_status":
                    "Public",

                "ncbi_status":
                    "Not checked (assembly only)"
            })

        print(
            f"\rChecked ENA sequence metadata batches "
            f"{number}/{len(batches)} "
            f"(components so far: {len(components)})",
            end="",
            flush=True
        )

    print()

    ambiguous = sum(
        1
        for row in components
        if row.get("mapping_status")
        == "ambiguous_multiple_assemblies_for_sample"
    )

    print(
        f"Retrieved {len(components)} ENA assembly component records"
    )

    print(
        f"Ambiguous component-to-GCA mappings: {ambiguous}"
    )

    print(
        "No assembly FASTA files were downloaded."
    )

    return components


def add_ena_component_status(component_rows):
    """
    Query ENA's public sequence result for GenBank/INSDC component accessions.
    If the accession is returned, it is public in ENA.
    """
    accession_to_rows = defaultdict(list)

    for row in component_rows:
        accession = row.get("genbank_accession", "") or row.get(
            "component_accession", ""
        )

        if accession:
            # ENA sequence searches normally use the base accession without
            # an NCBI RefSeq-only accession if a GenBank accession is present.
            accession_to_rows[accession].append(row)

    accessions = sorted(accession_to_rows)

    if not accessions:
        return component_rows

    print()
    print("Checking chromosome/contig accessions in ENA ...")

    public_in_ena = set()
    groups = list(chunks(accessions, 40))

    for number, group in enumerate(groups, start=1):
        query = " OR ".join(
            f'accession="{accession}"'
            for accession in group
        )

        try:
            returned = query_ena(
                result="sequence",
                query=query,
                fields=["accession", "last_updated"]
            )

            for record in returned:
                accession = record.get("accession", "")
                if accession:
                    public_in_ena.add(accession)

        except requests.RequestException as error:
            print(
                f"\nWARNING: ENA component status batch "
                f"{number} failed: {error}"
            )

        print(
            f"\rChecked ENA component batches "
            f"{number}/{len(groups)}",
            end="",
            flush=True
        )

    print()

    for accession, rows in accession_to_rows.items():
        status = (
            "Public"
            if accession in public_in_ena
            else "Not found / not public"
        )
        for row in rows:
            row["ena_status"] = status

    return component_rows


def build_assembly_table(
    assemblies,
    samples,
    ncbi_assembly_status=None
):
    """
    Build the assembly TSV from ENA result=assembly records.
    """

    rows = []

    if ncbi_assembly_status is None:
        ncbi_assembly_status = {}

    sample_lookup = {}

    for key, sample in samples.items():
        for identifier in (
            key,
            sample.get("accession", ""),
            sample.get("sample_accession", ""),
            sample.get("biosample_accession", ""),
        ):
            if identifier:
                sample_lookup[
                    str(identifier).strip()
                ] = sample

    for assembly in assemblies:

        sample_accession = (
            assembly.get(
                "sample_accession",
                ""
            )
            .strip()
        )

        sample = sample_lookup.get(
            sample_accession,
            {}
        )

        assembly_accession = assembly.get(
            "accession",
            ""
        )

        ncbi = ncbi_assembly_status.get(
            assembly_accession,
            {
                "ncbi_status": "Not checked",
                "ncbi_release_date": "",
                "ncbi_assembly_status": "",
                "ncbi_assembly_level": ""
            }
        )

        rows.append({
            "assembly_accession":
                assembly_accession,

            "description":
                assembly.get(
                    "description",
                    ""
                ),

            "sample_accession":
                sample_accession,

            "project_accession":
                sample.get(
                    "study_accession",
                    ""
                ),

            "project_name":
                sample.get(
                    "project_name",
                    ""
                ),

            "scientific_name":
                sample.get(
                    "scientific_name",
                    ""
                ),

            # result=assembly returns public ENA assembly records.
            "ena_status":
                "Public",

            "ncbi_status":
                ncbi.get(
                    "ncbi_status",
                    "Not checked"
                ),

            "ncbi_release_date":
                ncbi.get(
                    "ncbi_release_date",
                    ""
                ),

            "ncbi_assembly_status":
                ncbi.get(
                    "ncbi_assembly_status",
                    ""
                ),

            "ncbi_assembly_level":
                ncbi.get(
                    "ncbi_assembly_level",
                    ""
                ),

            "included_by":
                sample.get(
                    "included_by",
                    ""
                )
        })

    return rows

# ------------------------------------------------------------
# Step 6
# NCBI public visibility
# ------------------------------------------------------------

def ncbi_run_public(run_accession):
    try:
        response = ncbi_get(
            NCBI_ESEARCH,
            params={
                "db": "sra",
                "term": f'"{run_accession}"[Accession]',
                "retmode": "json",
                "retmax": "1"
            },
            timeout=60
        )

        response.raise_for_status()

        count = int(
            response.json()
            .get("esearchresult", {})
            .get("count", "0")
        )

        return count > 0

    except Exception:
        return None


def check_ncbi_runs(run_accessions, history=None):
    """
    Check ENA run accessions in NCBI using the SRA Data Locator v2 API.

    This avoids Entrez E-utilities, which is heavily rate-limited on shared
    HPC IP addresses.

    The SDL endpoint accepts multiple repeated `acc` parameters and returns
    a status per accession. With meta-only=yes no sequence files are
    downloaded.

    Status interpretation:
        200 -> public/resolvable in NCBI SRA
        403 -> protected / access denied (not public)
        404 -> accession cannot be resolved (not found / not public)

    Runs already confirmed public in history are not queried again.
    """

    run_accessions = sorted({
        accession
        for accession in run_accessions
        if accession
    })

    if args.skip_ncbi or args.skip_ncbi_runs:
        if args.skip_ncbi_runs and not args.skip_ncbi:
            print()
            print(
                "Skipping NCBI SRA run checks (--skip-ncbi-runs); "
                "NCBI assembly checks remain enabled."
            )

        return {
            accession: None
            for accession in run_accessions
        }

    if history is None:
        history = {}

    historical_runs = history.get(
        "runs",
        {}
    )

    results = {}
    to_check = []

    for accession in run_accessions:
        previous = historical_runs.get(
            accession,
            {}
        )

        if previous.get(
            "ncbi_first_seen_public"
        ):
            results[accession] = True
        else:
            to_check.append(accession)

    cached_count = len(run_accessions) - len(to_check)

    print()
    print(
        "Checking public visibility in NCBI SRA "
        "(SRA Data Locator, cached + batched) ..."
    )
    print(
        f"Runs already confirmed public from history: "
        f"{cached_count}"
    )
    print(
        f"Runs requiring an NCBI lookup this run: "
        f"{len(to_check)}"
    )

    if not to_check:
        return results

    SDL_URL = (
        "https://locate.ncbi.nlm.nih.gov/"
        "sdl/2/retrieve"
    )

    # 100 accessions per request keeps the URL/request size reasonable.
    groups = list(
        chunks(
            to_check,
            100
        )
    )

    for number, group in enumerate(
        groups,
        start=1
    ):
        # Unknown until the request succeeds.
        for accession in group:
            results[accession] = None

        params = [
            ("acc", accession)
            for accession in group
        ]
        params.append(
            ("meta-only", "yes")
        )

        try:
            response = requests.get(
                SDL_URL,
                params=params,
                timeout=180,
                headers={
                    "Accept": "application/json",
                    "User-Agent": (
                        "AEGIS-submission-tracker/7.0 "
                        f"({DEFAULT_EMAIL})"
                    )
                }
            )

            response.raise_for_status()
            payload = response.json()

            # SDL may return:
            #   {"version":"2","result":[...]}
            # or a direct list of per-accession status objects.
            if isinstance(payload, dict):
                records = payload.get(
                    "result",
                    []
                )
            elif isinstance(payload, list):
                records = payload
            else:
                records = []

            returned = set()

            for record in records:
                if not isinstance(record, dict):
                    continue

                accession = str(
                    record.get("bundle")
                    or record.get("accession")
                    or ""
                ).strip()

                if not accession:
                    continue

                returned.add(accession)

                try:
                    status = int(
                        record.get(
                            "status",
                            0
                        )
                    )
                except (
                    TypeError,
                    ValueError
                ):
                    status = 0

                if status == 200:
                    results[accession] = True

                elif status in {
                    403,
                    404
                }:
                    results[accession] = False

                else:
                    results[accession] = None

            # If SDL omitted an accession entirely, leave it unresolved.
            for accession in group:
                if accession not in returned:
                    results[accession] = None

        except (
            requests.RequestException,
            ValueError
        ) as error:
            print(
                f"\nWARNING: NCBI SDL run batch "
                f"{number} failed: {error}"
            )

        print(
            f"\rChecked NCBI SDL run batches "
            f"{number}/{len(groups)}",
            end="",
            flush=True
        )

        # Small courtesy delay; SDL is not Entrez E-utilities.
        time.sleep(0.1)

    print()

    public_count = sum(
        value is True
        for value in results.values()
    )

    not_public_count = sum(
        value is False
        for value in results.values()
    )

    unknown_count = sum(
        value is None
        for value in results.values()
    )

    print(
        f"NCBI-public runs after cache/SDL checks: "
        f"{public_count}/{len(run_accessions)}"
    )

    if not_public_count:
        print(
            f"NCBI runs not found/public: "
            f"{not_public_count}"
        )

    if unknown_count:
        print(
            f"NCBI run checks unresolved: "
            f"{unknown_count}"
        )

    return results


def load_history(filename):
    if not os.path.exists(filename):
        return {
            "runs": {},
            "snapshots": []
        }

    with open(
        filename,
        "r",
        encoding="utf-8"
    ) as handle:
        return json.load(handle)


def save_history(filename, history):
    with open(
        filename,
        "w",
        encoding="utf-8"
    ) as handle:
        json.dump(
            history,
            handle,
            indent=2,
            sort_keys=True
        )


# ------------------------------------------------------------
# Step 8
# Detailed accession table
# ------------------------------------------------------------

def build_record_table(
    samples,
    runs,
    ncbi_status,
    history
):
    records = []

    runs_by_sample = defaultdict(list)

    for run in runs:
        sample = run.get(
            "sample_accession",
            ""
        )

        runs_by_sample[sample].append(run)

    today = date.today().isoformat()

    for sample_accession, sample in samples.items():

        sample_runs = runs_by_sample.get(
            sample_accession,
            []
        )

        # Sample with no run yet
        if not sample_runs:
            ena_public_date = sample.get(
                "first_public",
                ""
            )

            records.append({
                "project_accession":
                    sample.get(
                        "study_accession",
                        ""
                    ),

                "sample_accession":
                    sample_accession,

                "experiment_accession":
                    "",

                "run_accession":
                    "",

                "submission_accession":
                    sample.get(
                        "submission_accession",
                        ""
                    ),

                "project_name":
                    sample.get(
                        "project_name",
                        ""
                    ),

                "scientific_name":
                    sample.get(
                        "scientific_name",
                        ""
                    ),

                "tax_id":
                    sample.get(
                        "tax_id",
                        ""
                    ),

                "first_created":
                    "",

                "ena_status":
                    (
                        "Public"
                        if ena_public_date
                        else "Not public / unknown"
                    ),

                "ena_public_date":
                    ena_public_date,

                "ena_days_since_public":
                    days_since(
                        ena_public_date
                    ),

                "ncbi_status":
                    "No run to check",

                "ncbi_first_seen_public":
                    "",

                "ncbi_days_since_first_seen":
                    "",

                "biosamples_release":
                    sample.get(
                        "biosamples_release",
                        ""
                    ),

                "biosamples_update":
                    sample.get(
                        "biosamples_update",
                        ""
                    ),

                "last_updated":
                    sample.get(
                        "last_updated",
                        ""
                    ),

                "included_by":
                    sample.get(
                        "included_by",
                        ""
                    )
            })

            continue

        # Sample has one or more runs
        for run in sample_runs:

            run_accession = run.get(
                "accession",
                ""
            )

            ena_public_date = (
                run.get(
                    "first_public",
                    ""
                )
                or sample.get(
                    "first_public",
                    ""
                )
            )

            if ena_public_date:
                ena_status = "Public"
            else:
                ena_status = (
                    "Not public / unknown"
                )

            ncbi_value = ncbi_status.get(
                run_accession
            )

            if ncbi_value is True:
                ncbi_label = "Public"

            elif ncbi_value is False:
                ncbi_label = (
                    "Not found / not public"
                )

            else:
                ncbi_label = (
                    "Not checked / check failed"
                )

            run_history = history.setdefault(
                "runs",
                {}
            ).setdefault(
                run_accession,
                {}
            )

            if (
                ena_public_date
                and not run_history.get(
                    "ena_first_public"
                )
            ):
                run_history[
                    "ena_first_public"
                ] = ena_public_date

            if (
                ncbi_value is True
                and not run_history.get(
                    "ncbi_first_seen_public"
                )
            ):
                run_history[
                    "ncbi_first_seen_public"
                ] = today

            ncbi_first_seen = (
                run_history.get(
                    "ncbi_first_seen_public",
                    ""
                )
            )

            records.append({
                "project_accession":
                    run.get(
                        "study_accession",
                        ""
                    )
                    or sample.get(
                        "study_accession",
                        ""
                    ),

                "sample_accession":
                    sample_accession,

                "experiment_accession":
                    run.get(
                        "experiment_accession",
                        ""
                    ),

                "run_accession":
                    run_accession,

                "submission_accession":
                    run.get(
                        "submission_accession",
                        ""
                    )
                    or sample.get(
                        "submission_accession",
                        ""
                    ),

                "project_name":
                    sample.get(
                        "project_name",
                        ""
                    ),

                "scientific_name":
                    run.get(
                        "scientific_name",
                        ""
                    )
                    or sample.get(
                        "scientific_name",
                        ""
                    ),

                "tax_id":
                    run.get(
                        "tax_id",
                        ""
                    )
                    or sample.get(
                        "tax_id",
                        ""
                    ),

                "first_created":
                    run.get(
                        "first_created",
                        ""
                    ),

                "ena_status":
                    ena_status,

                "ena_public_date":
                    ena_public_date,

                "ena_days_since_public":
                    days_since(
                        ena_public_date
                    ),

                "ncbi_status":
                    ncbi_label,

                "ncbi_first_seen_public":
                    ncbi_first_seen,

                "ncbi_days_since_first_seen":
                    days_since(
                        ncbi_first_seen
                    ),

                "biosamples_release":
                    sample.get(
                        "biosamples_release",
                        ""
                    ),

                "biosamples_update":
                    sample.get(
                        "biosamples_update",
                        ""
                    ),

                "last_updated":
                    run.get(
                        "last_updated",
                        ""
                    ),

                "included_by":
                    sample.get(
                        "included_by",
                        ""
                    )
            })

    return records


# ------------------------------------------------------------
# Step 9
# Project summary
# ------------------------------------------------------------

def build_project_summary(records):

    projects = defaultdict(
        lambda: {
            "samples": set(),
            "experiments": set(),
            "runs": set(),
            "submission_accessions": set(),
            "project_names": set(),
            "included_by": set(),
            "ena_public_runs": set(),
            "ena_not_public_runs": set(),
            "ncbi_public_runs": set(),
            "ncbi_not_public_runs": set()
        }
    )

    for row in records:

        project = row.get(
            "project_accession",
            ""
        )

        if not project:
            project = "NO_PROJECT_ACCESSION"

        p = projects[project]

        if row.get("sample_accession"):
            p["samples"].add(
                row["sample_accession"]
            )

        if row.get(
            "experiment_accession"
        ):
            p["experiments"].add(
                row[
                    "experiment_accession"
                ]
            )

        run_accession = row.get(
            "run_accession"
        )

        if run_accession:
            p["runs"].add(
                run_accession
            )

            if (
                row.get("ena_status")
                == "Public"
            ):
                p["ena_public_runs"].add(
                    run_accession
                )
            else:
                p["ena_not_public_runs"].add(
                    run_accession
                )

            if (
                row.get("ncbi_status")
                == "Public"
            ):
                p["ncbi_public_runs"].add(
                    run_accession
                )
            else:
                p["ncbi_not_public_runs"].add(
                    run_accession
                )

        if row.get(
            "submission_accession"
        ):
            p[
                "submission_accessions"
            ].add(
                row[
                    "submission_accession"
                ]
            )

        if row.get("project_name"):
            p["project_names"].add(
                row["project_name"]
            )

        if row.get("included_by"):
            p["included_by"].add(
                row["included_by"]
            )

    output = []

    for project, data in projects.items():

        output.append({
            "project_accession":
                project,

            "project_name":
                "; ".join(
                    sorted(
                        data[
                            "project_names"
                        ]
                    )
                ),

            "samples":
                len(
                    data["samples"]
                ),

            "experiments":
                len(
                    data[
                        "experiments"
                    ]
                ),

            "runs":
                len(
                    data["runs"]
                ),

            "submissions":
                len(
                    data[
                        "submission_accessions"
                    ]
                ),

            "ena_public_runs":
                len(
                    data[
                        "ena_public_runs"
                    ]
                ),

            "ena_not_public_or_unknown_runs":
                len(
                    data[
                        "ena_not_public_runs"
                    ]
                ),

            "ncbi_public_runs":
                len(
                    data[
                        "ncbi_public_runs"
                    ]
                ),

            "ncbi_not_found_or_not_public_runs":
                len(
                    data[
                        "ncbi_not_public_runs"
                    ]
                ),

            "included_by":
                "; ".join(
                    sorted(
                        data[
                            "included_by"
                        ]
                    )
                )
        })

    return sorted(
        output,
        key=lambda x: x[
            "project_accession"
        ]
    )


# ------------------------------------------------------------
# Step 10
# Write TSV
# ------------------------------------------------------------

def write_tsv(
    filename,
    rows,
    fieldnames
):

    with open(
        filename,
        "w",
        newline="",
        encoding="utf-8"
    ) as handle:

        writer = csv.DictWriter(
            handle,
            fieldnames=fieldnames,
            delimiter="\t",
            extrasaction="ignore"
        )

        writer.writeheader()
        writer.writerows(rows)


# ------------------------------------------------------------
# Step 11
# Console table and summary
# ------------------------------------------------------------

def print_project_table(projects):

    print()
    print("PROJECT SUMMARY")
    print("=" * 145)

    header = (
        f"{'Project':<15}"
        f"{'Project name':<28}"
        f"{'Samples':>10}"
        f"{'Experiments':>13}"
        f"{'Runs':>10}"
        f"{'Submissions':>13}"
        f"{'ENA public':>12}"
        f"{'NCBI public':>13}"
        f"  {'Included by'}"
    )

    print(header)
    print("-" * 145)

    for row in projects:

        project_name = row[
            "project_name"
        ][:26]

        print(
            f"{row['project_accession']:<15}"
            f"{project_name:<28}"
            f"{row['samples']:>10}"
            f"{row['experiments']:>13}"
            f"{row['runs']:>10}"
            f"{row['submissions']:>13}"
            f"{row['ena_public_runs']:>12}"
            f"{row['ncbi_public_runs']:>13}"
            f"  {row['included_by']}"
        )


def make_snapshot(
    samples,
    runs,
    project_summary,
    umbrella_projects,
    records,
    biosamples_count
):

    return {
        "date":
            date.today().isoformat(),

        "umbrella_child_projects":
            len(
                umbrella_projects
            ),

        "projects":
            len(
                project_summary
            ),

        "samples":
            len(
                samples
            ),

        "biosamples_aegis":
            biosamples_count,

        "experiments":
            len({
                row.get(
                    "experiment_accession"
                )
                for row in runs
                if row.get(
                    "experiment_accession"
                )
            }),

        "runs":
            len({
                row.get("accession")
                for row in runs
                if row.get("accession")
            }),

        "submissions":
            len({
                row.get(
                    "submission_accession"
                )
                for row in runs
                if row.get(
                    "submission_accession"
                )
            }),

        "ena_public_runs":
            len({
                row["run_accession"]
                for row in records
                if (
                    row.get(
                        "run_accession"
                    )
                    and row.get(
                        "ena_status"
                    ) == "Public"
                )
            }),

        "ncbi_public_runs":
            len({
                row["run_accession"]
                for row in records
                if (
                    row.get(
                        "run_accession"
                    )
                    and row.get(
                        "ncbi_status"
                    ) == "Public"
                )
            })
    }


def change_text(
    current,
    previous,
    key
):

    value = current[key]

    if not previous:
        return str(value)

    old = previous.get(
        key,
        0
    )

    change = value - old

    return (
        f"{value} ({change:+d})"
    )


def make_summary(
    samples,
    runs,
    assembly_rows,
    assembly_component_rows,
    publication_issues,
    project_status_rows,
    project_summary,
    umbrella_projects,
    records,
    current,
    previous
):

    ena_project_only = 0
    biosamples_only = 0
    umbrella_only = 0
    ena_plus_biosamples = 0
    biosamples_plus_umbrella = 0
    ena_plus_umbrella = 0
    all_three = 0

    for row in samples.values():

        ena = row.get(
            "_project_name_match",
            False
        )

        bio = row.get(
            "_biosamples_match",
            False
        )

        umbrella = row.get(
            "_umbrella_match",
            False
        )

        if ena and bio and umbrella:
            all_three += 1

        elif ena and bio:
            ena_plus_biosamples += 1

        elif bio and umbrella:
            biosamples_plus_umbrella += 1

        elif ena and umbrella:
            ena_plus_umbrella += 1

        elif ena:
            ena_project_only += 1

        elif bio:
            biosamples_only += 1

        elif umbrella:
            umbrella_only += 1

    samples_with_runs = {
        row.get(
            "sample_accession"
        )
        for row in runs
        if row.get(
            "sample_accession"
        )
    }

    samples_without_runs = (
        set(samples.keys())
        - samples_with_runs
    )

    ena_public = {
        row["run_accession"]
        for row in records
        if (
            row.get(
                "run_accession"
            )
            and row.get(
                "ena_status"
            ) == "Public"
        )
    }

    ncbi_public = {
        row["run_accession"]
        for row in records
        if (
            row.get(
                "run_accession"
            )
            and row.get(
                "ncbi_status"
            ) == "Public"
        )
    }

    ena_only = (
        ena_public
        - ncbi_public
    )

    both = (
        ena_public
        & ncbi_public
    )

    ena_public_projects = {
        row.get("project_accession")
        for row in project_status_rows
        if row.get("ena_project_status") == "Public"
    }

    ncbi_public_projects = {
        row.get("project_accession")
        for row in project_status_rows
        if row.get("ncbi_project_status") == "Public"
    }

    projects_both = ena_public_projects & ncbi_public_projects
    projects_ena_only = ena_public_projects - ncbi_public_projects
    projects_check_failed = {
        row.get("project_accession")
        for row in project_status_rows
        if row.get("ncbi_project_status") == "Check failed"
    }

    summary = (
        "\n"
        "AEGIS WEEKLY SUBMISSION TRACKING\n"
        "=================================================================\n"
        f"Generated:               {datetime.now():%Y-%m-%d %H:%M}\n"
        f"Umbrella:                {args.umbrella}\n"
        "\n"
        "TOTALS\n"
        "-----------------------------------------------------------------\n"
        f"Umbrella child projects: {change_text(current, previous, 'umbrella_child_projects')}\n"
        f"Projects with records:   {change_text(current, previous, 'projects')}\n"
        f"Samples:                 {change_text(current, previous, 'samples')}\n"
        f"BioSamples AEGIS:        {change_text(current, previous, 'biosamples_aegis')}\n"
        f"Experiments:             {change_text(current, previous, 'experiments')}\n"
        f"Runs:                    {change_text(current, previous, 'runs')}\n"
        f"Assemblies:              {len(assembly_rows)}\n"
        f"ENA public assemblies:   {sum(1 for row in assembly_rows if row.get('ena_status') == 'Public')}\n"
        f"NCBI public assemblies:  {sum(1 for row in assembly_rows if row.get('ncbi_status') == 'Public')}\n"
        f"Public assemblies both:  {sum(1 for row in assembly_rows if row.get('ena_status') == 'Public' and row.get('ncbi_status') == 'Public')}\n"
        f"Assembly components:     {len(assembly_component_rows)}\n"
        f"Chromosomes:             {sum(1 for row in assembly_component_rows if row.get('component_type') == 'chromosome')}\n"
        f"Contig/scaffold seqs:    {sum(1 for row in assembly_component_rows if row.get('component_type') != 'chromosome')}\n"
        f"Components public ENA:   {sum(1 for row in assembly_component_rows if row.get('ena_status') == 'Public')}\n"
        f"Components public NCBI:  Not checked (NCBI assembly-only)\n"
        f"Potential delays >={args.stale_days}d: {len(publication_issues)}\n"
        f"Submission accessions:   {change_text(current, previous, 'submissions')}\n"
        f"Samples without runs:    {len(samples_without_runs)}\n"
        "\n"
        "FOUND BY\n"
        "-----------------------------------------------------------------\n"
        f"ENA project_name only:       {ena_project_only}\n"
        f"BioSamples only:             {biosamples_only}\n"
        f"Umbrella only:               {umbrella_only}\n"
        f"ENA + BioSamples:            {ena_plus_biosamples}\n"
        f"ENA + umbrella:              {ena_plus_umbrella}\n"
        f"BioSamples + umbrella:       {biosamples_plus_umbrella}\n"
        f"All three sources:           {all_three}\n"
        "\n"
        "PROJECT PUBLIC STATUS\n"
        "-----------------------------------------------------------------\n"
        f"ENA projects checked:     {len(ena_public_projects)}\n"
        f"NCBI public projects:     {len(ncbi_public_projects)}\n"
        f"Public projects in both:  {len(projects_both)}\n"
        f"ENA public, NCBI absent:  {len(projects_ena_only)}\n"
        f"NCBI project check failed:{len(projects_check_failed)}\n"
        "\n"
        "RUN PUBLIC STATUS\n"
        "-----------------------------------------------------------------\n"
        f"ENA public runs:          {change_text(current, previous, 'ena_public_runs')}\n"
        f"NCBI public runs:         {change_text(current, previous, 'ncbi_public_runs')}\n"
        f"Public in both:           {len(both)}\n"
        f"ENA public, not at NCBI:  {len(ena_only)}\n"
        "\n"
        "PUBLICATION WATCH\n"
        "-----------------------------------------------------------------\n"
        f"Potential delayed samples: {len(publication_issues)} "
        f"(threshold {args.stale_days} days)\n"
        "These are potential delays, not confirmed processing errors.\n"
        "Confirmed FAILED/ACTIVE processing states require authenticated "
        "Webin Reports.\n"
        "\n"
        "The attached TSVs contain project, sample, experiment and run accessions,\n"
        "BioSamples release/update metadata, ENA public status and NCBI status.\n"
    )

    return summary




def build_publication_issues(
    samples,
    runs,
    assembly_rows,
    stale_days
):
    """
    Flag potential ENA publication delays using public archive evidence.

    This is deliberately NOT labelled an ENA processing error. Public Portal
    data cannot expose truly private/failed submissions. Instead, this flags
    cases where a sample/BioSample has already been public for >= stale_days
    but no public run or public assembly is detected for that sample.

    Genuine archival/processing failures should be obtained from the
    authenticated Webin Reports Service.
    """

    runs_by_sample = defaultdict(set)

    for run in runs:
        sample_accession = str(
            run.get(
                "sample_accession",
                ""
            )
        ).strip()

        run_accession = str(
            run.get(
                "accession",
                ""
            )
        ).strip()

        if (
            sample_accession
            and run_accession
        ):
            runs_by_sample[
                sample_accession
            ].add(
                run_accession
            )

    assemblies_by_sample = defaultdict(
        set
    )

    for assembly in assembly_rows:
        sample_accession = str(
            assembly.get(
                "sample_accession",
                ""
            )
        ).strip()

        assembly_accession = str(
            assembly.get(
                "assembly_accession",
                ""
            )
        ).strip()

        if (
            sample_accession
            and assembly_accession
        ):
            assemblies_by_sample[
                sample_accession
            ].add(
                assembly_accession
            )

    issues = []

    for sample_accession, sample in (
        samples.items()
    ):
        ena_public_date = (
            sample.get(
                "first_public",
                ""
            )
        )

        biosamples_release = (
            sample.get(
                "biosamples_release",
                ""
            )
        )

        reference_date = (
            ena_public_date
            or biosamples_release
        )

        parsed = parse_date(
            reference_date
        )

        if not parsed:
            continue

        age_days = (
            date.today()
            - parsed
        ).days

        if age_days < stale_days:
            continue

        public_runs = (
            runs_by_sample.get(
                sample_accession,
                set()
            )
        )

        public_assemblies = (
            assemblies_by_sample.get(
                sample_accession,
                set()
            )
        )

        if (
            public_runs
            or public_assemblies
        ):
            continue

        issues.append({
            "sample_accession":
                sample_accession,

            "project_accession":
                sample.get(
                    "study_accession",
                    ""
                ),

            "project_name":
                sample.get(
                    "project_name",
                    ""
                ),

            "reference_public_date":
                reference_date,

            "days_since_reference_public":
                age_days,

            "public_runs_found":
                0,

            "public_assemblies_found":
                0,

            "issue_type":
                "potential_publication_delay",

            "reason":
                (
                    f"Sample/BioSample has been public for "
                    f"{age_days} days but no public ENA run "
                    f"or assembly was detected."
                ),

            "included_by":
                sample.get(
                    "included_by",
                    ""
                )
        })

    return sorted(
        issues,
        key=lambda row: (
            -int(
                row.get(
                    "days_since_reference_public",
                    0
                )
            ),
            row.get(
                "sample_accession",
                ""
            )
        )
    )


# ------------------------------------------------------------
# Step 12
# Email
# ------------------------------------------------------------

def send_email(
    subject,
    body,
    attachments
):
    """
    Send by configured SMTP when SMTP_HOST is available.

    If SMTP_HOST is not configured, try a local sendmail executable
    (/usr/sbin/sendmail or equivalent). This is useful on institutional
    clusters that provide a local mail relay.

    The function is only called when --send-email is specified.
    """

    recipients = (
        args.email_to
        if args.email_to
        else [DEFAULT_EMAIL]
    )

    message = EmailMessage()

    message["Subject"] = subject
    message["From"] = args.smtp_from
    message["To"] = ", ".join(
        recipients
    )

    message.set_content(
        body
    )

    for filename in attachments:

        with open(
            filename,
            "rb"
        ) as handle:

            message.add_attachment(
                handle.read(),
                maintype="text",
                subtype="tab-separated-values",
                filename=os.path.basename(
                    filename
                )
            )

    # Preferred explicit SMTP configuration.
    if args.smtp_host:

        with smtplib.SMTP(
            args.smtp_host,
            args.smtp_port,
            timeout=60
        ) as server:

            if not args.smtp_no_tls:
                server.starttls()

            if (
                args.smtp_user
                and args.smtp_password
            ):
                server.login(
                    args.smtp_user,
                    args.smtp_password
                )

            server.send_message(
                message
            )

        print(
            "Email sent via SMTP to: "
            + ", ".join(
                recipients
            )
        )

        return

    # Fallback to a local sendmail binary if the cluster provides one.
    sendmail = (
        shutil.which(
            "sendmail"
        )
        or (
            "/usr/sbin/sendmail"
            if os.path.exists(
                "/usr/sbin/sendmail"
            )
            else None
        )
    )

    if sendmail:

        process = subprocess.run(
            [
                sendmail,
                "-t"
            ],
            input=message.as_bytes(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE
        )

        if process.returncode != 0:
            raise RuntimeError(
                "Local sendmail failed: "
                + process.stderr.decode(
                    "utf-8",
                    errors="replace"
                )
            )

        print(
            "Email sent via local sendmail to: "
            + ", ".join(
                recipients
            )
        )

        return

    raise RuntimeError(
        "Email was requested but no mail transport is configured. "
        "Set SMTP_HOST (and SMTP_USER/SMTP_PASSWORD if required), "
        "or use a machine that provides a local sendmail command."
    )

# ------------------------------------------------------------
# Main
# ------------------------------------------------------------

def main():

    print(
        "Starting AEGIS ENA submission tracker..."
    )

    history_file = (
        f"{args.output_prefix}_history.json"
    )

    history = load_history(
        history_file
    )

    previous = None

    if history.get(
        "snapshots"
    ):
        previous = history[
            "snapshots"
        ][-1]

    # 1. Resolve umbrella membership
    umbrella_projects = (
        get_umbrella_projects(
            args.umbrella
        )
    )

    # 2. ENA samples explicitly labelled AEGIS
    aegis_samples = (
        get_project_name_aegis_samples()
    )

    # 3. BioSamples samples explicitly labelled AEGIS
    biosamples_raw = (
        get_biosamples_aegis()
    )

    biosample_ena_samples = (
        get_ena_metadata_for_biosamples(
            biosamples_raw.keys()
        )
    )

    # 4. Samples underneath umbrella
    umbrella_samples = (
        get_umbrella_samples(
            umbrella_projects
        )
    )

    # 5. Merge / deduplicate all three routes
    samples = merge_samples(
        aegis_samples,
        umbrella_samples,
        biosample_ena_samples,
        biosamples_raw,
        umbrella_projects
    )

    print()
    print(
        f"Total unique AEGIS samples "
        f"after merging: {len(samples)}"
    )

    if not samples:
        print(
            "No AEGIS samples found."
        )
        return

    # 6. Runs and experiments
    runs = get_runs_for_samples(
        samples.keys()
    )

    # 6b. Assembly analyses
    assemblies = get_assemblies_for_samples(
        samples.keys(),
        samples,
        umbrella_projects
    )

    ncbi_assembly_status = get_ncbi_assembly_statuses(
        [
            row.get("accession", "")
            for row in assemblies
            if row.get("accession")
        ]
    )

    assembly_rows = build_assembly_table(
        assemblies,
        samples,
        ncbi_assembly_status
    )

    # Fetch components from ENA metadata only.
    # NCBI is intentionally assembly-level only.
    assembly_component_rows = (
        get_ena_assembly_components_metadata(
            assembly_rows,
            samples
        )
    )

    # 7. NCBI public visibility
    ncbi_status = check_ncbi_runs(
        [
            row.get(
                "accession",
                ""
            )
            for row in runs
            if row.get(
                "accession"
            )
        ],
        history=history
    )

    # 8. Detailed table
    records = build_record_table(
        samples,
        runs,
        ncbi_status,
        history
    )

    publication_issues = build_publication_issues(
        samples,
        runs,
        assembly_rows,
        args.stale_days
    )

    # 9. Project summary
    project_summary = (
        build_project_summary(
            records
        )
    )

    # Check NCBI BioProject visibility for all relevant ENA projects:
    # umbrella child projects plus any additional AEGIS projects represented
    # in the merged records. ENA remains authoritative for membership.
    all_relevant_projects = set(umbrella_projects)
    all_relevant_projects.update(
        row.get("project_accession", "")
        for row in project_summary
        if row.get("project_accession")
        and row.get("project_accession") != "NO_PROJECT_ACCESSION"
    )

    ncbi_project_status = get_ncbi_project_statuses(
        all_relevant_projects
    )

    project_status_rows = build_project_status_table(
        all_relevant_projects,
        ncbi_project_status
    )

    # Enrich the existing project summary TSV with project-level status.
    for row in project_summary:
        accession = row.get("project_accession", "")
        status = ncbi_project_status.get(accession, {})
        row["ena_project_status"] = (
            "Public" if accession != "NO_PROJECT_ACCESSION" else "Unknown"
        )
        row["ncbi_project_status"] = status.get(
            "ncbi_project_status",
            "Not checked"
        )
        row["ncbi_project_accession"] = status.get(
            "ncbi_project_accession",
            ""
        )
        row["ncbi_project_uid"] = status.get(
            "ncbi_project_uid",
            ""
        )

    # --------------------------------------------------------
    # Output filenames
    # --------------------------------------------------------

    detailed_file = (
        f"{args.output_prefix}_records.tsv"
    )

    projects_file = (
        f"{args.output_prefix}_projects.tsv"
    )

    project_status_file = (
        f"{args.output_prefix}_project_status.tsv"
    )

    samples_file = (
        f"{args.output_prefix}_samples.tsv"
    )

    biosamples_file = (
        f"{args.output_prefix}_biosamples.tsv"
    )

    assemblies_file = (
        f"{args.output_prefix}_assemblies.tsv"
    )

    assembly_components_file = (
        f"{args.output_prefix}_assembly_components.tsv"
    )

    publication_issues_file = (
        f"{args.output_prefix}_publication_issues.tsv"
    )

    # --------------------------------------------------------
    # Detailed records
    # --------------------------------------------------------

    detailed_fields = [
        "project_accession",
        "project_name",
        "sample_accession",
        "experiment_accession",
        "run_accession",
        "submission_accession",
        "scientific_name",
        "tax_id",
        "first_created",
        "ena_status",
        "ena_public_date",
        "ena_days_since_public",
        "ncbi_status",
        "ncbi_first_seen_public",
        "ncbi_days_since_first_seen",
        "biosamples_release",
        "biosamples_update",
        "last_updated",
        "included_by"
    ]

    write_tsv(
        detailed_file,
        records,
        detailed_fields
    )

    # --------------------------------------------------------
    # Project table
    # --------------------------------------------------------

    project_fields = [
        "project_accession",
        "project_name",
        "ena_project_status",
        "ncbi_project_status",
        "ncbi_project_accession",
        "ncbi_project_uid",
        "samples",
        "experiments",
        "runs",
        "submissions",
        "ena_public_runs",
        "ena_not_public_or_unknown_runs",
        "ncbi_public_runs",
        "ncbi_not_found_or_not_public_runs",
        "included_by"
    ]

    write_tsv(
        projects_file,
        project_summary,
        project_fields
    )

    project_status_fields = [
        "project_accession",
        "ena_project_status",
        "ncbi_project_status",
        "ncbi_project_accession",
        "ncbi_project_uid"
    ]

    write_tsv(
        project_status_file,
        project_status_rows,
        project_status_fields
    )

    # --------------------------------------------------------
    # Sample table
    # --------------------------------------------------------

    sample_rows = []

    for accession, row in sorted(
        samples.items()
    ):

        sample_rows.append({
            "sample_accession":
                accession,

            "secondary_sample_accession":
                row.get(
                    "secondary_sample_accession",
                    ""
                ),

            "sample_alias":
                row.get(
                    "sample_alias",
                    ""
                ),

            "project_accession":
                row.get(
                    "study_accession",
                    ""
                ),

            "project_name":
                row.get(
                    "project_name",
                    ""
                ),

            "submission_accession":
                row.get(
                    "submission_accession",
                    ""
                ),

            "scientific_name":
                row.get(
                    "scientific_name",
                    ""
                ),

            "tax_id":
                row.get(
                    "tax_id",
                    ""
                ),

            "first_public":
                row.get(
                    "first_public",
                    ""
                ),

            "biosamples_release":
                row.get(
                    "biosamples_release",
                    ""
                ),

            "biosamples_update":
                row.get(
                    "biosamples_update",
                    ""
                ),

            "last_updated":
                row.get(
                    "last_updated",
                    ""
                ),

            "included_by":
                row.get(
                    "included_by",
                    ""
                )
        })

    sample_fields = [
        "sample_accession",
        "secondary_sample_accession",
        "sample_alias",
        "project_accession",
        "project_name",
        "submission_accession",
        "scientific_name",
        "tax_id",
        "first_public",
        "biosamples_release",
        "biosamples_update",
        "last_updated",
        "included_by"
    ]

    write_tsv(
        samples_file,
        sample_rows,
        sample_fields
    )

    # --------------------------------------------------------
    # BioSamples-only source table
    # --------------------------------------------------------

    biosample_rows = []

    for accession, row in sorted(
        biosamples_raw.items()
    ):

        merged_row = samples.get(
            accession,
            {}
        )

        biosample_rows.append({
            "biosample_accession":
                accession,

            "biosamples_name":
                row.get(
                    "biosamples_name",
                    ""
                ),

            "project_name":
                "AEGIS",

            "biosamples_release":
                row.get(
                    "biosamples_release",
                    ""
                ),

            "biosamples_update":
                row.get(
                    "biosamples_update",
                    ""
                ),

            "ena_project_accession":
                merged_row.get(
                    "study_accession",
                    ""
                ),

            "ena_submission_accession":
                merged_row.get(
                    "submission_accession",
                    ""
                ),

            "ena_first_public":
                merged_row.get(
                    "first_public",
                    ""
                ),

            "included_by":
                merged_row.get(
                    "included_by",
                    "BioSamples"
                )
        })

    biosample_fields = [
        "biosample_accession",
        "biosamples_name",
        "project_name",
        "biosamples_release",
        "biosamples_update",
        "ena_project_accession",
        "ena_submission_accession",
        "ena_first_public",
        "included_by"
    ]

    write_tsv(
        biosamples_file,
        biosample_rows,
        biosample_fields
    )

    # --------------------------------------------------------
    # Assembly table
    # --------------------------------------------------------

    assembly_fields = [
        "assembly_accession",
        "description",
        "sample_accession",
        "project_accession",
        "project_name",
        "scientific_name",
        "ena_status",
        "ncbi_status",
        "ncbi_accession",
        "ncbi_release_date",
        "ncbi_assembly_status",
        "ncbi_assembly_level",
        "included_by"
    ]

    write_tsv(
        assemblies_file,
        assembly_rows,
        assembly_fields
    )

    assembly_component_fields = [
        "assembly_accession",
        "sample_accession",
        "project_accession",
        "project_name",
        "component_accession",
        "genbank_accession",
        "refseq_accession",
        "component_type",
        "chromosome_name",
        "sequence_name",
        "role",
        "molecule_location_type",
        "length",
        "description",
        "first_public",
        "last_updated",
        "mapping_status",
        "ena_status",
        "ncbi_status"
    ]

    write_tsv(
        assembly_components_file,
        assembly_component_rows,
        assembly_component_fields
    )

    publication_issue_fields = [
        "sample_accession",
        "project_accession",
        "project_name",
        "reference_public_date",
        "days_since_reference_public",
        "public_runs_found",
        "public_assemblies_found",
        "issue_type",
        "reason",
        "included_by"
    ]

    write_tsv(
        publication_issues_file,
        publication_issues,
        publication_issue_fields
    )

    # --------------------------------------------------------
    # Weekly snapshot and console report
    # --------------------------------------------------------

    current = make_snapshot(
        samples,
        runs,
        project_summary,
        umbrella_projects,
        records,
        len(biosamples_raw)
    )

    summary = make_summary(
        samples,
        runs,
        assembly_rows,
        assembly_component_rows,
        publication_issues,
        project_status_rows,
        project_summary,
        umbrella_projects,
        records,
        current,
        previous
    )

    print(
        summary
    )

    print_project_table(
        project_summary
    )

    history.setdefault(
        "snapshots",
        []
    ).append(
        current
    )

    history["snapshots"] = (
        history["snapshots"][-104:]
    )

    save_history(
        history_file,
        history
    )

    print()
    print("OUTPUT FILES")
    print("-" * 65)
    print(
        f"Project summary:  {projects_file}"
    )
    print(
        f"Project status:   {project_status_file}"
    )
    print(
        f"Sample details:   {samples_file}"
    )
    print(
        f"BioSamples:       {biosamples_file}"
    )
    print(
        f"Assemblies:       {assemblies_file}"
    )
    print(
        f"Assembly parts:   {assembly_components_file}"
    )
    print(
        f"Publication watch:{publication_issues_file}"
    )
    print(
        f"All accessions:   {detailed_file}"
    )
    print(
        f"Weekly history:   {history_file}"
    )
    print()

    # --------------------------------------------------------
    # Optional email
    # --------------------------------------------------------

    if args.send_email:

        send_email(
            subject=(
                "AEGIS weekly ENA/NCBI "
                "submission update - "
                f"{date.today().isoformat()}"
            ),
            body=summary,
            attachments=[
                projects_file,
                project_status_file,
                samples_file,
                biosamples_file,
                assemblies_file,
                assembly_components_file,
                publication_issues_file,
                detailed_file
            ]
        )

    else:
        print()
        print(
            "Email was NOT sent because --send-email was not specified."
        )
        print(
            "To send it, rerun with --send-email and configure SMTP_HOST."
        )


if __name__ == "__main__":
    main()
