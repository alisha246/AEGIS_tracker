# AEGIS ENA/NCBI Submission Tracker

A production-oriented tracker for monitoring AEGIS submissions across ENA, BioSamples, and NCBI.

The tracker discovers AEGIS-related projects and samples, retrieves runs and assemblies, compares public visibility between ENA and NCBI, retrieves ENA assembly-component metadata, tracks weekly changes, identifies potential publication delays, writes TSV reports, and can email a weekly summary.

## Overview

The tracker treats **ENA as the source for AEGIS datasets**.

AEGIS records are discovered through three routes:

1. **ENA sample metadata** where `project_name = "AEGIS"`.
2. **BioSamples** where the characteristic `project name` is `AEGIS`, `Aegis`, or `aegis`.
3. **ENA umbrella components** under umbrella project `PRJEB80366`.

These sources are merged and deduplicated by sample accession. The tracker records why each sample was included.

NCBI is used as an independent public-status comparison for:

- BioProjects
- SRA runs
- genome assemblies

ENA remains authoritative for deciding what belongs to AEGIS.

---

## High-level workflow

```text
ENA umbrella PRJEB80366
        |
        +--> child PRJEB projects ----------------------+
        |                                               |
ENA samples with project_name = AEGIS ----------------+--> merge/deduplicate
        |                                               |        |
BioSamples project name = AEGIS ----------------------+        |
                                                                 v
                                                        all AEGIS samples
                                                           /          \
                                                          /            \
                                                         v              v
                                                     ENA runs      ENA assemblies
                                                         |              |
                                                         v              v
                                                    NCBI SRA      NCBI Datasets
                                                                        |
                                                                        v
                                                          ENA assembly components
                                                           metadata only, no FASTA
                                                                        |
                                                                        v
                                                         publication-delay watch
                                                                        |
                                                                        v
                                                           history + TSV + email
```

## What the tracker monitors

| Record type | ENA | NCBI | Purpose |
|---|---|---|---|
| Projects / BioProjects | Discover and track | Check BioProject visibility | ENA ↔ NCBI propagation |
| Samples | Discover and merge | — | Build complete AEGIS population |
| Experiments | Retrieve through ENA runs | — | Track experiment accessions |
| Runs | Retrieve / public status | SRA visibility check | ENA ↔ NCBI propagation |
| Assemblies | Retrieve / public status | NCBI Datasets check | ENA ↔ NCBI propagation |
| Assembly components | Retrieve sequence metadata | Not checked individually | Chromosomes / scaffolds / contigs |
| Publication watch | Detect potential delays | — | Identify records needing review |

## Requirements

### Python

The tracker requires Python 3 and the packages used by the script, including:

```text
requests
numpy
```

On the Codon environment currently used for production, the Python interpreter is:

```bash
/hps/software/users/cochrane/ena/conda/bin/python
```

On another machine, use an appropriate local Python or virtual environment.

### Network access

The script requires access to:

- ENA Portal API
- ENA Browser API
- BioSamples API
- NCBI BioProject
- NCBI Datasets
- NCBI SRA Data Locator

### Email

The script can send email using either configured SMTP or a local `sendmail` executable.

## Command-line options

| Option | Meaning |
|---|---|
| `-u`, `--umbrella` | Umbrella project accession. Default: `PRJEB80366`. |
| `-o`, `--output-prefix` | Prefix for generated files. Default: `aegis`. |
| `--username` | Optional ENA/Webin username. |
| `--password` | Optional ENA/Webin password. |
| `--send-email` | Send the summary email. |
| `--email-to ADDRESS` | Email recipient. Repeat for multiple recipients. |
| `--smtp-host` | SMTP server, or set `SMTP_HOST`. |
| `--smtp-port` | SMTP port. Default: `587`. |
| `--smtp-user` | SMTP username, or set `SMTP_USER`. |
| `--smtp-password` | SMTP password, or set `SMTP_PASSWORD`. |
| `--smtp-from` | From address. |
| `--smtp-no-tls` | Disable STARTTLS. |
| `--skip-ncbi` | Skip all NCBI public-status checks. |
| `--skip-ncbi-runs` | Skip NCBI SRA run checks but still check NCBI assemblies/projects. |
| `--stale-days N` | Publication-delay threshold in days. Default: `30`. |

## Running manually

A typical production-style run is:

```bash
/hps/software/users/cochrane/ena/conda/bin/python \
  /nfs/production/cochrane/ena/ena-curators/scripts/aegis_submission_tracker.py \
  --output-prefix "$HOME/aegis" \
  --send-email \
  --email-to ahamed@ebi.ac.uk \
  --email-to alishaahamed49@gmail.com
```

For a test without email:

```bash
python aegis_submission_tracker.py --output-prefix aegis
```

## AEGIS sample discovery

### 1. ENA `project_name = "AEGIS"`

ENA is queried for samples explicitly labelled with `project_name = "AEGIS"`.

### 2. BioSamples `project name = AEGIS`

BioSamples is queried using:

```text
attr:project name:AEGIS
attr:project name:Aegis
attr:project name:aegis
```

The tracker verifies the returned value locally using a case-insensitive comparison.

### 3. ENA umbrella components

The tracker retrieves the ENA umbrella project XML for `PRJEB80366` and reads the child project relationships from the umbrella metadata.

The child-project list is cached as:

```text
<output_prefix>_umbrella_projects.txt
```

If the live lookup temporarily fails, the tracker can reuse the last successful cached list.

## Sample merge logic

All discovery routes are combined and deduplicated by sample accession. Each sample retains an `included_by` value showing why it was included, for example:

```text
ENA project_name
BioSamples
umbrella
ENA project_name + BioSamples
ENA project_name + umbrella
ENA project_name + BioSamples + umbrella
```

## Run and experiment retrieval

For every detected AEGIS sample, the script queries ENA `read_run` records.

ENA queries are batched. Temporary `500`, `502`, `503`, or `504` responses are retried automatically. If a large run batch still fails, the tracker splits it into smaller groups and retries rather than aborting the whole run.

## NCBI SRA run visibility

The tracker checks whether ENA run accessions are also resolvable/public in NCBI SRA.

The current implementation uses the **NCBI SRA Data Locator** in batches. Previous successful NCBI results are reused from history so already-confirmed public runs do not need to be checked again every week.

## Assembly retrieval

Assemblies are fetched from ENA using detected BioSamples-style accessions such as `SAMEA...`, `SAMN...`, `SAMD...`, and `SAMS...`.

Assemblies are deduplicated by public assembly accession.

## NCBI assembly status

ENA assemblies are checked independently in **NCBI Datasets**.

The tracker normalizes GCA/GCF assembly accessions by removing an optional version suffix before comparison, so values such as:

```text
GCA_977066645
GCA_977066645.1
GCA_977066645.2
```

can be matched on the same base accession.

## NCBI project status

Projects discovered from ENA are checked against **NCBI BioProject**.

ENA remains the source for AEGIS datasets. The NCBI check only asks whether each ENA-derived project is also visible/public in NCBI.

Typical project-status fields include:

```text
project_accession
ena_project_status
ncbi_project_status
ncbi_project_accession
ncbi_project_uid
```

## Assembly components

Assembly components are retrieved from ENA using **sequence metadata only**. The tracker does not download assembly FASTA files.

This keeps the component stage much faster and avoids large memory/network costs.

If a sample has exactly one assembly, component rows can be mapped directly to that GCA. If a sample has multiple assemblies, the tracker reports:

```text
mapping_status = ambiguous_multiple_assemblies_for_sample
```

rather than guessing.

## Publication-delay watch

The default watch condition is:

```text
sample/BioSample has a public date
AND
days since that date >= --stale-days
AND
no public ENA run exists
AND
no public ENA assembly exists
```

Default threshold: `30` days.

These are **potential publication delays, not confirmed ENA submission errors**.

## History and weekly change tracking

The script loads:

```text
<output_prefix>_history.json
```

at startup.

If a previous snapshot exists, it becomes the comparison baseline. At the end of a successful run, a new snapshot is appended.

The tracker keeps the most recent **104 snapshots**, corresponding to approximately two years of weekly history.

When previous history is available, the summary can display values such as:

```text
Runs: 2409 (+3)
Samples: 3799 (+0)
Assemblies: 509 (+2)
```

History is normally per user because the file location depends on `--output-prefix`.

For example:

```bash
--output-prefix "$HOME/aegis"
```

means each user has their own:

```text
$HOME/aegis_history.json
```

If another user runs the tracker for testing and has no previous history file, their first report will not contain `(+N)` values. This is expected.

Deleting or replacing the history file resets the comparison history.

History is also used to reduce NCBI traffic because runs previously confirmed public do not need to be queried again on every execution.

## Output files

With `--output-prefix "$HOME/aegis"`, typical outputs include:

```text
$HOME/aegis_history.json
$HOME/aegis_umbrella_projects.txt
$HOME/aegis_projects.tsv
$HOME/aegis_project_status.tsv
$HOME/aegis_samples.tsv
$HOME/aegis_biosamples.tsv
$HOME/aegis_assemblies.tsv
$HOME/aegis_assembly_components.tsv
$HOME/aegis_publication_issues.tsv
$HOME/aegis_detailed.tsv
```

## Email reporting

Email is sent only when `--send-email` is provided.

Multiple recipients are supported by repeating `--email-to`, for example:

```bash
--email-to ahamed@ebi.ac.uk \
--email-to alishaahamed49@gmail.com
```

SMTP is used when configured; otherwise the script can fall back to a local `sendmail` executable.

Slurm job notifications are separate from the AEGIS report email. Python report recipients are controlled by `--email-to`; Slurm BEGIN/END/FAIL mail is controlled by `#SBATCH --mail-user`.

## Slurm execution

The shared production wrapper is:

```text
/nfs/production/cochrane/ena/ena-curators/scripts/aegis_weekly.sh
```

A portable wrapper should use the submitting user's home directory:

```bash
#!/bin/bash
#SBATCH --job-name=aegis_tracker
#SBATCH --output=%x_%j.out
#SBATCH --error=%x_%j.err

cd "$HOME" || exit 1

/hps/software/users/cochrane/ena/conda/bin/python \
  /nfs/production/cochrane/ena/ena-curators/scripts/aegis_submission_tracker.py \
  --output-prefix "$HOME/aegis" \
  --send-email
```

Submit with:

```bash
sbatch /nfs/production/cochrane/ena/ena-curators/scripts/aegis_weekly.sh
```

Useful Slurm commands:

```bash
squeue -u $USER
sacct -u $USER --starttime today
scancel <JOBID>
```

## Weekly cron automation

Not every user needs their own cron entry.

The official weekly scheduler can live on one machine/account.

Current weekly schedule:

```cron
0 9 * * 1 /bin/bash -c '/ebi/slurm/codon/bin/squeue -h -u "$USER" -n aegis_tracker | /usr/bin/grep -q . || /ebi/slurm/codon/bin/sbatch /nfs/production/cochrane/ena/ena-curators/scripts/aegis_weekly.sh'
```

This runs every Monday at 09:00 and avoids overlapping `aegis_tracker` jobs.

## Troubleshooting

### ENA `500 Server Error`

Usually a transient server or query-size issue. The current tracker retries automatically, and run batches can be split into smaller queries.

### `RemoteDisconnected`

The remote service closed the connection unexpectedly. Retrying usually resolves it. The umbrella lookup also has a cache fallback.

### NCBI `429 Too Many Requests`

This is NCBI rate limiting, especially on a shared HPC IP. The current run checker uses batching and history caching to reduce requests.

### `python: command not found` in Slurm

Use the full interpreter path:

```bash
/hps/software/users/cochrane/ena/conda/bin/python
```

### Permission denied writing `aegis_umbrella_projects.txt`

The job is running from a non-writable directory.

Use:

```bash
cd "$HOME" || exit 1
```

and preferably:

```bash
--output-prefix "$HOME/aegis"
```

### Missing `(+0)` / `(+N)` values

The current user has no previous history snapshot. This is expected on another user's first test run.

### Many ambiguous component mappings

The sample has multiple assemblies and metadata-only component retrieval cannot uniquely assign each component to one GCA. This is expected behavior.

## Key functions

| Function | Responsibility |
|---|---|
| `chunks()` | Split values into API-sized groups |
| `make_auth()` | Optional ENA/Webin authentication |
| `query_ena()` | Generic ENA Portal query with retries |
| `get_biosamples_aegis()` | Find AEGIS BioSamples |
| `get_umbrella_projects()` | Resolve ENA umbrella child projects |
| `get_project_name_aegis_samples()` | Find ENA `project_name=AEGIS` samples |
| `get_umbrella_samples()` | Find samples under umbrella child projects |
| `get_ena_metadata_for_biosamples()` | Enrich BioSamples records from ENA |
| `merge_samples()` | Merge and deduplicate all sample sources |
| `get_runs_for_samples()` | Retrieve ENA runs and experiments |
| `check_ncbi_runs()` | Check NCBI SRA run visibility |
| `get_assemblies_for_samples()` | Retrieve ENA assemblies |
| `get_ncbi_assembly_statuses()` | Check assemblies in NCBI Datasets |
| `get_ena_assembly_components_metadata()` | Retrieve ENA sequence/component metadata |
| `build_publication_issues()` | Flag potential publication delays |
| `load_history()` / `save_history()` | Persist weekly history |
| `build_project_summary()` | Aggregate project-level data |
| `make_summary()` | Build weekly text summary |
| `send_email()` | Send summary and attachments |
| `main()` | Orchestrate the workflow |

## Maintenance notes

- ENA should remain the authority for defining AEGIS ownership.
- NCBI should remain a comparison/propagation check rather than the discovery source.
- Do not interpret API failures as evidence that records are private.
- Preserve the production history and umbrella cache between weekly runs.
- Test API-field changes with a small known accession before running the full workflow.
- Avoid reintroducing full assembly FASTA downloads for component discovery unless necessary.
- Treat publication-watch rows as investigation candidates, not confirmed submission errors.

## Git workflow

After updating the tracker or README:

```bash
git status
git add aegis_submission_tracker.py README.md
git commit -m "Update AEGIS submission tracker"
git push
```

If using SSH with GitHub:

```bash
git clone git@github.com:alisha246/AEGIS_tracker.git
```

## Important interpretation note

This tracker is a **monitoring and reporting tool**.

It reads metadata, compares public visibility, stores history, and highlights possible delays.

It does **not** modify ENA, BioSamples, or NCBI records.
