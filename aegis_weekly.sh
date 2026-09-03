#!/bin/bash
#SBATCH --job-name=aegis_weekly
#SBATCH --output=aegis_weekly_%j.out
#SBATCH --error=aegis_weekly_%j.err
#SBATCH --time=06:00:00
#SBATCH --mem=8G

set -euo pipefail

cd "$HOME" || exit 1

# If needed, activate the environment that has requests/numpy installed.
# source .venv/bin/activate

/hps/software/users/cochrane/ena/conda/bin/python /nfs/production/cochrane/ena/ena-curators/scripts/aegis_submission_tracker.py \
    --output-prefix "$HOME/aegis" \
    --send-email \
    --email-to ahamed@ebi.ac.uk \
    --email-to alishaahamed49@gmail.com \
    --email-to ihsan@ebi.ac.uk \
    --email-to joanap@ebi.ac.uk \
    --email-to woollard@ebi.ac.uk \
    --email-to mar@ebi.ac.uk
