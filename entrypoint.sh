#! /bin/bash --login

## entrypoint.sh is modified from https://pythonspeed.com/articles/activate-conda-dockerfile/

# The --login ensures the bash configuration is loaded,
# enabling Conda.

# Enable strict mode.
set -euo pipefail
# ... Run whatever commands ...

# Temporarily disable strict mode and activate conda:
set +euo pipefail
source $(conda info --base)/etc/profile.d/conda.sh
conda activate starscope_scRNAseq_env

# Re-enable strict mode:
set -euo pipefail

# exec the final command:
exec "$@"
