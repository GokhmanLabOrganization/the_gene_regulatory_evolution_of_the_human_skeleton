#!/usr/bin/env bash
set -euo pipefail

if [ "$#" -lt 3 ]; then
    echo "Usage: $0 GFF_PATH EXCLUDE_SPEC ROOT_DIR [ROOT_DIR2 ...]" >&2
    echo "Example: $0 genes.gtf \"chrM,chrUn\" /data/project1 /data/project2" >&2
    exit 1
fi

GFF="$1"
EXCLUDE="$2"
shift 2  # now "$@" is the list of root paths

for ROOT in "$@"; do
    # find all matching files under this root
    find "$ROOT" -type f -name '*_ase_by_reads_merged.txt' | while IFS= read -r COUNTS; do
        DIR=$(dirname "$COUNTS")
        FNAME=$(basename "$COUNTS")

        # sample_name is everything before the suffix
        SAMPLE_NAME=${FNAME%_ase_by_reads_merged.txt}

        OUT="${DIR}/${SAMPLE_NAME}_filtered_counts_with_tpm.txt"

        echo "Processing: $COUNTS"
        echo "Output    : $OUT"

        python TPM_calculation.py \
            --counts "$COUNTS" \
            --gff "$GFF" \
            --exclude "$EXCLUDE" \
            --out "$OUT"
    done
done
echo "All done!"


# bash run_tpm_all.sh \
#     /path/to/genes.gtf \
#     "chrM,chrUn" \
#     /proj/ase_runs/batch1 \
#     /proj/ase_runs/batch2