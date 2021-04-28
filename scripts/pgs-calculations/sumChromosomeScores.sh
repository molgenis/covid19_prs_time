#!/bin/bash
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --ntasks=1

# Load R module
ml R

# Configuration
# A run identifier that is used as a marker for this run.
RUN_IDENTIFIER="<run-identifier>"
# File that contains the summary statistics name corresponding to a GWAS in the first column, and
# the sample size for that GWAS in the second column (comma-separated).
GWASES="<gwas-details-file>"
# Should the trait be skipped if output is already present
SKIP_TRAIT_IF_OUTPUT_PRESENT=true
# The script in this same directory
MERGER_RSCRIPT="./sum-plink-profiles.R"

echo -e "\n[INFO] Configuration:"
echo "--------------"
echo "Run: ${RUN_IDENTIFIER}"
echo "Script to merge profiles: ${MERGER_RSCRIPT}"
echo "Skip trait if the merged profile is already present? ${SKIP_TRAIT_IF_OUTPUT_PRESENT}"
echo -e "--------------\n"

# Specify the output directory
OUTDIR="/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/${RUN_IDENTIFIER}/"

# Loop throug all GWASes to merge the plink profiles for every GWAS.
while read p; do
# Split input line into the summary statistic file and the GWAS sample size.
INPUT=$(echo $p | cut -f1 -d,)
N_GWAS=$(echo $p | cut -f2 -d,)

FILE="${INPUT%.*}"
OUTPUT="${OUTDIR}/${FILE}/full.UGLI.pgs.profile"

echo "-------------------------------------------------------------"
echo "[INFO]	Merging profiles for '$FILE'..."
echo "-------------------------------------------------------------"

# Test if the merged output is already present. If this is the case,
# either skip this or overwrite depending on $SKIP_TRAIT_IF_OUTPUT_PRESENT.

if [ -f $OUTPUT ]; then
  if [[ $SKIP_TRAIT_IF_OUTPUT_PRESENT == true ]]; then
  echo "[WARN] Destination '${FILE}/full.UGLI.pgs.profile' already present. Skipping..."
  continue
  else
    echo "[WARN] Destination '${FILE}/full.UGLI.pgs.profile' already present. Will overwrite final scores..."
  fi
fi

# Check if there is a profile present for every chromosome.
CHROMOSOME_PROFILES_ALL_PRESENT=true

for CHROM in {1..22}; do
  if [[ -f "${OUTDIR}/${FILE}/${CHROM}.UGLI.pgs.sscore" ]]; then
    true
    else
      echo "[WARN] Could not find '${FILE}/${CHROM}.UGLI.pgs.sscore'. Skipping trait..."
    CHROMOSOME_PROFILES_ALL_PRESENT=false
    break
  fi
done

if [[ $CHROMOSOME_PROFILES_ALL_PRESENT == true ]]; then
  echo "[INFO] Rscript $MERGER_RSCRIPT --plink-profiles ${OUTDIR}/${FILE}/*.UGLI.pgs.sscore --out $OUTPUT"
  echo "[INFO] Starting..."
  Rscript $MERGER_RSCRIPT \
  --plink-profiles ${OUTDIR}/${FILE}/*.UGLI.pgs.sscore \
  --scoresum-column "SCORE1_SUM" \
  --out $OUTPUT
fi

done <$GWASES
