#!/bin/bash
#SBATCH --time=0-23:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=196gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --ntasks=9

# Configuration
SKIP_TRAIT_IF_OUT_DIR_PRESENT=false
SKIP_PRSCS_FOR_CHROM_IF_OUTPUT_PRESENT=true
# A run identifier that is used as a marker for this run.
RUN_IDENTIFIER="<run-identifier>"
SINGLE_CHROMOSOME_SCRIPT="./calculatePolygenicScoresForChromosome.sh"

echo -e "\nConfiguration:"
echo "--------------"
echo "Run: ${RUN_IDENTIFIER}"
echo "Skip trait if the corresponding output directory is already present? ${SKIP_TRAIT_IF_OUT_DIR_PRESENT}"
echo "Skip chromosome if there is already output found for the specific chromosome? ${SKIP_PRSCS_FOR_CHROM_IF_OUTPUT_PRESENT}"
echo -e "--------------\n"

BFILE_PATH="/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/data/UGLI/20200811/plink/"
SUMMARY_DIR="/groups/umcg-lld/tmp01/other-users/umcg-rwarmerdam/pgs_based_mixup_correction/data/summarystatistics/PRScs/"
# File that contains the summary statistics name corresponding to a GWAS in the first column, and
# the sample size for that GWAS in the second column (comma-separated).
GWASES="<gwas-details-file>"
OUTDIR="/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/output/PRScs/${RUN_IDENTIFIER}/"
LOGS="/groups/umcg-lifelines/tmp01/projects/ugli_blood_gsa/pgs_based_mixup_correction/logs/PRScs/${RUN_IDENTIFIER}/"

while read p; do
# Split input line into the summary statistic file and the GWAS sample size
INPUT=$(echo $p | cut -f1 -d,)
N_GWAS=$(echo $p | cut -f2 -d,)

FILE="${INPUT%.*}"
SUMMARY=$SUMMARY_DIR/$INPUT

echo "-------------------------------------------------------------"
echo "[INFO]	Submitting jobs for '$FILE' (n = ${N_GWAS})..."
echo "-------------------------------------------------------------"
if [ "$(ls -A $OUTDIR/$FILE)" ]; then
if [[ $SKIP_TRAIT_IF_OUT_DIR_PRESENT == true ]]; then
echo "[WARN] $FILE is not empty. skipping..."
continue
else
  echo "[WARN] $FILE is not empty; will overwrite final scores."
fi
fi

mkdir -p $OUTDIR/$FILE
mkdir -p $LOGS/$FILE

for CHROM in {1..22}; do

# When there are already files matching the PRScs output for this chromosome, skip this step
if [[ -f "${OUTDIR}/${FILE}/${CHROM}.UGLI.pgs.sscore" && $SKIP_PRSCS_FOR_CHROM_IF_OUTPUT_PRESENT == true ]]; then
echo "File '${FILE}/${CHROM}.UGLI.pgs.sscore' is already present! Skipping..."
else
  echo "Starting processing chromosome ${CHROM}, summary statistics '${FILE}'..."

sbatch --time=00:59:00 --mem=16gb --ntasks=1 \
-J "${RUN_IDENTIFIER}.chr${CHROM}.${FILE}.calculate-polygenic-scores" \
-o "${LOGS}/${FILE}/${CHROM}.calculate-polygenic-scores.out" \
-e "${LOGS}/${FILE}/${CHROM}.calculate-polygenic-scores.err" \
$SINGLE_CHROMOSOME_SCRIPT \
$BFILE_PATH \
$CHROM \
$SUMMARY \
$N_GWAS \
$OUTDIR \
$FILE </dev/null \

fi

#sleep 0.5
done

done <$GWASES

wait
