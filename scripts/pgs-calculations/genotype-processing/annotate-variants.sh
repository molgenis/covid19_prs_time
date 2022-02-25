#!/bin/bash
#SBATCH --job-name=annotate-variants-from-cyto
#SBATCH --output=annotate/%A_%a.out
#SBATCH --error=annotate/%A_%a.err
#SBATCH --time=01:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --ntasks=1
#SBATCH --array=1-22

module load Python/3.7.4-GCCcore-7.3.0-bare

CHR=${SLURM_ARRAY_TASK_ID}

variantsDir="<directory_with_queried_variant_files>"
variantsFile="${variantsDir}/${CHR}_variants.tmp.txt"
annotatedVariantsFile="${variantsDir}/${CHR}_variants+RSIDs.tmp.txt"
variantIdsFile="${variantsDir}/${CHR}_RSIDs.tmp.txt"

python -u ./add_rsids/rsid-annotation-small_V2.py \
-d "<snpdb_vcf_file>" \
-s $variantsFile \
-CHR 1 -BP 2 -A1 3 -A2 4 \
> $variantIdsFile

#exit 0

# Paste the variant list together with the matched RSIDs.
# make sure that the file with matched rsids does not contain 2 header lines

# Thereafter pick the RSID, effect allele (A1), ref allele (A2), beta, p-val columns
# Use sed to rename the header to something that PRScs can understand

paste $variantsFile $variantIdsFile > $annotatedVariantsFile
