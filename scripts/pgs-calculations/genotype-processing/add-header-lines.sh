#!/bin/bash
#SBATCH --job-name=annotate-variants-from-cyto
#SBATCH --output=../../logs/cyto-processing/20210217/annotate/%A_%a.out
#SBATCH --error=../../logs/cyto-processing/20210217/annotate/%A_%a.err
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
runIdentifier="20210217"

variantsDir="/home/umcg-rwarmerdam/lifelines_covid_research/analysis/risky_behaviour/data/cyto/${runIdentifier}/annotations"
variantsFile="${variantsDir}/${CHR}_variants.tmp.txt"
annotatedVariantsFile="${variantsDir}/${CHR}_variants+RSIDs.tmp.txt"
variantIdsFile="${variantsDir}/${CHR}_RSIDs.tmp.txt"

#sed  -i '1i\CHROM\tPOS\tREF\tALT\tSNP' $annotatedVariantsFile > "${variantsDir}/${CHR}_variants+RSIDs+header.tmp.txt"
#paste $variantsFile $variantIdsFile > $annotatedVariantsFile

echo -e "CHROM\tPOS\tREF\tALT\tID" | cat - $annotatedVariantsFile > "${variantsDir}/${CHR}_variants+RSIDs+header.tmp.txt"

#exit 0

# Paste the variant list together with the matched RSIDs.
# make sure that the file with matched rsids does not contain 2 header lines

# Thereafter pick the RSID, effect allele (A1), ref allele (A2), beta, p-val columns
# Use sed to rename the header to something that PRScs can understand

#paste $variantsFile $variantIdsFile > $annotatedVariantsFile
