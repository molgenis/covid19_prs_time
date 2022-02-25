#!/bin/bash
#SBATCH --job-name=chr%a.query-variants-from-cyto
#SBATCH --output=query-variants-from-cyto/%A_%a.out
#SBATCH --error=query-variants-from-cyto/%A_%a.err
#SBATCH --time=01:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --ntasks=1
#SBATCH --array=1-22

ml BCFtools/1.11-GCCcore-7.3.0

CHR=${SLURM_ARRAY_TASK_ID}
runIdentifier="20210217"

outDir="<directory_to_write_variants_to>"
outFile="${outDir}/${CHR}_variants.tmp.txt"
inputVcf="/groups/umcg-lifelines/tmp01/releases/cytosnp_imputed/v5/sanger_vcf/${CHR}.pbwt_reference_impute_renamed.vcf.gz"

mkdir -p ${outDir}

bcftools query -f '%CHROM\t%POS\t%REF\t%ALT{0}\n' ${inputVcf} > ${outFile}
