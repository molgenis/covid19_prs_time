#!/bin/bash
#SBATCH --job-name=annotate-cyto-with-dbSNP
#SBATCH --output=../../logs/cyto-processing/20210217/annotate/%x_%a.out
#SBATCH --error=../../logs/cyto-processing/20210217/annotate/%x_%a.err
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=16gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --ntasks=1
#SBATCH --array=1-22

ml HTSlib/1.9-foss-2018b

CHR=${SLURM_ARRAY_TASK_ID}
#CHR=1
runIdentifier="20210217"

outDir="/home/umcg-rwarmerdam/lifelines_covid_research/analysis/risky_behaviour/data/cyto/${runIdentifier}/annotated"
outFile="${outDir}/${CHR}.pbwt_reference_impute_renamed_annotated.vcf.gz"
inputVcf="/groups/umcg-lifelines/tmp01/releases/cytosnp_imputed/v5/sanger_vcf/${CHR}.pbwt_reference_impute_renamed.vcf.gz"
annotations="/home/umcg-rwarmerdam/lifelines_covid_research/analysis/risky_behaviour/data/cyto/${runIdentifier}/annotations/${CHR}_variants+RSIDs+header.tmp.txt"

bgzip -c ${annotations} > "${annotations}.gz"
tabix -s1 -b2 -e2 "${annotations}.gz"

module purge
ml BCFtools/1.11-GCCcore-7.3.0

mkdir -p ${outDir}

bcftools annotate \
--annotations "${annotations}.gz" \
--columns CHROM,POS,REF,ALT,ID \
--output-type z ${inputVcf} > ${outFile}
