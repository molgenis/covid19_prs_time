#!/bin/bash
#SBATCH --time=01:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=8gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --ntasks=1

CHR=${1}
runIdentifier=${2}

OUT_DIR="/home/umcg-rwarmerdam/lifelines_covid_research/analysis/risky_behaviour/data/cyto/${runIdentifier}/ambiguous"
OUT_RSID_LIST="${OUT_DIR}/${CHR}.cyto.rsids.txt"
inputVcf="/home/umcg-rwarmerdam/lifelines_covid_research/analysis/risky_behaviour/data/cyto/20210217/annotated/${CHR}.pbwt_reference_impute_renamed_annotated.vcf.gz"

mkdir -p ${OUT_DIR}

gzip -cd ${inputVcf} | \
awk '$1 ~ /^#/ {next} {if (($4 ~ /C|G/ && $5 ~ /C|G/) || ($4 ~ /A|T/ && $5 ~ /A|T/)) print $3}' > $OUT_RSID_LIST
