#!/bin/bash
#SBATCH --time=05:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=96gb
#SBATCH --nodes=1
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L
#SBATCH --ntasks=1

ml PLINK/2.0-alpha2-20191006

CHR=${1}
runIdentifier=${2}

PLINK_OUTPUT_DIR="/home/umcg-rwarmerdam/lifelines_covid_research/analysis/risky_behaviour/data/cyto/${runIdentifier}/plink"
VCF_PATH="/home/umcg-rwarmerdam/lifelines_covid_research/analysis/risky_behaviour/data/cyto/20210217/annotated/${CHR}.pbwt_reference_impute_renamed_annotated.vcf.gz"

AMBIGUOUS_RSIDS="/home/umcg-rwarmerdam/lifelines_covid_research/analysis/risky_behaviour/data/cyto/${runIdentifier}/ambiguous/${CHR}.cyto.rsids.txt"

mkdir -p $PLINK_OUTPUT_DIR

PLINK_OUTPUT_PREFIX="${PLINK_OUTPUT_DIR}/${CHR}.cyto.imputed.qc.plink"

plink2 --vcf $VCF_PATH dosage=DS \
--maf 0.01 \
--geno 0.25 \
--extract-if-info "INFO >= 0.3" \
--exclude $AMBIGUOUS_RSIDS \
--make-bpgen --out $PLINK_OUTPUT_PREFIX
