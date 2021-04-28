#!/bin/bash

# Script that calculates polygenic scores for a single trait and a single chromosome.
# This should be called using './calculatePolygenicScores.sh'

# Load modules
ml PLINK/2.0-alpha2-20191006
ml PythonPlus/2.7.16-foss-2018b-v20.02.1
ml HDF5/1.10.6-foss-2018b

# Activate the virtual environment that contains the PRScs dependencies.
source /groups/umcg-lld/tmp01/other-users/umcg-rwarmerdam/pgs_based_mixup_correction/tools/PRScs_venv/bin/activate

set -e

BFILE_PATH=${1}
CHROM=${2}
SUMMARY=${3}
N_GWAS=${4}
OUTDIR=${5}
FILE=${6}
$LD_REFERENCE_PANEL="/groups/umcg-lld/tmp01/other-users/umcg-rwarmerdam/pgs_based_mixup_correction/data/PRScs/reference_panel/ldblk_1kg_eur/"

echo -e "\nStarting:"
echo "--------------"
echo "Processing chromosome ${CHROM}, for trait '$FILE', using the following paramters:"
echo "BED file prefix = ${BFILE_PATH}"
echo "Summary statistics file = ${SUMMARY}"
echo "GWAS Sample size = ${N_GWAS}"
echo "Output directory = ${OUTDIR}"
echo -e "--------------\n"

# BFILE
BFILE="${BFILE_PATH}/${CHROM}.UGLI.imputed.qc.plink"

python ~/pgs_based_mixup_correction/tools/PRScs/PRScs.py \
--ref_dir=$LD_REFERENCE_PANEL \
--bim_prefix=$BFILE \
--sst_file=$SUMMARY \
--n_gwas=$N_GWAS \
--out_dir=$OUTDIR/$FILE/ \
--chrom=$CHROM

echo "Finished PRScs step"

PRScs_OUTPUT="${OUTDIR}/${FILE}/*_chr${CHROM}.txt"

echo "Starting PLINK allelic scoring command using '${PRScs_OUTPUT}'..."

# Perform plink scoring for this chromosome
plink2 --bpfile $BFILE \
--score $PRScs_OUTPUT 2 4 6 ignore-dup-ids cols=+scoresums \
--out "${OUTDIR}/${FILE}/${CHROM}.UGLI.pgs"

echo "Finished PLINK allelic scoring command"
echo "DONE, finished calculating PGSes for chromosome ${CHROM} for trait '${FILE}'"

deactivate
