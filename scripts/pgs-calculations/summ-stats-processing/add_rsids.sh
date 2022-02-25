#!/bin/bash
#SBATCH --job-name=add_rsids
#SBATCH --output=/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/logs/add_rsids/COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.b37/20210211.out
#SBATCH --error=/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/logs/add_rsids/COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.b37/20210211.err
#SBATCH --time=00:59:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=20gb
#SBATCH --nodes=1
#SBATCH --qos=priority
#SBATCH --open-mode=append
#SBATCH --export=NONE
#SBATCH --get-user-env=L

module load Python/3.7.4-GCCcore-7.3.0-bare

sumStatsBaseDir="/groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/data/summarystatistics"
sumStatsFile="COVID19_HGI_B2_ALL_eur_leave_23andme_20210107.b37.txt"

chromosomeIndex=1
basepairIndex=2
effectAlleleIndex=4
nonEffectAlleleIndex=3
effectSizeIndex=7
pValueIndex=9
newSnpIndex=$(expr 1 + $(awk -F'\t' '{print NF; exit}' $sumStatsBaseDir/raw/$sumStatsFile))

# First grep the RSIDs corresponding to each SNP
python -u /groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/jobs/add_rsids/rsid-annotation-small_V2.py \
-d /groups/umcg-lifelines/tmp01/projects/ov20_0554/analysis/risky_behaviour/data/dbsnp/GCF_000001405.25 \
-s $sumStatsBaseDir/raw/$sumStatsFile \
-CHR $chromosomeIndex -BP $basepairIndex -A1 $effectAlleleIndex -A2 $nonEffectAlleleIndex \
> $sumStatsBaseDir/matched_rsids/$sumStatsFile

#exit 0

# Paste the summary statistics together with the matched RSIDs.
# make sure that the file with matched rsids does not contain 2 header lines

# Thereafter pick the RSID, effect allele (A1), ref allele (A2), beta, p-val columns
# Use sed to rename the header to something that PRScs can understand

paste $sumStatsBaseDir/raw/$sumStatsFile \
$sumStatsBaseDir/matched_rsids/$sumStatsFile \
| awk -v FS='\t' -v OFS='\t' -v ES="$effectSizeIndex" -v P="$pValueIndex" -v A1="$effectAlleleIndex" -v A2="$nonEffectAlleleIndex" -v SNP="$newSnpIndex" \
'{print $SNP,$A1,$A2,$ES,$P}' \
| sed '1 s/^.*$/SNP\tA1\tA2\tBETA\tP/' > \
$sumStatsBaseDir/PRScs/$sumStatsFile
