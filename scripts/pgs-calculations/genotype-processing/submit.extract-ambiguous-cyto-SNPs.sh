runIdentifier="20210217"
jobpath="/home/umcg-rwarmerdam/lifelines_covid_research/analysis/risky_behaviour/jobs/cyto-processing/extract-ambiguous-cyto-SNPs.sh"
logs=$(readlink -f "/home/umcg-rwarmerdam/lifelines_covid_research/analysis/risky_behaviour/logs/cyto-processing/${runIdentifier}/ambiguous-cyto-SNPs/")

mkdir -p ${logs}

for chr in {1..22}
  do

  sbatch -J "chr${chr}.extract-ambiguous-SNPs" -o "${logs}/${chr}.extract-ambiguous-SNPs.out" -e "${logs}/${chr}.extract-ambiguous-SNPs.err" -v ${jobpath} ${chr} ${runIdentifier}

  sleep 0.5
done
