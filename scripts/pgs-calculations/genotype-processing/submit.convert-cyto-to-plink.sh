runIdentifier="20210217"
jobpath="/home/umcg-rwarmerdam/lifelines_covid_research/analysis/risky_behaviour/jobs/cyto-processing/convert-cyto-to-plink.sh"
logs=$(readlink -f "/home/umcg-rwarmerdam/lifelines_covid_research/analysis/risky_behaviour/logs/cyto-processing/${runIdentifier}/convert_cyto_to_plink")

mkdir -p ${logs}

for chr in {1..22}
  do

  sbatch -J "chr${chr}.convert_cyto_to_plink" -o "${logs}/${chr}.convert_cyto_to_plink.out" -e "${logs}/{chr}.convert_cyto_to_plink.err" -v ${jobpath} ${chr} ${runIdentifier}

  sleep 0.5
done
