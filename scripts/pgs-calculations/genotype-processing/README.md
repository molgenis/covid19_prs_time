# Readme
The genetic data was preprocessed by excluding ambiguous SNPs and converting 
the genotype data from VCF to the hybrid PLINK 2.0 bpgen format, maintaining 
allelic dosage information. Variants were removed if either their minor 
allele frequency (MAF) was less than 0.01, imputation score was below 0.3, 
or missing call rates exceeded 0.25.

- `query-vcf-variants.sh`: Queries variants from vcf file
- `annotate-variants.sh`: Matches variants from dbsnp to queried variants
- `add-header-lines.sh`: Adds header lines to the output from the previous step
- `annotate-vcf-with-dbSnp.sh`: Replaces variant IDs with matched variants ids from dbsnp
- `submit.extract-ambiguous-cyto-SNPs.sh`: Extracts ambiguous variants
- `submit.convert-cyto-to-plink.sh`: Removes the ambiguous variants and performs other variant filtering steps
