# direct to the Imperial hpc server
ssh xxxxx@login.hpc.imperial.ac.uk
# load the software plink
module load plink
# Recode vcf to plink
plink --vcf b36_c1_extraction_gz.vcf --recode --out b36_c1_extraction_plink
# convert the PED/MAP files to PLINK binary format
plink --file b36_c1_extraction_plink --make-bed --out b36_c1_extraction_plink
plink --vcf b36_c1_extraction_gz.vcf --make-bed --out b36_c1_test

# Check frequencies of the variants using  --freq and/or --frqx
# generate allele counts
plink --bfile pcsk9_extraction --freq counts --out allele_counts

# generate a file with minor allele frequencies
plink --bfile pcsk9_extraction --freq --out variant_frequencies
plink --bfile pcsk9_extraction --hardy --out hwe_results

