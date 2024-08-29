# these codes should be supported by plink software
# direct to the Imperial hpc server
ssh xxxxx@login.hpc.imperial.ac.uk
# load the software plink
module load plink

# Perform association analysis and Interpret the results
plink --bfile 	pcsk9_extraction --extract missense_variant.txt --make-bed --out missense_dataset
# perform linear and logistic regression
plink --bfile missense_dataset --linear --hide-covar --pheno T2D_depr_BMItransformed_RG_final.sample --pheno-name RG_ln --ci 0.95  --covar T2D_depr_BMItransformed_RG_final.sample --covar-name PC1,PC2,PC3,PC4,PC5,PC6 --out rg_pcsk9_only6PCs
plink --bfile missense_dataset --logistic --hide-covar --pheno T2D_depr_BMItransformed_RG_final.sample --pheno-name T2D --ci 0.95 --covar T2D_depr_BMItransformed_RG_final.sample --covar-name PC1,PC2,PC3,PC4,PC5,PC6 --out t2d_pcsk9_only6PCs
# save the results to tsv format
awk '{$1=$1}1' OFS='\t' t2d_pcsk9_only6PCs.assoc.logistic > t2d_pcsk9_only6PCs.assoc.tsv
awk '{$1=$1}1' OFS='\t' rg_pcsk9_only6PCs.assoc.linear > t2d_pcsk9_only6PCs.assoc.tsv

