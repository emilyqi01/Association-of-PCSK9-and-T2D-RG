import pandas as pd
import numpy as np
# Load the allele frequencies file
freq_file = 'allele_frequencies.frq' 
freq_df = pd.read_csv(freq_file, delim_whitespace=True)
# Extract minor and major allele information
allele_dict = {row['SNP']: (row['A1'], row['A2']) for _, row in freq_df.iterrows()}
# load the phenotype data
phenotype_file = 'phenotype.sample'  
phenotype_df = pd.read_csv(phenotype_file, sep=r"\s+")
# Load the MAP file
map_file = 'output_missense.map' 
map_df = pd.read_csv(map_file, delim_whitespace=True, header=None)
map_df.columns = ['Chromosome', 'SNP', 'Genetic_Distance', 'Base_Pair_Position']
# Load the PED file
ped_file = 'output_missense.ped' 
ped_df = pd.read_csv(ped_file, delim_whitespace=True, header=None)
# Define the column names for PED file
n_snps = (ped_df.shape[1] - 6) // 2
ped_columns = ['Family_ID', 'Individual_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phenotype']
for i in range(n_snps):
    ped_columns.append(f'SNP_{i+1}_A1')
    ped_columns.append(f'SNP_{i+1}_A2')
ped_df.columns = ped_columns
# Merge SNP_A1 and SNP_A2 columns into a single column SNP_{i+1}
for i in range(n_snps):
    ped_df[f'SNP_{i+1}'] = ped_df[f'SNP_{i+1}_A1'] + ped_df[f'SNP_{i+1}_A2']

# Drop the original SNP_A1 and SNP_A2 columns
ped_df = ped_df[['Family_ID', 'Individual_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phenotype'] + [f'SNP_{i+1}' for i in range(n_snps)]]

# Vectorized function to encode genotype data based on minor and major alleles
def encode_genotype_vectorized(genotype, minor_allele, major_allele):
    encoding = {
        f'{minor_allele}{minor_allele}': 0,
        f'{minor_allele}{major_allele}': 1,
        f'{major_allele}{minor_allele}': 1,
        f'{major_allele}{major_allele}': 2,
        '00': np.nan, '0A': np.nan, '0C': np.nan, '0G': np.nan, '0T': np.nan,
        f'{minor_allele}0': np.nan, f'{major_allele}0': np.nan,
        'A0': np.nan, 'C0': np.nan, 'G0': np.nan, 'T0': np.nan,
        f'0{minor_allele}': np.nan, f'0{major_allele}': np.nan
    }
    # Debug: Print the genotype and its encoding
   #print(f'Genotype: {genotype}, Minor Allele: {minor_allele}, Major Allele: {major_allele}, Encoding: {encoding.get(genotype, np.nan)}')
    return encoding.get(genotype, np.nan)

# Combine and encode SNP columns
encoded_snps = []
for i in range(n_snps):
    snp_col = f'SNP_{i+1}'
    snp_name = map_df.iloc[i, 1]
    alleles = allele_dict.get(snp_name, None)
    if alleles:
        minor_allele, major_allele = alleles
        encoded_col = f'Encoded_{snp_col}'
        # Debug: Print the SNP column and minor allele
        ped_df[encoded_col] = ped_df[snp_col].apply(lambda x: encode_genotype_vectorized(x, minor_allele, major_allele))
        encoded_snps.append(encoded_col)

# Drop the original merged SNP columns
ped_df = ped_df[['Family_ID', 'Individual_ID', 'Paternal_ID', 'Maternal_ID', 'Sex', 'Phenotype'] + encoded_snps]

# Merge phenotype data with PED data
merged_df = pd.merge(ped_df, phenotype_df, left_on=['Family_ID', 'Individual_ID'], right_on=['FID', 'IID'])

# Selecting the first 230 columns from the merged DataFrame
selected_columns = list(merged_df.columns[:230])

# Adding the specific columns requested
additional_columns = ['RG_ln_r', 'T2D_r', 'age', 'sex', 'BMI', 'T2D', 'RG', 'PC1', 'PC2', 'PC3', 'PC4', 'PC5', 'PC6', 'RG_ln']

# Ensure we are selecting the correct columns from the merged DataFrame
final_columns = selected_columns + additional_columns

# Creating the final DataFrame with the selected columns
final_df = merged_df[final_columns]
final_df.to_csv("merged_data_with_phenotypes.csv", index = False)
