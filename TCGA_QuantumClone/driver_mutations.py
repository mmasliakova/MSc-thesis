import os
import pandas as pd
from collections import defaultdict
import csv
import sys
import glob

if len(sys.argv) != 2:
    print("Usage: python3 results_postprocessing.py <TUMOR_TYPE>")
    sys.exit(1)
    
tumor_type = sys.argv[1]

quantum_path = "/masterthesis_marina/QuantumClone_TCGA"
vcf_path = "/masterthesis_marina/data/GDC/snv" 

tumor_path = os.path.join(quantum_path, tumor_type)
samples = [
    dir for dir in os.listdir(tumor_path)
    if os.path.isdir(os.path.join(tumor_path, dir)) and 
    len(os.listdir(os.path.join(tumor_path, dir))) > 2
]
vcf_tumor_path = os.path.join(vcf_path, tumor_type)

for sample in samples:
    print(sample)
    prefix = sample[:-8]
    vcf_files = glob.glob(os.path.join(vcf_tumor_path, f"{prefix}*.vcf"))
    vcf_file_path = vcf_files[0]
    with open(vcf_file_path, 'r') as file:
        lines = []
        
        for line in file:
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                header = line.strip().split('\t')
            else:
               lines.append(line.strip().split('\t'))

    vcf = pd.DataFrame(lines, columns=header)
    vcf["POS"] = vcf["POS"].astype(int)
    vcf["#CHROM"] = vcf["#CHROM"].astype(str)
    vcf["#CHROM"] = vcf["#CHROM"].str.replace('chr', '', regex=True)
    
    sample_path = os.path.join(tumor_path, sample, "clonalmutations.txt")
    df = pd.read_csv(sample_path, sep="\t", index_col=False)
    df["POS"] = df["POS"].astype(int)
    df["#CHROM"] = df["#CHROM"].astype(str)
    df["REF"] = None
    df["ALT"] = None
    
    for idx, mutation in df.iterrows():
        chr = mutation["#CHROM"]
        pos = mutation["POS"]
        match = vcf[(vcf["#CHROM"] == chr) & (vcf["POS"] == pos)]
        
        if not match.empty:
            df.at[idx, "REF"] = match.iloc[0]["REF"]
            df.at[idx, "ALT"] = match.iloc[0]["ALT"]

    updated_sample_path = os.path.join(tumor_path, sample, "clonalmutations_with_vcf.txt")
    df.to_csv(updated_sample_path, sep="\t", index=False)

        

mutation_counts = defaultdict(int)

for sample in samples:
    sample_path = os.path.join(tumor_path, sample, "clonalmutations_with_vcf.txt")
    df = pd.read_csv(sample_path, sep="\t", index_col=False)
    unique_mutations = set(zip(df["#CHROM"], df["POS"], df["REF"], df["ALT"]))
    for mutation in unique_mutations:
        mutation_counts[mutation] += 1

mutation_df = pd.DataFrame(list(mutation_counts.items()), columns=["Mutation", "Sample_Count"])
mutation_df[["CHROM", "POS", "REF", "ALT"]] = pd.DataFrame(mutation_df["Mutation"].tolist(), index=mutation_df.index)
mutation_df.drop(columns=["Mutation"], inplace=True)
mutation_df = mutation_df[mutation_df["Sample_Count"] >= 2]
mutation_df = mutation_df.sort_values(by=['Sample_Count'], ascending=False)

mutation_df.to_csv(f"{tumor_path}/common_mutations.tsv", sep="\t", index=False)
