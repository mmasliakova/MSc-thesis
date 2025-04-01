import pandas as pd
import numpy as np
import os
import csv
import sys
import glob

if len(sys.argv) != 2:
    print("Usage: python3 results_postprocessing.py <TUMOR_TYPE>")
    sys.exit(1)
    
tumor_type = sys.argv[1]

quantum_dir = "/masterthesis_marina/QuantumClone_TCGA"
vcf_path = os.path.join('/masterthesis_marina/data/GDC/snv', tumor_type)

tumor_path = os.path.join(quantum_dir, tumor_type)
samples = [
    dir for dir in os.listdir(tumor_path)
    if os.path.isdir(os.path.join(tumor_path, dir)) and 
    len(os.listdir(os.path.join(tumor_path, dir))) > 2
]

def get_mutation_type(info_column):
    mutation_type = "SNV"
    
    csq_field = [field.split('|') for field in info_column.split(';') if field.startswith('CSQ=')]
    
    if csq_field:
        for csq in csq_field[0]:
            if 'insertion' in csq or 'deletion' in csq:
                mutation_type = "Indel"
    
    return mutation_type


for sample in samples:
    sample_path = os.path.join(tumor_path, sample)
    filtered_path = os.path.join(sample_path, "filtered.csv")
    cluster_path = os.path.join(sample_path, "clustering.csv")
    centers_path = os.path.join(sample_path, "centers.csv")

    filtered = pd.read_csv(filtered_path, sep=',', header=0, index_col=0)
    clusters = pd.read_csv(cluster_path, sep=',', header=0, index_col=0)
    centers = pd.read_csv(centers_path, sep=',', header=0, index_col=0)
    
    clonal_cluster = int(centers["X..i.."].idxmax())
    filtered = filtered.sort_values(by="id")
    filtered['Cluster'] = clusters['Number']
    
    quantum_clonal = filtered[filtered['Cluster'] == clonal_cluster]
    
    quantum_clonal = quantum_clonal.rename(columns={'Chr': '#CHROM', 'Start': 'POS'})
    
    snv_files = os.listdir(vcf_path)
    vcf_files = [f for f in snv_files if f.endswith(".vcf") and f.startswith(sample)]
    
    vcf_file_path = os.path.join(vcf_path, vcf_files[0])
    vcf_data = {}
    
    with open(vcf_file_path, 'r') as vcf_file:
        for line in vcf_file:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            chrom, pos = fields[0], fields[1]
            chrom = chrom[3:]
            info_column = fields[7]
            
            #get the mutation type from the INFO field
            mutation_type = get_mutation_type(info_column)
            
            #store the mutation type by chromosome and position
            vcf_data[(chrom, pos)] = mutation_type
    
    #add mutation type column directly to the quantum_clonal DataFrame
    mutation_types = []
    for _, row in quantum_clonal.iterrows():
        chrom = row['#CHROM']
        pos = str(row['POS'])

        mutation_type = vcf_data.get((chrom, pos), "SNV")
        mutation_types.append(mutation_type)
    
    quantum_clonal['MutationType'] = mutation_types
    
    
    quantumclone_mut = f"{sample_path}/clonalmutations.txt"
    quantum_clonal.to_csv(quantumclone_mut, sep='\t', index=False)
    
    stat = len(quantum_clonal) - 1

    print(f"Processed {sample}: {stat} clonal mutations")
