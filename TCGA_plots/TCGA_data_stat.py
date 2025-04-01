import pandas as pd
import numpy as np
import os
import csv
import sys
import glob
import matplotlib.pyplot as plt
import seaborn as sns

samples_base_path = '/masterthesis_marina/QuantumClone_TCGA'
vcf_base_path = '/masterthesis_marina/data/GDC/snv'

mutation_counts = {}

def get_mutation_type(info_column):
    mutation_type = "SNV"
    
    csq_field = [field.split('|') for field in info_column.split(';') if field.startswith('CSQ=')]
    
    if csq_field:
        for csq in csq_field[0]:
            if 'insertion' in csq or 'deletion' in csq:
                mutation_type = "Indel"
    
    return mutation_type

                
for tumor in os.listdir(vcf_base_path):
    tumor_path = os.path.join(vcf_base_path, tumor)
    if os.path.isdir(tumor_path) and tumor != 'metadata':
        samples_path = os.path.join(samples_base_path, tumor)
        
        total_snv = 0
        total_indel = 0
        sample_count = 0

        for sample in os.listdir(samples_path):
            sample_path = os.path.join(samples_path, sample)
            if not os.path.isdir(sample_path):
                continue

            snv_files = os.listdir(tumor_path)
            vcf_files = [f for f in snv_files if f.endswith(".vcf") and f.startswith(sample)]

            vcf_file_path = os.path.join(vcf_base_path, tumor, vcf_files[0])
            snv_count = 0
            indel_count = 0
            
            with open(vcf_file_path, 'r') as vcf_file:
                for line in vcf_file:
                    if line.startswith('#'):
                        continue
                    fields = line.strip().split('\t')
                    info_column = fields[7]
                    mutation_type = get_mutation_type(info_column)

                    if mutation_type == "SNV":
                        snv_count += 1
                    elif mutation_type == "Indel":
                        indel_count += 1

            total_snv += snv_count
            total_indel += indel_count
            sample_count += 1

        if sample_count > 0:
            avg_snv = (total_snv / (total_snv + total_indel)) * 100 if (total_snv + total_indel) > 0 else 0
            avg_indel = (total_indel / (total_snv + total_indel)) * 100 if (total_snv + total_indel) > 0 else 0
            mutation_counts[tumor] = {"SNV": avg_snv, "Indel": avg_indel}


mutation_df = pd.DataFrame(mutation_counts).T
mutation_df.index.name = "Tumor Type"


mutation_df.sort_index().plot(kind="bar", stacked=True, color=["#1f77b4", "#ffcc00"], figsize=(12, 6))

plt.ylabel("Average Percentage (%)")
plt.xlabel("Tumor Type")
plt.title("Average Percentage of SNVs and Indels per Tumor Type")
plt.xticks(rotation=45, ha="right")
plt.legend(title="Mutation Type")
stacked_plot_path = "TCGA_data_Indel_SNV_barplot.png"
plt.savefig(stacked_plot_path)
plt.close()
