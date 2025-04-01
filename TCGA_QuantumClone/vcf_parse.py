import pandas as pd
import glob
import os
import csv
import sys
import random

if len(sys.argv) != 2:
    print("Usage: python3 vcf_parse.py <TUMOR_TYPE>")
    sys.exit(1)

tumor_type = sys.argv[1]

snv_data_path = "/masterthesis_marina/data/GDC/snv" 
cnv_data_path = "/masterthesis_marina/data/GDC/cnv"
sample_list_file = f"/masterthesis_marina/data/GDC/{tumor_type}_matched.txt"
base_output_dir = "/masterthesis_marina/QuantumClone_TCGA"

tumor_path = os.path.join(base_output_dir, tumor_type)

if not os.path.exists(tumor_path):
        os.makedirs(tumor_path)


tumor_cnv = os.path.join(cnv_data_path, tumor_type)
tumor_snv = os.path.join(snv_data_path, tumor_type)

if not os.path.isfile(sample_list_file):
    print(f"Error: Sample list file not found: {sample_list_file}")
    sys.exit(1)

if not os.path.isdir(tumor_cnv) or not os.path.isdir(tumor_snv):
    print(f"Error: Tumor folder not found: {tumor_type}")
    sys.exit(1)

with open(sample_list_file, "r") as f:
    sample_list = [line.strip() for line in f]
    
sample_groups = {}
for sample in sample_list:
    prefix = sample[:-8]
    if prefix not in sample_groups:
        sample_groups[prefix] = []
    sample_groups[prefix].append(sample)

paired_samples = [samples for samples in sample_groups.values() if len(samples) == 2]

if len(paired_samples) > 50:
    paired_samples = random.sample(paired_samples, 50)
print(paired_samples)

final_sample_list = [sample for pair in paired_samples for sample in pair]

samples = set(final_sample_list)
print(samples)

def get_sample_id(filename):
    parts = filename.split('.')
    return parts[0]

cnv_files = os.listdir(tumor_cnv)
snv_files = os.listdir(tumor_snv)

cnv_files = [f for f in cnv_files if get_sample_id(f) in samples]
snv_files = [f for f in snv_files if f.endswith(".vcf") and get_sample_id(f) in samples]

def parse_format_field(row, field_name):
    format_fields = row['FORMAT'].split(':')
    tumor_values = row['TUMOR'].split(':')
    field_index = format_fields.index(field_name) 
    return tumor_values[field_index]
    
    
def remove_zeroCN(cnv_file):
    cnv_df = pd.read_csv(cnv_file, sep='\t')
    cn0_regions = cnv_df[cnv_df["Copy_Number"] == 0]
    return cn0_regions


def convert_cnv_to_freec(cnv_file):
    
    header = ["Chromosome", "Start", "Ratio", "MedianRatio", "CopyNumber", "BAF", "EstimatedBAF", "Genotype", "UncertaintyOfGT"]

    cnv_df = pd.read_csv(cnv_file, sep='\t')
    
    cn0_regions = remove_zeroCN(cnv_file)
    filtered_cnv = cnv_df[~cnv_df.apply(
        lambda row: any(
            (row["Chromosome"] == cn_row["Chromosome"]) and
            (cn_row["Start"] <= row["Start"] <= cn_row["End"])
            for _, cn_row in cn0_regions.iterrows()
        ),
        axis=1
    )]
    

    freec = []
    for _, row in filtered_cnv.iterrows():
        ratio = int(row['Copy_Number']) / 2

        BAF = int(row['Minor_Copy_Number']) / int(row['Copy_Number'])
        corr_BAF = abs(round(BAF, 2) - 0.5)

        #assign genotype
        if int(row['Major_Copy_Number']) == 0 and int(row['Minor_Copy_Number']) == 0:
            genotype = "0"  # Homozygous deletion
        elif int(row['Major_Copy_Number']) == 1 and int(row['Minor_Copy_Number']) == 0:
            genotype = "A"  # Hemizygous deletion
        else:
            genotype = "A" * int(row['Major_Copy_Number']) + "B" * int(row['Minor_Copy_Number'])

        freec.append([row['Chromosome'], row['Start'], ratio, ratio, row['Copy_Number'], corr_BAF, round(BAF, 2), genotype, 0])
        
        freec_file = pd.DataFrame(freec, columns=header)
        freec_file['Chromosome'] = freec_file['Chromosome'].str.replace('chr', '', regex=True)
    
    return freec_file


cn0_dict = {}
for filename in cnv_files:
    cn0_dict[filename] = remove_zeroCN(os.path.join(tumor_cnv, filename))


for filename in snv_files:
    sample_id = get_sample_id(filename)[:-8]
    print(f"Processing SNV: {filename}")
    
    cn0_regions = cn0_dict.get(filename.replace('.vcf', '.cnv'), pd.DataFrame())
    
    with open(os.path.join(tumor_snv, filename), 'r') as file:
        lines = []
        
        for line in file:
            if line.startswith('##'):
                continue
            elif line.startswith('#'):
                header = line.strip().split('\t')
            else:
               lines.append(line.strip().split('\t'))

    vcf = pd.DataFrame(lines, columns=header)
    
    vcf['POS'] = pd.to_numeric(vcf['POS'])
    drop_indices = []
    
    for index, record in vcf.iterrows():
        chrom = record['#CHROM']
        pos = record['POS']

        #check if the SNV falls in any CN=0 region
        in_cn0_region = any(
            (chrom == row["Chromosome"]) and (row["Start"] <= pos <= row["End"])
            for _, row in cn0_regions.iterrows()
        )

        #if it's in CN=0 regions, keep it
        if in_cn0_region:
            drop_indices.append(index)
        
    vcf.drop(index=drop_indices, inplace=True)
        
    filtered_vcf = vcf[vcf['FILTER'] == 'PASS']
    
    filtered_vcf = vcf[vcf['FILTER'] == 'PASS']

    if filtered_vcf.empty:
        print(f"Warning: Skipping {sample_id}.")
        continue
    
    output = pd.DataFrame({
        'Sample': sample_id,
        'Chr': filtered_vcf['#CHROM'],
        'Start': filtered_vcf['POS']
    })

    
    output['Depth'] = filtered_vcf.apply(lambda row: parse_format_field(row, 'DP'), axis=1)
    output['Alt'] = filtered_vcf.apply(lambda row: parse_format_field(row, 'AD').split(',')[1], axis=1)
    
    output['Depth'] = output['Depth'].astype(int)
    output['Alt'] = output['Alt'].astype(int)
    
    output['Chr'] = output['Chr'].str.replace('chr', '', regex=True)
    
    output_dir = os.path.join(tumor_path, f"{sample_id}")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    output_filename = os.path.join(output_dir, f"{sample_id}_SNVlist.txt")
    
    output.to_csv(output_filename, sep='\t', index=False)


for filename in cnv_files:
    
    sample_id = get_sample_id(filename)[:-8]
    print(f"Processing CNV: {filename}")

    freec = convert_cnv_to_freec(os.path.join(tumor_cnv, filename))
    
    output_dir = os.path.join(tumor_path, f"{sample_id}")
    
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)
    
    freec_filename = os.path.join(output_dir, f"{sample_id}_freec.txt")
    
    freec.to_csv(freec_filename, sep='\t', index=False)



