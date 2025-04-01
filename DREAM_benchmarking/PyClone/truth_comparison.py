import pandas as pd
import numpy as np
import os
import glob


pyclone_dir = "/masterthesis_marina/DREAM_benchmarking/pyclone"
truthfiles_dir = "/masterthesis_marina/data/dreams_challenge/truthfiles"
samples = [dir for dir in os.listdir(pyclone_dir) if os.path.isdir(os.path.join(pyclone_dir, dir))]
output_file = "/masterthesis_marina/DREAM_benchmarking/pyclone/statistics_pyclone_allsamples.txt"

with open(output_file, 'w') as out_f:
    for sample in samples:
        sample_id = sample.replace("sample", "")
        print(sample_id)
        
        truth_file_path = os.path.join(truthfiles_dir, f"{sample_id}-noXY", f"{sample_id}_clonal_truth.txt")
        pyclone_res_path = os.path.join(pyclone_dir, f"sample{sample_id}", f"sample{sample_id}_result.tsv")

        pyclone_res = pd.read_csv(pyclone_res_path, sep='\t', header=0, index_col=False)
        max_value = pyclone_res["cellular_prevalence"].max() 
        max_rows = pyclone_res[pyclone_res["cellular_prevalence"] == max_value]
        
        second_max_value = pyclone_res[pyclone_res["cellular_prevalence"] < max_value]["cellular_prevalence"].max()
        second_max_rows = pyclone_res[pyclone_res["cellular_prevalence"] == second_max_value]
        
        #combine both subsets
        pyclone_clones = pd.concat([max_rows, second_max_rows])
        #max_cell = pyclone_res["cellular_prevalence"].max()
        #pyclone_clones = pyclone_res[pyclone_res["cellular_prevalence"] == max_cell]
        
        pyclone_clones[['#CHROM', 'POS']] = pyclone_clones['mutation_id'].str.split(':', n=2, expand=True).iloc[:, :2]
        pyclone_clones['#CHROM'] = pyclone_clones['#CHROM'].astype(int)
        pyclone_clones['POS'] = pyclone_clones['POS'].astype(int)
        pyclone_clones = pyclone_clones.sort_values(by=['#CHROM', 'POS'])
        
        pyclone_mut = f"{pyclone_dir}/sample{sample_id}/{sample_id}_clonalmutations.txt"
        pyclone_clones.to_csv(pyclone_mut, sep='\t', index=False)

        
        df_filtered = pd.read_csv(truth_file_path, sep='\t', header=0, index_col=False)

        #find common mutations
        common_mut = pyclone_clones.merge(df_filtered, on=['#CHROM', 'POS'])
        false_positives = pyclone_clones.merge(df_filtered, on=['#CHROM', 'POS'], how='left', indicator=True)
        false_positives = false_positives[false_positives['_merge'] == 'left_only']
        fp_count = len(false_positives)
        fpr = (fp_count / len(pyclone_clones)) * 100
        stat = len(common_mut) / len(df_filtered) * 100

        out_f.write(f"Percent of identified mutations by Pyclone-VI in {sample_id}: {stat:.2f}%\n")
        out_f.write(f"Percent of false positives in {sample_id}: {fpr:.2f}%\n")
        print(f"Processed {sample_id}: {stat:.2f}%")