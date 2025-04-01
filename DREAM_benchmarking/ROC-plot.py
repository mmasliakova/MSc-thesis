import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

pyclone_file = "/masterthesis_marina/DREAM_benchmarking/pyclone/statistics_pyclone_partsamples.txt"
phylo_file = "/masterthesis_marina/DREAM_benchmarking/phylowgs/statistics.txt"
quantum_file = "/masterthesis_marina/DREAM_benchmarking/QuantumClone/statistics_quantumclone_partsamples.txt"
fastclone_file = "/masterthesis_marina/DREAM_benchmarking/fastclone/statistics_fastclone_partsamples.txt"
quantum_file2 = "/masterthesis_marina/DREAM_benchmarking/QuantumClone_run2/statistics_quantumclone_partsamples.txt"


def load_data(file_path, method_name):
    df = pd.read_csv(file_path, sep=":", header=None, names=["desc", "percent"])
    df["sample_id"] = df["desc"].str.extract(r'in (\S+)')
    df["percent"] = df["percent"].str.replace("%", "", regex=False).astype(float)
    
    false_positive_df = df[df["desc"].str.contains("false positives", na=False)]
    false_positive_df = false_positive_df.rename(columns={"percent": "false_positive_rate", "desc": "desc_false_positive"})
    
    mutation_df = df[df["desc"].str.contains("identified mutations", na=False)]
    mutation_df = mutation_df.rename(columns={"percent": "identified_mutation_rate", "desc": "desc_mutation"})
    
    merged_df = false_positive_df.merge(mutation_df, on="sample_id", suffixes=("", ""))
    
    merged_df["method"] = method_name
    
    return merged_df[["sample_id", "false_positive_rate", "identified_mutation_rate", "method"]]
    

pyclone_df = load_data(pyclone_file, "PyClone-VI")
phylowgs_df = load_data(phylo_file, "PhyloWGS")
quantum_df = load_data(quantum_file, "QuantumClone")
fastclone_df = load_data(fastclone_file, "FastClone")
quantum_df2 = load_data(quantum_file2, "QuantumClone_Run2")

#combine all methods into one DataFrame
merged_df = pd.concat([pyclone_df, phylowgs_df, quantum_df, fastclone_df, quantum_df2])

summarized_df = merged_df.groupby('method').agg(
    mean_TPR=('identified_mutation_rate', 'mean'),
    std_TPR=('identified_mutation_rate', 'std'),
    mean_FPR=('false_positive_rate', 'mean'),
    std_FPR=('false_positive_rate', 'std')
).reset_index()

tool_colors = {
    'PyClone-VI': 'pink',
    'PhyloWGS': 'darkblue',
    'QuantumClone': 'brown',
    'FastClone': 'orange',
    'QuantumClone_Run2': 'grey',
}

plt.figure(figsize=(10, 6))
sns.scatterplot(x='mean_FPR', y='mean_TPR', hue='method', data=summarized_df, 
              s=100, palette=["orange", "darkblue", "pink", "brown", "grey"], edgecolor="black")
              
for idx, row in summarized_df.iterrows():
    plt.errorbar(
        row['mean_FPR'], 
        row['mean_TPR'], 
        xerr=row['std_FPR'], 
        yerr=row['std_TPR'], 
        fmt='o', 
        linestyle='none', 
        color=tool_colors[row['method']],
        capsize=3
    )
    
plt.xlim(-10, 100) 
plt.ylim(0, 110)

plt.title('Summarized True Positive Rate vs False Positive Rate', fontsize=14)
plt.xlabel('1 - Specificity (False Positive Rate)', fontsize=12)
plt.ylabel('Recall (True Positive Rate)', fontsize=12)
plt.legend(title="Method", bbox_to_anchor=(1.05, 1), loc='upper left')
plt.grid(True)

plt.savefig(f"summary_plot.png", dpi=300, bbox_inches="tight")
plt.close()