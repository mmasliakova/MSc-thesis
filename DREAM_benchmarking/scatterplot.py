import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

pyclone_file = "/masterthesis_marina/DREAM_benchmarking/pyclone/statistics_pyclone_partsamples.txt"
phylo_file = "/masterthesis_marina/DREAM_benchmarking/phylowgs/statistics.txt"
quantum_file = "/masterthesis_marina/DREAM_benchmarking/QuantumClone/statistics_quantumclone_part.txt"
fastclone_file = "/masterthesis_marina/DREAM_benchmarking/fastclone/statistics_fastclone_part.txt"

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

#combine all methods into one DataFrame
merged_df = pd.concat([pyclone_df, phylowgs_df, quantum_df, fastclone_df])

plt.figure(figsize=(10, 6))
sns.scatterplot(
    data=merged_df,
    x="identified_mutation_rate",
    y="false_positive_rate",
    hue="method",
    palette=["orange", "darkblue", "pink", "brown"],
    s=100,
    edgecolor="black"
)

plt.xlabel("Identified Mutations (%)")
plt.ylabel("False Positive Rate (%)")
plt.title("False Positive Rate vs Identified Mutations")
plt.grid(True)
plt.legend(title="Method")

plt.savefig("scatterplot.png", dpi=300, bbox_inches="tight")
plt.close()