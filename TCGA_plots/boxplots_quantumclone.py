import pandas as pd
import numpy as np
import os
import csv
import glob
import re
import matplotlib.pyplot as plt
import seaborn as sns

stat_dir = "/masterthesis_marina"

tumortype_pattern = re.compile(r"statistics_fastclone-quantum_TCGA_(\w+)\.txt")
mutation_pattern = re.compile(r"Number of QuantumClone clonal mutations: (\d+\.\d+)")

data = {}

for filename in os.listdir(stat_dir):
    match = tumortype_pattern.match(filename)
    if match:
        tumor_type = match.group(1)
        file_path = os.path.join(stat_dir, filename)

        with open(file_path, "r") as file:
            mutations = [float(m.group(1)) for line in file if (m := mutation_pattern.search(line))]
            
            if mutations:
                data[tumor_type] = mutations

top_30_tumors = sorted(data.keys())[:30]
data = {t: data[t] for t in top_30_tumors}

mutation_df = pd.DataFrame([(t, m) for t, ms in data.items() for m in ms], columns=["Tumor Type", "Mutations"])

plt.figure(figsize=(12, 6))
sns.boxplot(x="Tumor Type", y="Mutations", data=mutation_df)
plt.xticks(rotation=45, ha="right")
plt.xlabel("Tumor Type")
plt.ylabel("Number of QuantumClone Clonal Mutations")
plt.title("Distribution of QuantumClone Clonal Mutations Across Tumor Types")
plt.grid(axis="y", linestyle="--", alpha=0.7)


plot_path = "QuantumClone_boxplots.png"
plt.savefig(plot_path)
plt.close()

print(f"Plot saved as {plot_path}")


