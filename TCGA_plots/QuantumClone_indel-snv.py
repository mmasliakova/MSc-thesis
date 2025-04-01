import os
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

stat_dir = "/masterthesis_marina"

tumortype_pattern = re.compile(r"statistics_fastclone-quantum_TCGA_(\w+)\.txt")
indel_pattern = re.compile(r"Percentage of QuantumClone Indel clonal mutations: (\d+\.\d+)%")
snv_pattern = re.compile(r"Percentage of QuantumClone SNV clonal mutations: (\d+\.\d+)%")

indel_snv_data = {}

for filename in os.listdir(stat_dir):
    match = tumortype_pattern.match(filename)
    if match:
        tumor_type = match.group(1)
        file_path = os.path.join(stat_dir, filename)

        indels, snvs = [], []
        with open(file_path, "r") as file:
            for line in file:
                if (i := indel_pattern.search(line)):
                    indels.append(float(i.group(1)))
                if (s := snv_pattern.search(line)):
                    snvs.append(float(s.group(1)))

        if indels and snvs:
            indel_snv_data[tumor_type] = {
                "SNV": np.median(snvs), 
                "Indel": np.median(indels) 
            }

indel_snv_df = pd.DataFrame.from_dict(indel_snv_data, orient="index").reset_index().rename(columns={"index": "Tumor Type"})

indel_snv_df.set_index("Tumor Type").sort_index().plot(kind="bar", stacked=True, color=["#1f77b4", "#ffcc00"], figsize=(12, 6))

plt.xticks(rotation=45, ha="right")
plt.xlabel("Tumor Type")
plt.ylabel("Average Percentage")
plt.title("Average Percentage of Indel vs. SNV in QuantumClone Mutations")
plt.legend(title="Mutation Type")
plt.grid(axis="y", linestyle="--", alpha=0.7)

stacked_plot_path = "QuantumClone_Indel_SNV_barplot.png"
plt.savefig(stacked_plot_path)
plt.close()