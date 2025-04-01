import glob
import sys
import os

tumor_type = sys.argv[1]

base_path = "masterthesis_marina/QuantumClone_TCGA"
samples_path = os.path.join(base_path, tumor_type, "*")

jobs_dir = os.path.join(base_path, f"{tumor_type}_jobs")
os.makedirs(jobs_dir, exist_ok=True)

sample_folders = [f for f in glob.glob(samples_path)]

for sample in sample_folders:
    sample_name = os.path.basename(sample)
    
    sample_dir = os.path.join(base_path, tumor_type, sample_name)

    file_snv = glob.glob(os.path.join(sample_dir, "*_SNVlist.txt"))
    file_freec = glob.glob(os.path.join(sample_dir, "*_freec.txt"))

    if not file_snv or not file_freec:
        print(f"Skipping {sample_name} - Required files not found.")
        continue
    
    job_script = os.path.join(jobs_dir, f"job_{sample_name}.sh")

    print(f"Generating job script for {sample_name} in {jobs_dir}:")
    with open(job_script, "w") as fh:
        fh.write(f"""#!/bin/sh

#PBS -N QC_{sample_name}
#PBS -l nodes=1:ppn=12
#PBS -l walltime=3:00:00
#PBS -l mem=200gb
#PBS -m abe
#PBS -o {jobs_dir}/RunClustering_{sample_name}.out
#PBS -e {jobs_dir}/RunClustering_{sample_name}.err
#PBS -d {base_path}/{tumor_type}

ml R/4.4.1-gfbf-2023b
module load R-bundle-CRAN/2024.06-foss-2023b

cd {base_path}/{tumor_type}

Rscript "{base_path}/clustering.R" {tumor_type} {sample_name}
""")

    os.system(f"qsub {job_script}")

print(f"All jobs for tumor type {tumor_type} submitted.")