#!/bin/bash
#SBATCH --job-name=cellrangercount
#SBATCH --output=cellranger-count%j.out
#SBATCH --error=cellranger-count%j.err
#SBATCH  -p cpu-nodes
#SBATCH  -N 1
#SBATCH  -n 16
#SBATCH --mem=64G
#SBATCH --time=4:00:00
#SBATCH --mail-type=ALL
#SBATCH --mail-user=your-jphiri@mlw.mw
#SBATCH --account=laboratory
#SBATCH --array=0-11 # 0-11 for 12 samples

# Load modules
module load cellranger

# Define sample names from your FAST directory
SAMPLES=("CFU12Y-GE" "CFU130-GE" "CFU131-GE" "CFU134-GE" "CUF135-GE" "CUF136-GE" \
"CUF137-GE" "CUF13I-GE" "CUF13J-GE" "CUF13K-GE" "CUG11X" "CUH124")
SAMPLE=${SAMPLES[$SLURM_ARRAY_TASK_ID]}

# Run  Cell Ranger
cellranger count --id=${SAMPLE} \
                 --transcriptome=/head/NFS/jphiri/refdata-gex-GRCh38-2024-A \
                 --fastqs=/head/NFS/jphiri/fastq/ \
                 --sample=${SAMPLE} \
                 --localcores=8 \
                 --localmem=64 \
--create-bam=true

mkdir -p /head/NFS/jphiri/data/raw
mv ${SAMPLE}/outs/raw_feature_bc_matrix /head/NFS/jphiri/data/raw/${SAMPLE}

