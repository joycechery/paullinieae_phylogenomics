#!/bin/bash

#SBATCH --job-name=01.hybpiper_batch1.sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=6
#SBATCH --mem=128GB
#SBATCH --time=30:00:00
#SBATCH --output=job.%A_%a.out
#SBATCH --error=job.%A_%a.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jgo5750@nyu.edu
#SBATCH --array=0-799

#Script written by Dr.Isaac H. Lichter Marck and modified by Dr. Joyce G. Onyenedum

# Get the SLURM_ARRAY_TASK_ID
index=${SLURM_ARRAY_TASK_ID}

# Read the corresponding line (based on SLURM_ARRAY_TASK_ID) from namelist_batch1.txt
namelist=$(sed -n "$((index + 1))p" /scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/scripts/namelist/namelist_batch1.txt)

# Define input and output directories
INPUT_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/00.trimmomatic/clean_data/batch1"
HYB_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/01.hybpiper"
OUTPUT_DIR="/vast/jgo5750/phylogeny_pau/hybpiper_out"
EXPECTED_OUTPUT="${OUTPUT_DIR}/${namelist}_assembly"

# Check if the output for the current sample already exists
if [[ -d "$EXPECTED_OUTPUT" ]]; then
  echo "Output for sample $namelist already exists. Skipping..."
  exit 0
fi

# Print debug information
echo "Processing sample: $namelist"
echo "Current working directory: $(pwd)"

# Run HybPiper assemble
singularity exec --overlay /vast/jgo5750/cs7566/singularity/overlay-50G-10M.ext3:ro \
/scratch/work/public/singularity/ubuntu-20.04.4.sif \
/bin/bash -c "source /ext3/env.sh && conda activate hybpiper && \
hybpiper assemble -t_dna $HYB_DIR/Paullinieae_140s_352g_target_fixed.fasta \
-r $INPUT_DIR/${namelist}_R*_clean.fastq.gz --bwa --cpu 32 --cov_cutoff 4 \
--hybpiper_output $OUTPUT_DIR"

# Mark the process completion
echo "Sample $namelist processed successfully."
