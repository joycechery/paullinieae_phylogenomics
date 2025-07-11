#!/bin/bash

#SBATCH --job-name=03.hybpiper_heatmap.sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64GB
#SBATCH --time=5:00:00
#SBATCH --output=job.%A.out
#SBATCH --error=job.%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jgo5750@nyu.edu

#Script written by Dr.Isaac H. Lichter Marck and modified by Dr. Joyce G. Onyenedum

#set home directory
HYB_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/01.hybpiper/"

#set directory for sequence files
TARGET="/vast/jgo5750/phylogeny_pau/hybpiper_out"

cd $HYB_DIR

# Run hybpiper recovery_heatmap
singularity exec --overlay /vast/jgo5750/cs7566/singularity/overlay-50G-10M.ext3:ro \
    /scratch/work/public/singularity/ubuntu-20.04.4.sif \
    /bin/bash -c "source /ext3/env.sh; conda activate hybpiper; \
    hybpiper recovery_heatmap $TARGET/seq_lengths.tsv --heatmap_dpi 200"
