#!/bin/bash

#SBATCH --job-name=02.hybpiper_stats.sh 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64GB
#SBATCH --time=10:00:00
#SBATCH --output=job.%A.out
#SBATCH --error=job.%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jgo5750@nyu.edu

#Script written by Dr.Isaac H. Lichter Marck and modified by Dr. Joyce G. Onyenedum

#set home directory
HYB_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/01.hybpiper"
NAME_LIST="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/scripts/namelist"
#set directory for sequence files
TARGET="/vast/jgo5750/phylogeny_pau/hybpiper_out"

cd $TARGET

# Run hybpiper stats
singularity exec --overlay /vast/jgo5750/cs7566/singularity/overlay-50G-10M.ext3:ro \
    /scratch/work/public/singularity/ubuntu-20.04.4.sif \
    /bin/bash -c "source /ext3/env.sh; conda activate hybpiper; \
    hybpiper stats -t_dna $HYB_DIR/Paullinieae_140s_352g_target_fixed.fasta gene $NAME_LIST/namelist_all.txt"

