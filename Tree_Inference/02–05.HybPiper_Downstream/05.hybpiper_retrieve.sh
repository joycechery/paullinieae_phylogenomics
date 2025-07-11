#!/bin/bash

#SBATCH --job-name=05.hybpiper_retrieve.sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=64GB
#SBATCH --time=4:00:00
#SBATCH --output=job.%A.out
#SBATCH --error=job.%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jgo5750@nyu.edu

#Script written by Dr.Isaac H. Lichter Marck and modified by Dr. Joyce G. Onyenedum

#set directory for sequence files
TARGET="/vast/jgo5750/phylogeny_pau/hybpiper_out"

#set home directory
HYB_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/01.hybpiper/"

#namelist
NAMELIST="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/scripts/namelist"

#set directory for tree inference                                                                                                           
                                                                
ALIGN_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/02.alignments/pau_332"

mkdir $ALIGN_DIR/exons_pau
mkdir $ALIGN_DIR/supercontigs_pau
cd $TARGET

singularity exec --overlay /vast/jgo5750/cs7566/singularity/overlay-50G-10M.ext3:ro \
        /scratch/work/public/singularity/ubuntu-20.04.4.sif \
        /bin/bash -c "source /ext3/env.sh; conda activate hybpiper; hybpiper retrieve_sequences dna -t_dna $HYB_DIR/Paullinieae_140s_352g_target_fixed.fasta --sample_names $NAMELIST/namelist_all.txt"

singularity exec --overlay /vast/jgo5750/cs7566/singularity/overlay-50G-10M.ext3:ro \
        /scratch/work/public/singularity/ubuntu-20.04.4.sif \
        /bin/bash -c "source /ext3/env.sh; conda activate hybpiper; hybpiper retrieve_sequences supercontig -t_dna $HYB_DIR/Paullinieae_140s_352g_target_fixed.fasta --sample_names $NAMELIST/namelist_all.txt"

mv *.FNA $ALIGN_DIR/exons_pau/
mv *.fasta $ALIGN_DIR/supercontigs_pau/
