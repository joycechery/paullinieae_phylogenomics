#!/bin/bash

#SBATCH --job-name=06.mafft_exons_pau333.sh 
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64GB
#SBATCH --time=10:00:00
#SBATCH --output=job.%A.out
#SBATCH --error=job.%A.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jgo5750@nyu.edu

#Script written by Dr.Isaac H. Lichter Marck and modified by Dr. Joyce G. Onyenedum

#set directory for tree inference 

HOME_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/02.alignments/"
ALIGN_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/02.alignments/pau_333"

#set directory for exon fasta files
EXONS_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/02.alignments/pau_333/exons_pau"

#Step 6: align matrices using MAFFT
module load mafft/intel/7.475 

mkdir -p $ALIGN_DIR/aligned_exons_pau  # Create directory if it doesn't exist

while read name; do
    echo "Aligning ${name}.fasta"
    mafft --thread -8 --auto $EXONS_DIR/"${name}.FNA" > $ALIGN_DIR/aligned_exons_pau/aligned_"${name}".fasta
done < $HOME_DIR/namelist_regions_order.txt

