#!/bin/bash

#SBATCH --job-name=07.trimal_exons.sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64GB
#SBATCH --time=3:00:00
#SBATCH --output=job.%A.out
#SBATCH --error=job.%A.err
#SBATCH --mail-type=END
#SBATCH --mail-user=jgo5750@nyu.edu

#Script written by Dr.Isaac H. Lichter Marck and modified by Dr. Joyce G. Onyenedum

#set directory for tree inference

HOME_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/02.alignments/"

ALIGN_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/02.alignments/pau_333"

#set directory for exon fasta files
EXONS_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/02.alignments/pau_333/aligned_exons_pau"

#Step 7: trim alignments to remove trash using the trimal utility
module load trimal/intel/1.4.1

cd $ALIGN_DIR
mkdir -p trimmed_exons_pau  # Create directory if it doesn't exist

while read name; do
    trimal -in $EXONS_DIR/aligned_${name}.fasta -out $ALIGN_DIR/trimmed_exons_pau/clean_${name}.fasta -gt 0.85 -cons 60 
done < $HOME_DIR/namelist_regions_order.txt
