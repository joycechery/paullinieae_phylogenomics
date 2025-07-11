#!/bin/bash                                                                                                                              

#SBATCH --job-name=09.astral_exons_pau333s_351g.sh
#SBATCH --nodes=1 
#SBATCH --ntasks-per-node=8               
#SBATCH --mem=32GB
#SBATCH --time=100:00:00  
#SBATCH --output=job.%A.out
#SBATCH --error=job.%A.err
#SBATCH --mail-type=BEGIN,END
#SBATCH --mail-user=jgo5750@nyu.edu                                                                                      

#Script written by Dr.Isaac H. Lichter Marck and modified by Dr. Joyce G. Onyenedum

HOME="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/02.alignments"
EXONS="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/02.alignments/pau_333/trimmed_exons_pau/"
ASTRAL="/vast/jgo5750/phylogeny_pau/astral"

module load iqtree/1.6.12

cd $ASTRAL

# Process each name from the list once
while read name; do
    iqtree -s $EXONS/clean_${name}.fasta -m GTR -bb 1000 -nt 8 -mem 32G -pre $ASTRAL/clean_${name}

    # Move all output files related to this run into pau_333_test
    mv $ASTRAL/clean_${name}* $ASTRAL/pau_333/

done < $HOME/namelist_regions_order.txt

echo "$(date) job $JOB_NAME done"
