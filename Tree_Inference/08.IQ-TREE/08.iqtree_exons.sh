#!/bin/bash

#SBATCH --job-name=08.iqtree_exons_pau333s_351g.sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=8
#SBATCH --mem=64GB
#SBATCH --time=100:00:00
#SBATCH --output=job.%J.out
#SBATCH --error=job.%J.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jgo5750@nyu.edu

#Script written by Dr.Isaac H. Lichter Marck and modified by Dr. Joyce G. Onyenedum

#set directory for tree inference 
TREE_DIR="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/03.tree_inference/iqtree"

#Step 7: partitioned iqtree ML tree inference
module purge
module load iqtree/1.6.12  

iqtree -s $TREE_DIR/pau_333s_351g.fasta -m GTR -bb 1000 -spp $TREE_DIR/pau_333s_351g_partitions

echo = `date` job $JOB_NAME done

