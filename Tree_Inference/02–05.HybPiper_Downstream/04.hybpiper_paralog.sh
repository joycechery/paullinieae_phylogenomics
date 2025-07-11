#!/bin/bash

#SBATCH --job-name=04.hybpiper_paralog.sh
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --mem=30GB
#SBATCH --time=2:00:00
#SBATCH --output=job.%A.out
#SBATCH --error=job.%A.err
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --mail-user=jgo5750@nyu.edu

#Script written by Dr.Isaac H. Lichter Marck and modified by Dr. Joyce G. Onyenedum

module load python/intel/3.8.6

#set home directory
NAME_LIST="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/scripts/namelist"
SCRIPTS="/scratch/projects/onyenedumlab/Joyce_Onyenedum/paullinieae_phylogeny/scripts/"

#set directory for sequence files
TARGET="/vast/jgo5750/phylogeny_pau/hybpiper_out"


# Execute Python script to filter samples

cd $TARGET

while read name
do grep -o '....$' ${name}/${name}_genes_with_long_paralog_warnings.txt > ${name}/genes_with_paralog_warnings.txt
sed "s/[^,]*,\([^,]*\),.*/\1/" ${name}/${name}_genes_derived_from_putative_chimeric_stitched_contig.csv >> ${name}/genes_with_paralog_warnings.txt

while read directory

do rm -r  ${name}/${directory}/${name}/sequences

done < ${name}/genes_with_paralog_warnings.txt

done < $NAME_LIST/namelist_all.txt
