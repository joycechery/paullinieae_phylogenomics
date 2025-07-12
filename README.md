This README document describes the bioinformatic workflow for Cunha Neto et al. “Rampant convergent evolution of vascular oddities and a synnovation characterize the rapid radiation of Paullinieae lianas”.
MAIN FOLDERS.
-	Tree Inference. This folder houses subfolders accounting for the full tree inference workflow: cleaning raw data--> HybPiper Assembly-->HybPiper Downstream Analyses --> Sequence Alignments --> Maximum Likelihood Tree Building -->Chronogram. 
-	Phylogenetic Comparative Methods. This folder houses subfolders accounting for all phylogenetic comparative methods analyses: tests of correlated evolution, state-dependent diversification analyses, ancestral state reconstructions of vascular variant evolution. 
General notes about the content of each folder
-	Files ending in “.sh” are bash scripts
-	Bash scripts are ordered from 00–09, indicating the order in which each script was performed throughout the bioinformatic pipeline
-	Files ending in “.txt” are tab delimited text files
-	Files ending in “.fa” and “.FNA” are fasta files
-	Files ending in “.tre” are phylogenetic tree files
See below for a description of the contents and details of the workflow for each of the two main folders.
Tree Inference Subfolders 
-	00–01.Batch1_Trimmomatic&HybPiper_Assembly
-	00–01.Batch2_Trimmomatic&HybPiper_Assembly
-	00–01.Batch3_Trimmomatic&HybPiper_Assembly
Description: Given the computational constraints on the NYU High Performance Computer (HPC), the cleaning of raw data and the subsequent HybPiper assemblies of those cleaned data were processed separated into three batches. These three batches are contained in the following subfolders within the main folder “Tree Inference”.
-	02–05.HybPiper_Downstream
Description: This folder contains scripts and results to evaluate the HybPiper assembly, including creating a heatmap, sequence statistics, removing paralogs for downstream analysis, and the retrieval of sequences to generate alignments for tree inference.
-	06–07.Sequence_Alignments
Description: This folder contains scripts to align and trim exons. Also within this subfolder is the concatenated alignment that is used for the iqtree analysis (pau_333s_351g.fasta).

-	08.IQtree
Description: This folder contains the script to run a partitioned maximum likelihood tree inference using IQ-TREE v. 1.6.12 on the mafft & post-trimal exon concatenated alignment. The iqtree results are presented in a subfolder named “iqtree_results”. A script to change the tree tip names from sample ID (e.g., “PAU24”) to species names (“Serjania_atrolineata_b”) is available (ChangeSampleNames_iqtree.sh)
-	09.ASTRAL
Description: This folder contains the script to run individual gene trees for each of the 351 mafft + trimal exon alignments. The subfolders content: “trimmed_exons” contains 351 post mafft +trimal alignments; “exon_trees” contain 351 gene trees; “ASTRAL_analysis” contains the input trees and result of an ASTRAL analysis. 
-	10.CHRONOGRAM
Given the size of the dataset (351 genes and 333 samples) a BEAST analysis was not feasible, so we opted to use treePL. Here we used the iqtree consensus tree as a topological constraint to make 1000 bootstrap trees with varying branch lengths. Then we subsampled 100 random bootstrap trees to individually time calibrate using treePL. Those 100 chronograms were collapsed using TreeAnnotator to generate a chronogram where the height 95%HPD is equivalent to the uncertainty given the variability of branch lengths among the 100 bootstrap chronogram trees. 


00–01.Batch1_Trimmomatic&HybPiper_Assembly
Scripts:
00.trimmomatic_batch1.sh
-	This bash script used trimmomatic 0.39 to remove adapters and low-quality sequences from the raw data in batch 1.
01.hybpiper_batch1.sh
-	This bash script used HybPiper to generate consensus assemblies of Angiosperms353 loci data in batch 1.
Files:
namelist_batch1.txt
-	this is the list of sample names within batch 1.
TruSeq3-PE.fa
-	 This is the TruSeq adapter fasta file necessary to run 00.trimmomatic_batch1.sh.


00–01.Batch2_Trimmomatic&HybPiper_Assembly
Scripts:
00.trimmomatic_batch2.sh
-	This bash script used trimmomatic 0.39 to remove adapters and low-quality sequences from the raw data in batch 2.
01.hybpiper_batch2.sh
-	This bash script used HybPiper to generate consensus assemblies of Angiosperms353 loci data in batch 2.
Files:
namelist_batch2.txt
-	This is the list of sample names within batch 2.
TruSeq3-PE.fa
-	 This is the TruSeq adapter fasta file necessary to run 00.trimmomatic_batch2.sh.

00–01.Batch3_Trimmomatic&HybPiper_Assembly
Scripts:
00.trimmomatic_batch3.sh
-	This bash script used trimmomatic 0.39 to remove adapters and low-quality sequences from the raw data in batch 3.
01.hybpiper_batch3.sh
-	This bash script used HybPiper to generate consensus assemblies of Angiosperms353 loci data in batch 3.
Files:
namelist_batch3.txt
-	This is the list of sample names within batch 3.
TruSeq3-PE.fa
-	This is the TruSeq adapter fasta file necessary to run 00.trimmomatic_batch3.sh.
