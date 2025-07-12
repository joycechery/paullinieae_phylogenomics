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

### TREE_INFERENCE SUBFOLDERS ####
-	00–01.Batch1_Trimmomatic&HybPiper_Assembly
-	00–01.Batch2_Trimmomatic&HybPiper_Assembly
-	00–01.Batch3_Trimmomatic&HybPiper_Assembly

Given the computational constraints on the NYU High Performance Computer (HPC), the cleaning of raw data and the subsequent HybPiper assemblies of those cleaned data were processed separated into three batches. These three batches are contained in the following subfolders within the main folder “Tree Inference”.

-	02–05.HybPiper_Downstream: This folder contains scripts and results to evaluate the HybPiper assembly, including creating a heatmap, sequence statistics, removing paralogs for downstream analysis, and the retrieval of sequences to generate alignments for tree inference.

-	06–07.Sequence_Alignments: This folder contains scripts to align and trim exons. Also within this subfolder is the concatenated alignment that is used for the iqtree analysis (pau_333s_351g.fasta).

-	08.IQ-TREE: This folder contains the script to run a partitioned maximum likelihood tree inference using IQ-TREE v. 1.6.12 on the mafft & post-trimal exon concatenated alignment. The iqtree results are presented in a subfolder named “iqtree_results”. A script to change the tree tip names from sample ID (e.g., “PAU24”) to species names (“Serjania_atrolineata_b”) is available (ChangeSampleNames_iqtree.sh)

-	09.ASTRAL: This folder contains the script to run individual gene trees for each of the 351 mafft + trimal exon alignments. The subfolders content: “trimmed_exons” contains 351 post mafft +trimal alignments; “exon_trees” contain 351 gene trees; “ASTRAL_analysis” contains the input trees and result of an ASTRAL analysis. 

-	10.CHRONOGRAM: Given the size of the dataset (351 genes and 333 samples) a BEAST analysis was not feasible, so we opted to use treePL. Here we used the iqtree consensus tree as a topological constraint to make 1000 bootstrap trees with varying branch lengths. Then we subsampled 100 random bootstrap trees to individually time calibrate using treePL. Those 100 chronograms were collapsed using TreeAnnotator to generate a chronogram where the height 95%HPD is equivalent to the uncertainty given the variability of branch lengths among the 100 bootstrap chronogram trees. 

###########
00–01.Batch1_Trimmomatic&HybPiper_Assembly
###########
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

###########
00–01.Batch2_Trimmomatic&HybPiper_Assembly
###########

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

###########
00–01.Batch3_Trimmomatic&HybPiper_Assembly
###########

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


###########
02–05.HybPiper_Downstream
###########

Scripts:

02.hybpiper_stats.sh
-	This command will summarize target enrichment and gene recovery efficiency for a set of samples. creating seq_length.csv and hybpiper_stats.tsv.

04.hybpiper_heatmap.sh
-	This bash script inputs seq_length.csv to generate recovery_heatmap.png of gene coverage per sample.

04.hybpiper_paralog.sh
-	This bash script searches for sequences flagged as paralogs and deletes them from the folders containing putative orthologous loci, so they are not included in downstream analyses

05.hybpiper_retrieve.sh
-	This bash script extracts the exons and supercontigs (introns and exons) and places them into separate directories.

Files:

seq_length.tsv

-	This file is the product of 02.hybpiper_stats.sh, detailing the length of each gene per sample.
hybpiper_stats.tsv

-	This file is the product of 02.hybpiper_stats.sh, detailing read mapping statistics per sample. 
recovery_heatmap.png

-	This image file is the product of 04.hybpiper_heatmap.sh, showing the percentage coverage for each sample for each gene.
namelist_all.txt

-	This is the sample namelist file needed to run 02.hybpiper_stats.sh, 04.hybpiper_paralog.sh, and 05.hybpiper_retrieve.sh.
Paullinieae_140s_352g_target_fixed.fasta

-	This is a custom fasta file composed of published assembled Angiosperms353 genes from Sapindaceae/Paullinieae used by HybPiper for sequence assembly.

###########
06–07.Sequence_Alignments
###########
Scripts:

06.mafft_exons.sh
-	This bash script inputs the unaligned exons resulting from 05.hybpiper_retrieve.sh, then produces maaft alignments for each exon. 

07.trimal_exons.sh
-	This bash script inputs the aligned exons from 06.mafft_exons.sh, then exports trimal-trimmed alignments.

Files:

namelist_regions_order.txt
-	This text file lists all gene names. This file is needed to run 06.mafft_exons.sh and 07.trimal_exons.sh

pau_333s_351g.fasta	
-	This is the concatenated alignment of 333 samples and 351 mafft aligned and post trimal exons. To create this file, the alignments produced by 07.trimal_exons.sh were imported into Geneious and concatenated using the “Concatenate Sequence or Alignment” tool, where exons were order from lowest (clean_4691) to highest number (clean_7628).

###########
08.IQ-TREE
###########

Scripts:

08.iqtree_exons.sh
-	This bash script will generate the maximum likelihood tree of the concatenated alignment, using the partition gene file (pau_333s_351g_partitions) 

iqtree.tre
-	This is the maximum likelihood tree with bootstrap support at the nodes with the tips assigned their species names (e.g., Paullinia_atrolineata_b)

pau_333s_351g_partitions.contree
-	This is the maximum likelihood tree with bootstrap support at the nodes with the sample ID numbers (PAU24)

pau_333s_351g.fasta
-	This is the concatenated alignment of 351 genes across 333 samples. The sample ID names are still present.

pau_333s_351g_partitions
-	This is the partition gene file, showing the coordinates of each of 351 genes across the concatenated alignment.

ChangeSampleNames_iqtree.sh
-	This script uses the function “sed” to find and replace sample ID (e.g., PAU24) in pau_333s_351g_partitions.contree species names (e.g., “Serjania_atrolineata_b”) in a file called iqtree.tre

###########
09.ASTRAL
###########

Scripts:

09.astral_exons.sh
-	This bash script will generate one gene tree for each of the alignments that was created by 07.trimal_exons.sh. 

Files:

namelist_regions_order.txt
-	This text file lists all gene names. This file is needed to run 09.astral_exons.sh

Subfolders: ASTRAL_analysis: This subfolder has files and scripts needed to generate a species tree using ASTRAL 5.7.8.

Tree Preparation scripts and files:

exon_trees.tre
-	 this is all “.contrees” (maximum likelihood consensus trees generated for each exon gene alignment)

exon_trees_names.tre
-	this is the exon_trees.tre, but with the names corrected from sample ID (e.g.,“PAU24”) to species name (e.g., “Serjania_atrolineata_b”) by the ChangeSampleNames_astral.txt script

ChangeSampleNames_astral.sh
-	This uses the function “sed” to find and replace sample ID (e.g.,“PAU24”) in exon_trees.tre to species name (e.g., “Serjania_atrolineata_b”) in a file called exon_trees_names.tre

ASTRAL Analysis files
-	astral.5.7.8.jar –– this is the astral program
-	lib –– these are the libraries associated with the astral.5.7.8.jar 
-	README.md–– this is the readme for astral.5.7.8.jar
-	astral.tre–– this is the result of the astral analysis
-	astral.pdf–– this is the result of the astral analysis as a pdf

-	 The workflow to generate the ASTRAL tree is outlined below:

1.	Tree Preparation for ASTRAL:
In terminal….
	-	navigate to the directory of exon trees
		-	cd exon_trees
	-	combine all consensus trees 
		-	cat *.contree > exon_trees.tre
	-	move exon_trees.tre to ASTRAL_analysis subfolder
		-	mv exon_trees.tre ../ ASTRAL_analysis
	-	Change sample ID names to species names to make exon_trees_names.tre
		-	./ChangeSampleNames_astral.sh”  

2.	Run ASTRAL:
In terminal…
	-	In the same directory as the exon_trees_names.tre, place the astral.5.7.8.jar and the associated lib/ folder
	-	Then run astral:
		-	java -jar astral.5.7.8.jar -i exon_trees_names.tre -o astral.tre
			-	This is the astral tree

###########
10.CHRONOGRAM 
###########

Subfolder: iqtree_bootstrapTrees. This subfolder contains input and output files for the iqtree constrained topology bootstrap analysis to generate 1000 trees with the same topology, but varying branch lengths.

-	Important Input files/scripts:
	-	iqtree2 
		-	 this is the iqtree2 program which has the constrained topology function
	-	iqtree_bootstraps.sh
		-	this is the script to execute the iqtree bootstrap analysis to generate 1000 trees with varying branch lengths 
	-	pau_333s_351g.fasta
		-	concatenated alignment of 333 samples and 351 exons alignments
	-	pau_333s_351g_partitions
		-	this file contains the gene partition information
	-	pau_333s_351g_partitions.contree
		-	 this is the consensus ML tree generated by 08.iqtree_exons.sh, which will serve here as a topological constraint

-	Important Output files:
	-	pau_333s_351g_partitions.ufboot
		-	result of this constrained iqtree boostrap analysis, containing 1000 bootstrap trees with the same topology, but varying branch lengths.

Subfolder: treepl_chronogram_nosubsp. This subfolder contains the input and output files associated with a treePL time-calibration analysis on the 100 randomly sampled bootstrap trees from the 1000 trees created through the constrained topology iqtree analysis.

Scripts:
	-	1. extract_bootstrap.py
		-	this script will randomly sample 100 bootstrap trees from the 1000 bootstrap tree file “pau_333s_351g_partitions.ufboot” to create a file called “pau_333s_351g_partitions.ufboot_100.tre”

	-	1a.ChangeSampleNames_treePL.sh
		-	This file uses the function “sed” to find and replace sample ID (e.g., PAU24) in pau_333s_351g_partitions.ufboot_100.tre  into species names (“Paullinia_pinnata_a”). 

	-	2.DropTips&RootTreesForChronogram.R
		-	This script uses the chronogram-taxa.txt namelist to prune each of the 100 input tree (pau_333s_351g_partitions.ufboot_100.tre), and exports each pruned tree into folder called “prunedTrees”

	-	3.TreePL_loop.py
		-	This script uses treePL and the “config_template.txt” file to loop through each of the pruned trees to time calibrate each of them and export each chronogram into a subfolder called “chronograms”
	
	-	4.Chronogram.R
		-	This script was used to take the raw chronogram tree file and make it into a figure for publication.

Files:
	-	chronogram-taxa.txt
		-	This file list all the names to be included in the chronogram, which are one per taxa. All subspecies and varieties were collapsed under the species name. This file is needed by DropTips&RootTreesForChronogram.R

	-	config.txt 
		-	This file tabulates the parameters for the penalized likelihood analysis on each of the 100 trees, including three calibration points. This file is needed 3.TreePL_loop.py

	-	pau_333s_351g_partitions.ufboot_100.tre
		-	random subset of 100 of the 1000 boostrap trees extracted via 1.extract_bootstrap.py to be used for treePL analyses.

	-	pau_333s_351g_partitions.ufboot_100_names.tre
		-	random subset of 100 of the 1000 boostrap trees extracted via 1.extract_bootstrap.py to be used for treePL analyses with sample  ID (e.g.,“PAU24”) replaced by species name (e.g., “Serjania_atrolineata_b”)

The workflow to generate the chronogram is outlined below:

1.	Generate bootstrap iqtrees with constrated topology and varying branch lengths
	-	In terminal…
		-	Place the following files within a single directory:
			-	Iqtree2
			-	pau_333s_351g.fasta,
			-	pau_333s_351g_partitions 
			-	pau_333s_351g_partitions.contree
	-	Run iqtree_bootstraps.sh 
		-	This will create the output “pau_333s_351g_partitions.ufboot”, which is the bootstrap trees

2.	Sample 100 randomly selected trees to root and prune
	-	In terminal
		i.	python3 1.extract_bootstrap.py
			-	this will  randomly subset 100 trees from the 1000 trees within pau_333s_351g_partitions.ufboot 
		ii.	./1a.ChangeSampleNames_treePL.sh
			-	This will change sample ID to species names in each of the 100 trees
	-	In R
		i.	Run the script 2.DropTips&RootTreesForChronogram.R
			1.	This will prune the tips down to the taxa listed in chronogam-taxa.txt, and generate the “prunedTrees” folder

3.	Time calibrate trees
	-	In terminal (navigate to the directory with “prunedTrees” )
		i.	python3 3.TreePL_loop.py
			-	this will use treePL to time calibrate each of the pruned trees and deposit them into a new folder called “chronograms”

4.	Generate chronogram
	-	In terminal (within the newly generated “chronograms” folder)
		i.	cat *.tre > combined_chronograms.tre
		ii.	Open combined_chronograms.tre in Figtree and export as a nexus tree  combined_chronograms.nexus

-	Download and open TreeAnnotator (https://beast.community/treeannotator) 
		i.	Settings:
			1.	Input Tree “combined_chronograms.nexus”  
			2.	Target tree type “Maximum clade credibility tree”
			3.	Node heights “Median Heights”
		ii.	Output Tree CI-combined_chronograms.tre
			-	This is is the chronogram.

## PHYLOGENETIC COMPARATIVE METHODS SUBFOLDERS ###

Subfolder: normal-simmaps. This subfolder contains the code and data for the normal stochastic mapping analyses (Fig. S11 and S12).

-	Script: simmaps_vascular-categories.R (will need to source plot_simmap.R)
-	Data: DatasetS1.txt
-	Input Tree: CI-combined_chronograms
-	Simmap results: simmap.trees_vc.RDS

-	Script: simmaps_vascular-patterns.R  (will need to source plot_simmap.R)
-	Data: DatasetS1.txt
-	Input Tree: CI-combined_chronograms
-	Simmap results: simmap.trees_vp.RDS

Subfolder: paramo-stochastic-maps: This subfolder contains the code and data for ancestral state reconstruction using PARAMO (Fig. 4).

-	Script: PARAMO_stochastic-maps.R
-	Data: data_PARAMO.txt
-	Input Tree: CI-combined_chronograms
-	Simmap results: results in subfolder called simmaps

Subfolder: Phylogenetic signal. This subfolder contains the code and data to estimate phylogenetic signal with vascular variants. 

-	Script: Phylogenetic_signal.R
-	Data: Dataset_PhylogeneticSignal.csv
	-	Dataset of the presence and absence of vascular variants. 
-	Input Tree: CI-combined_chronograms
	-	Original chronogram

Subfolder: MEDUSA. This subfolder contains the code, data, and input files for estimating diversification rate shifts using MEDUSA.

-	Script: MEDUSA.R
-	Data: Richness.csv
	-	The total number of species per genus (Table S5) is equally distributed across the species from the respective genus in the tree. This is the primary input.
-	Input Tree: CI-combined_chronograms_duplicated
	-	The tree was duplicated using a text editor to create a “multiPhylo” object (for plotting reasons)
-	Additional files/Outputs: Medusa_results.rdata: contains the objects of the analysis. 

Subfolder: BAMM. This subfolder contains the code, data, and input files for estimating diversification rate shifts using BAMM.

-	Script: BAMM.R
-	Data: sample_probs.txt
	-	Sampling fraction is similar to Richness in MEDUSA but formatted for BAMM. This is the primary input. 
-	Input Tree: CI-combined_chronograms
	-	Original chronogram
-	Additional files:
	-	Intermediate outputs/inputs
		-	Paullinieae_V1.tree
			-	A new Newick tree is created at the beginning of the analysis. 
	-	myPriors.txt
		-	Output of the BAMM analysis; used to indicate the parameter of the phylogeny in the BAMM_control_Paullinieae.txt. 
	-	BAMM_control_Paullinieae.txt  
		-	File used to indicate parameters for BAMM analysis; it describes the name and outputs of the following intermediate files:
			-	chain_swap.txt
			-	Paullinieae_event_data.txt
			-	Paullinieae_mcmc_out.txt

Subfolder: HiSSE. This subfolder contains the code and data for estimating trait-dependent diversification using HiSSE.

-	Script: HiSSE_VascularVariants.R
-	Data: Dataset_HiSSE.csv
		-	Dataset of the presence and absence of vascular variants, tendrils, lianas, and zygomorphic flowers. 
-	Input Tree: CI-combined_chronograms
		-	Original chronogram
-	Additional files: 
	-	Intermediate outputs/inputs
		-	Paullinieae_hisse.rdata 
			-	Contains the objects from the model selection step. 
		-	 Recon_objects.rdata
			-	Contains the objects from the model averaging ste




