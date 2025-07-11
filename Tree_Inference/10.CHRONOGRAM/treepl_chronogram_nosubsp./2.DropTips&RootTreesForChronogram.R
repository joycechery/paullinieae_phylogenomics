# Load necessary packages
library(ape)
library(phytools)

# Set the working directory
setwd("/Users/joyce.onyenedum/Library/CloudStorage/Box-Box/Onyenedum_Lab/Joyce_Onyenedum/Projects/Paullinieae_Phylogeny/*MANUSCRIPT/11.SCRIPTS/Tree_Inference/10.CHRONOGRAM/treepl_chronogram_nosubsp.")

# Load the list of taxa to retain
data <- read.delim("chronogram-taxa.txt")
data_names <- as.character(data$chronogram.taxa)

# Import the tree file containing multiple trees
trees <- read.tree("pau_333s_351g_partitions.ufboot_100_names.tre")

# Define the outgroup for rooting the trees
outgroup <- c("Cupania_guatemalensis","Cupania_racemosa")

# Create a directory to save pruned trees if it doesn't exist
output_dir <- "prunedTrees"
if (!dir.exists(output_dir)) {
  dir.create(output_dir)
}
setwd(output_dir)

# Process each tree in the list
for (i in seq_along(trees)) {
  # Root the tree using the specified outgroup
  trees[[i]] <- root(trees[[i]], outgroup, resolve.root = TRUE)
  
  # Ladderize the tree to ensure a consistent branch order
  trees[[i]] <- ladderize(trees[[i]])
  
  # Prune the tree to keep only the specified taxa
  trees[[i]] <- keep.tip(trees[[i]], data_names)
  
  # Optionally, plot the pruned tree
  plotTree(trees[[i]], main = paste("Tree", i))
  
  # Save the pruned tree to a file
  output_filename <- paste0("tree_", i, ".tre")
  write.tree(trees[[i]], file = output_filename)
}

# Optionally, check the tip labels of the last tree to verify pruning
print(trees[[length(trees)]]$tip.label)

