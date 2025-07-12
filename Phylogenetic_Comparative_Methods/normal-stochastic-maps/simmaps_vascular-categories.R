#Script to generate stochastic maps for vascular categories by Dr. Joyce G. Onyenedum

library(phytools)
library(corHMM)
library(grDevices)
library(geiger)

setwd("/Users/joyce.onyenedum/Library/CloudStorage/Box-Box/Onyenedum_Lab/Joyce_Onyenedum/Projects/Paullinieae_Phylogeny/*MANUSCRIPT/11.SCRIPTS/Phylogenetic_Comparative_Methods/jgo/normal-stochastic-maps")
# Load and prepare tree
tree <- read.nexus("CI-combined_chronograms.tre")
root <- c("Cupania_guatemalensis")
tree <- root(tree, root, resolve.root = TRUE)
tree <- ladderize(tree)
plotTree(tree)

# Load data
data <- read.delim("DatasetS1.txt", row.names = 1)

# Select trait
trait_column <- "Variant_Category"
datum <- data[, c(trait_column), drop = FALSE]
datum <- datum[complete.cases(datum), , drop = FALSE]

# Prepare data for corHMM
species <- rownames(datum)
datas <- datum[[trait_column]]
names(datas) <- species
pruned_tree <- keep.tip(tree, species)

# Format data for corHMM: a data frame with Taxa and State columns
cor_data <- data.frame(Taxa = species, State = datas)
table(datas)

# Factorize the trait data
state_labels <-c("Typical" ,"Procambial" ,"Cambial" ,"Ectopic-cambia" ,"Procambial+Cambial" ,"Procambial+Cambial+Ectopic-cambia" ,"Procambial+Ectopic-cambia" ,"Cambial+Ectopic-cambia")
factor_states <- factor(datas, levels = state_levels)
numeric_states <- as.numeric(factor_states)

# Define color mapping
col_vc <- c("grey", "tomato", "#0099FF", "gold2", "deeppink3", "mediumpurple", "darkgreen", "darkblue")
names(col_vc) <- state_labels  # helpful for downstream plotting

# Create cor_data using numeric states
cor_data <- data.frame(Taxa = species, State = numeric_states)

#To get the simmaps, you have two options:
#OPTION 1: if you want to just read in the saved RDS file with the 1000 simmaps for categories, run line 49
#####simmap.trees_vc <- readRDS("simmap.trees_vc.RDS")

#OPTION 2: If you want to actually make the simmmaps of vascular categories, then run lines  54-73

# Fit models using corHMM
models <- c("ER", "SYM", "ARD")
fit_corHMM <- lapply(models, function(model) {
  corHMM(phy = pruned_tree, data = cor_data, model = model, rate.cat = 1, root.p = "yang")
})

# Choose best model using AIC weights
aic_weights <- aicw(sapply(fit_corHMM, function(x) x$AICc))[,3]
best_model_index <- which.max(aic_weights)
Q <- fit_corHMM[[best_model_index]]$solution

cat("\nBest-fit model:", models[best_model_index], "\nRate matrix (Q):\n")
print(Q)

# Simulate stochastic character maps using best-fit Q
set.seed(123)
simmap.trees_vc <- makeSimmap(tree = pruned_tree, data = cor_data, model = Q, rate.cat = 1, nSim = 1000)
cor_data$State

# Save results
saveRDS(simmap.trees_vc, file = "simmap.trees_vc.RDS")

# >>> Summarize simmap results <<<
summary_stm_vc <- describe.simmap(simmap.trees_vc)
print(summary_stm_vc)

#Create node labels for posterior density
foo<-function(x){
  y<-sapply(x$maps,function(x) names(x)[1])
  names(y)<-x$edge[,1]
  y<-y[as.character(length(x$tip)+1:x$Nnode)]
  return(y)
}


Nsim <- length(simmap.trees_vc)
# XX: a matrix, rows = internal nodes, cols = each simmap’s first state at the parent
XX <- sapply(simmap.trees_vc, foo)

# Build the pie fractions
pies <- t(
  apply(XX, 1, function(x) {
    tab <- table(factor(x, levels = states))
    as.numeric(tab) / Nsim
  })
)

# Quick sanity checks:
range(pies)        # should be from 0 up to ≤1
all(rowSums(pies) > 0.999 & rowSums(pies) < 1.001)  # each row sums to ~1

#generate summary of stochastic maps with pies of posterior at nodes..
source("plot_simmap.R")
plot_simmap(time_tree = simmap.trees_vc[[1]], 
            tree = simmap.trees_vc[[1]], 
            simmaps = simmap.trees_vc, 
            states = states,
            show.tip.label = T, label.offset = .45,
            lwd = 2.5,
            label.cex = .45,
            colors = col_vc, edge.width=0, nt=10001)

nodelabels(pie=pies,cex=0.16,piecol=col_vec, lwd=1)
add.simmap.legend(leg = state_labels, colors = col_vc, prompt = TRUE, vertical = TRUE)
title(main = "Vascular Variant Categories", font.main = 2, line = -1)


