#Script to generate stochastic maps for vascular patterns by Dr. Joyce G. Onyenedum

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
trait_column <- "Variant_Pattern"
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

# Factorize the trait data & Define color mapping
col_vp <- c("grey", "tomato", "violet", "gold2", "purple3", "lightblue","blue", "mediumpurple1", "forestgreen", "orange","brown4", "#CCFA00", "deeppink3", "tan")   
state_levels<-c("Typical", "Compound", "Divided", "Ectopic-cambia", "Phloem-wedges", "Fissured", "Lobed", "Compound+Ectopic-cambia", "Compound+Fissured+Ectopic-cambia", "Compound+Phloem-wedges", "Compound+Phloem-wedges+Ectopic-cambia", "Divided+Ectopic-cambia", "Lobed+Phloem-wedges", "Phloem-wedges+Ectopic-cambia")
factor_states <- factor(datas, levels = state_levels)
numeric_states <- as.numeric(factor_states)
state_labels <- levels(factor_states)
names(col_vp) <- state_labels  # helpful for downstream plotting

# Create cor_data using numeric states
cor_data <- data.frame(Taxa = species, State = numeric_states)

#To get the simmaps, you have two options:
#OPTION 1: if you want to just read in the saved RDS file with the 1000 simmaps for categories, run line 49
#######simmap.trees_vp <- readRDS("simmap.trees_vp.RDS")

#OPTION 2: If you want to actually make the simmmaps of vascular patterns, then run lines  54-73

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
simmap.trees_vp <- makeSimmap(tree = pruned_tree, data = cor_data, model = Q, rate.cat = 1, nSim = 1000)
cor_data$State

# Save results
saveRDS(simmap.trees_vp, file = "simmap.trees_vp.RDS")

# >>> Summarize simmap results <<<
summary_stm_vp <- describe.simmap(simmap.trees_vp)
print(summary_stm_vp)

#Create node labels for posterior density
foo<-function(x){
  y<-sapply(x$maps,function(x) names(x)[1])
  names(y)<-x$edge[,1]
  y<-y[as.character(length(x$tip)+1:x$Nnode)]
  return(y)
}

Nsim <- length(simmap.trees_vp)
# XX: a matrix, rows = internal nodes, cols = each simmap’s first state at the parent
XX <- sapply(simmap.trees_vp, foo)

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
plot_simmap(time_tree = simmap.trees_vp[[1]], 
            tree = simmap.trees_vp[[1]], 
            simmaps = simmap.trees_vp, 
            states = states,
            show.tip.label = T, label.offset = .45,
            lwd = 2.5,
            label.cex = .45,
            colors = col_vp, edge.width=0, nt=10001)

nodelabels(pie=pies,cex=0.16,piecol=col_vp, lwd=1)
add.simmap.legend(leg = state_labels, colors = col_vp, prompt = TRUE, vertical = TRUE)
title(main = "Vascular Variant Patterns", font.main = 2, line = -1)


