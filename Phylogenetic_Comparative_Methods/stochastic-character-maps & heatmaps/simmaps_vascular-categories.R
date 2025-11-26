###############################################################################
# Stochastic Character Mapping for Variant Categories
# Joyce G. Onyenedum with help from ChatGPT
###############################################################################

library(phytools)
library(corHMM)
library(grDevices)
library(geiger)

# ---- Working dir & tree/data ----
setwd("/Users/joyce.onyenedum/Desktop")

tree <- read.nexus("CI-combined_chronograms.tre")
tree <- root(tree, outgroup = "Cupania_guatemalensis", resolve.root = TRUE)
tree <- ladderize(tree)

data <- read.delim("DatasetS1.txt", row.names = 1)

# ---- Trait selection & prep ----
trait_column <- "Variant_Category"
datum <- data[, trait_column, drop = FALSE]
datum <- datum[complete.cases(datum), , drop = FALSE]

species <- rownames(datum)
datas <- datum[[trait_column]]
names(datas) <- species
pruned_tree <- keep.tip(tree, species)

# ---- Explicit state ordering (biological order) ----
state_levels <- c(
  "Typical", "Procambial", "Cambial", "Ectopic-cambia",
  "Procambial+Cambial", "Procambial+Cambial+Ectopic-cambia",
  "Procambial+Ectopic-cambia", "Cambial+Ectopic-cambia"
)

# Force factor levels in the intended order (important!)
factor_states <- factor(datas, levels = state_levels)
if (any(is.na(factor_states))) {
  stop("Some values in 'datas' are not present in your 'state_levels' vector. Fix levels.")
}
numeric_states <- as.numeric(factor_states)

cor_data <- data.frame(Taxa = species, State = numeric_states)

# Print mapping for clarity
cat("\nNumeric coding of states (1 = first level):\n")
print(data.frame(State_Number = seq_along(state_levels), State_Name = state_levels))

#To get the simmaps, you have two options:
###OPTION 1: if you want to just read in the saved RDS file with the 1000 simmaps for categories, run line 52, THEN SKIP OVER to line 79
simmap.trees <- readRDS("simmap.trees_vc.RDS")

###OPTION 2: If you want to actually make the simmmaps of vascular categories, then run lines  57-76

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
simmap.trees <- makeSimmap(tree = pruned_tree, data = cor_data, model = Q, rate.cat = 1, nSim = 1000)
cor_data$State

# Save results
saveRDS(simmap.trees, file = "simmap.trees_vc.RDS")

######## PLOT SIMMAPS
state_levels <- c(
  "Typical", "Procambial", "Cambial", "Ectopic-cambia",
  "Procambial+Cambial", "Procambial+Cambial+Ectopic-cambia",
  "Procambial+Ectopic-cambia", "Cambial+Ectopic-cambia"
)
sim_states <- as.character(1:length(state_levels))
col_vc <- c("grey", "tomato", "#0099FF", "gold2",
            "deeppink3", "mediumpurple", "darkgreen", "darkblue")
names(col_vc) <- as.character(1:length(state_levels))
#check it worked
print(data.frame(State_Number = names(col_vc),
                 State_Name = state_levels,
                 Color = col_vc))

foo <- function(x) {
  # first state on each edge's map
  y <- sapply(x$maps, function(z) names(z)[1])
  names(y) <- x$edge[, 1]
  # internal node numbers are (ntips + 1) : (ntips + Nnode)
  node_ids <- (length(x$tip.label) + 1):(length(x$tip.label) + x$Nnode)
  y[as.character(node_ids)]
}

#Create node labels for posterior density
Nsim <- length(simmap.trees)
XX <- sapply(simmap.trees, foo)

# Build the pie fractions
pies <- t(
  apply(XX, 1, function(x) {
    tab <- table(factor(x, levels = sim_states))
    as.numeric(tab) / Nsim
  })
)

# Quick sanity checks:
range(pies)        # should be from 0 up to ≤1
all(rowSums(pies) > 0.999 & rowSums(pies) < 1.001)  # each row sums to ~1

#generate summary of stochastic maps with pies of posterior at nodes..
source("plot_simmap.R")
plot_simmap(time_tree = simmap.trees[[1]], 
            tree = simmap.trees[[1]], 
            simmaps = simmap.trees, 
            states = sim_states,
            show.tip.label = T, label.offset = .45,
            lwd = 2.5,
            label.cex = .45,
            colors = col_vc, edge.width=0, nt=10001)

nodelabels(pie=pies,cex=0.16,piecol=col_vc, lwd=1)
add.simmap.legend(leg = state_levels, colors = col_vc, prompt = TRUE, vertical = TRUE)
title(main = "Vascular Variant Categories", font.main = 2, line = -1)
