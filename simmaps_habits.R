###############################################################################
# Stochastic Character Mapping for Habits
# By Joyce G. Onyenedum and Israel L. Cunha Neto with help from ChatGPT
###############################################################################

library(phytools)
library(corHMM)
library(grDevices)
library(geiger)

setwd("/Users/joyce.onyenedum/Desktop")

# Load and prepare tree
tree <- read.nexus("CI-combined_chronograms.tre")
root <- c("Cupania_guatemalensis")
tree <- root(tree, root, resolve.root = TRUE)
tree <- ladderize(tree)
plotTree(tree)

# Load data
data <- read.delim("DatasetS1.txt", row.names = 1)

# Select trait
trait_column <- "Habit"
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
state_levels <-c("Tree" ,"Liana" ,"Shrub")
factor_states <- factor(datas, levels = state_levels)
numeric_states <- as.numeric(factor_states)

# Define color mapping
col_habits <- c("grey", "tomato", "#0099FF")
names(col_habits) <- state_levels  

# Create cor_data using numeric states
cor_data <- data.frame(Taxa = species, State = datas) 

#To get the simmaps, you have two options:
#OPTION 1: if you want to just read in the saved RDS file with the 1000 simmaps for categories, run line 52 THEN SKIP TO line 79
simmap.trees_habits <- readRDS("simmap.trees_habits.RDS")

#OPTION 2: If you want to actually make the simmmaps of habits from scratch, then run lines  57 onward

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
simmap.trees_habits <- makeSimmap(tree = pruned_tree, data = cor_data, model = Q, rate.cat = 1, nSim = 1000)
cor_data$State

# Save results
saveRDS(simmap.trees_habits, file = "simmap.trees_habits.RDS")

# >>> Summarize simmap results <<<
summary_stm_habits <- describe.simmap(simmap.trees_habits)
print(summary_stm_habits)

#Create node labels for posterior density
foo<-function(x){
  y<-sapply(x$maps,function(x) names(x)[1])
  names(y)<-x$edge[,1]
  y<-y[as.character(length(x$tip)+1:x$Nnode)]
  return(y)
}


Nsim <- length(simmap.trees_habits)
# XX: a matrix, rows = internal nodes, cols = each simmap’s first state at the parent
XX <- sapply(simmap.trees_habits, foo)

# Build the pie fractions
pies <- t(
  apply(XX, 1, function(x) {
    tab <- table(factor(x, levels = state_levels)) 
    as.numeric(tab) / Nsim
  })
)

# Quick sanity checks:
range(pies)        # should be from 0 up to ≤1
all(rowSums(pies) > 0.999 & rowSums(pies) < 1.001)  # each row sums to ~1

#generate summary of stochastic maps with pies of posterior at nodes..
source("plot_simmap.R")
plot_simmap(time_tree = simmap.trees_habits[[1]], 
            tree = simmap.trees_habits[[1]], 
            simmaps = simmap.trees_habits, 
            states = state_levels,
            show.tip.label = T, label.offset = .45,
            lwd = 2.5,
            label.cex = .45,
            colors = col_habits, edge.width=0, nt=10001)

nodelabels(pie=pies,cex=0.16,piecol=col_habits, lwd=1) 
add.simmap.legend(leg = state_levels, colors = col_habits, prompt = TRUE, vertical = TRUE)
title(main = "Habits", font.main = 2, line = -1)


