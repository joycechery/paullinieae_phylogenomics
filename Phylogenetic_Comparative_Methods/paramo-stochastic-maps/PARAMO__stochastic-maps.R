#This is modified from the tutorial 2 here:  https://github.com/diegosasso/ontophylo_tutorials/blob/main/ontophylo_tutorial2_paramo.Rmd
#final version by Dr. Joyce G. Onyenedum with the aid of ChatGPT

## Load packages.

library(ontophylo)
library(ape)
library(phytools)
library(geiger)
library(corHMM)
library(grDevices)
library(tidyverse)
library(ggplot2)
library(cowplot)
library(grid)


#set directory
setwd("/Users/joyce.onyenedum/Library/CloudStorage/Box-Box/Onyenedum_Lab/Joyce_Onyenedum/Projects/Paullinieae_Phylogeny/*MANUSCRIPT/6.PCM/paramo-based-ASR/pau_333s_351g/")

## Load data; #this is a dataset reorganized to account for dependencies
paul_mat<- read.delim("data_PARAMO.txt")
head(paul_mat)
paul_mat<- na.omit(paul_mat)
paul_mat$Taxa

# Import tree, root, and ladderize tree
paul_tree <- read.nexus("CI-combined_chronograms.tre")
root<- c("Cupania_guatemalensis")
paul_tree<- root(paul_tree, root, resolve.root=TRUE)
paul_tree<-ladderize(paul_tree)
plotTree(paul_tree)

# Now prune the tree based on the data available
data_names<-paul_mat$Taxa
data_names<-as.character(data_names)
paul_tree<-keep.tip(paul_tree, data_names)
plotTree(paul_tree)
paul_tree$tip.label

# Set some parameters.
n_stm = 1000
res = 500
set.seed(42)

# Create folder to store simmap objects.
dir.create("/Users/joyce.onyenedum/Library/CloudStorage/Box-Box/Onyenedum_Lab/Joyce_Onyenedum/Projects/Paullinieae_Phylogeny/*MANUSCRIPT/6.PCM/paramo-based-ASR/pau_333s_351g/simmaps")

# Set a vector with all character names.
PAUL_ANAT <- list(Typical = c("Typical"), Procambial_Variants= c("Procambial_Variants"), Cambial_Variants= c("Cambial_Variants"), Ectopic_cambia= c("Ectopic_cambia"))

#make my dataframe into a list object
br_chars <- unlist(PAUL_ANAT, use.names = FALSE)

for (i in 1:length(br_chars)) {
  
  cat(paste0("\n", "Working on: ", br_chars[i], ": ", Sys.time(), "\n"))
  
  # Get character vector.
  char <- cbind(paul_mat$Taxa, paul_mat[[br_chars[i]]])
  
  # Set candidate models.
  models <- c("ER", "SYM", "ARD")
  
  fit_corHMM <- vector(mode = "list", length = length(models))
  
  for (j in 1:length(models)) {
    # Fit model with corHMM.
    fit_corHMM[[j]] <- corHMM(phy = paul_tree, data = char, model = models[[j]], 
                              rate.cat = 1, root.p = "yang")
  }
  
  # Get best model.
  w <- aicw(sapply(fit_corHMM, function(x) x$AICc))[,3]
  
  # >>> Print selected model <<<
  selected_model <- models[which.max(w)]
  cat("Selected model for", br_chars[i], "is:", selected_model, "\n")
  
  # Set Q matrix.
  Q <- fit_corHMM[[min(which(w == max(w)))]]$solution
  
  # >>> Print rate matrix <<< 
  cat("\nRate matrix (Q) for", br_chars[i], ":\n")
  print(Q)
  
  # Simulate stochastic maps.
  stm <- makeSimmap(tree = paul_tree, data = char, model = Q, rate.cat = 1, nSim = n_stm)
  
  # >>> Summarize simmap results <<<
  summary_stm <- describe.simmap(stm)
  cat("\nSimmap summary for", br_chars[i], ":\n")
  print(summary_stm)
  
  # Discretize trees.
  stm_discr <- lapply(stm, function(x) discr_Simmap_all(x, res = res) )
  stm_discr <- do.call(c, stm_discr)
  
  # Save RDS files.
  saveRDS(stm_discr, file = paste0("simmaps/", br_chars[i], ".RDS"))
  
}

# STEP 2. PARAMO: amalgamating stochastic character maps.
# Set a vector with all character names.
br_chars <- unlist(PAUL_ANAT, use.names = FALSE)

# Create temporary list to store discretized maps from individual characters.
MAPS <- setNames(vector(mode = "list", length = length(br_chars)), br_chars)

# Import all discretized maps from characters of a given anatomical region.
for (k in 1:length(br_chars)) {
  MAPS[[k]] <- readRDS(paste0("simmaps/", br_chars[[k]], ".RDS"))
}

# Amalgamate all individual characters as a single complex character.
amalg_simmaps <- paramo.list(br_chars, tree.list = MAPS, ntrees = n_stm)
saveRDS(amalg_simmaps, file = "stm_amalg_anato.rds")
amalg_simmaps<- readRDS("stm_amalg_anato.rds")
get_rough_state_cols(amalg_simmaps[[1]])

 #####RECTANGULAR TREE OF ALL SIMMAPS

#Set character states and color vector
col_vec<- c("grey80", "tomato",  "violet", "gold2", "#3300FF", "lightblue", "#0099FF",  "#CCFA00", "mediumpurple2", "forestgreen", "darkorange" , "brown", "tan" , "deeppink3",  "lightpink" ,"black", "mediumspringgreen" ,"magenta3", "#A8B200","red", "#AA606F" , "#7C7298" , "maroon",  "#F9F333", "#8B7388" , "#95569D" , "#7A8C71" , "#CB982F" , "darkgreen" ,"#CC8FAD" , "cyan")
states <- c("2111" ,"1311" ,"1211" ,"1112" ,"1121" ,"1131" ,"1141" ,"1312" ,"1322" ,"1321" ,"1212" ,"1122" ,"1151" ,"1332" ,"1111" ,"2112" ,"2121" ,"2122" ,"2131" ,"2132" ,"2141" ,"2142" ,"2151" ,"2152" ,"2211" ,"2311" ,"2312" ,"1231" ,"1331" ,"1341" ,"1142")

foo<-function(x){
  y<-sapply(x$maps,function(x) names(x)[1])
  names(y)<-x$edge[,1]
  y<-y[as.character(length(x$tip)+1:x$Nnode)]
  return(y)
}

XX<-sapply(amalg_simmaps,foo)
pies<-t(apply(XX,1,function(x,levels,Nsim) summary(factor(x,levels))/Nsim,levels=states,Nsim=1000))

Nsim <- length(amalg_simmaps)
# XX: a matrix, rows = internal nodes, cols = each simmap’s first state at the parent
XX <- sapply(amalg_simmaps, foo)

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
plot_simmap(time_tree = amalg_simmaps[[1]], 
            tree = amalg_simmaps[[1]], 
            simmaps = amalg_simmaps, 
            states = states,
            show.tip.label = T, label.offset = .45,
            lwd = 2.5,
            label.cex = .45,
            colors = col_vec, edge.width=0, nt=10001)

#=nodelabels(pie=pies,cex=0.16,piecol=col_vec, lwd=1)
add.simmap.legend(leg=states,
                  colors=col_vec,
                  prompt=TRUE,
                  vertical=TRUE)

title(main = "PARAMO (plot_simmap.R) vascular variants in Paullinieae", font.main = 2, line = -1)


##### PLOT THIS AS A CIRCULAR TREE USING REVGADGETS
devtools::install_github("revbayes/revgadgets@stochastic_map")
library(RevGadgets)
library(phytools)

# get a phylo object for one tree 
dat<- readRDS("stm_amalg_anato.rds")
tree <- ape::as.phylo(dat[[1]])
ape::write.tree(tree, "temp.tre")
tree <- RevGadgets::readTrees("temp.tre")

#how many states do I have?
get_rough_state_cols(amalg_simmaps[[1]])

#Set character states and color vector
col_vec<- c("grey80", "tomato",  "violet", "gold2", "#3300FF", "lightblue", "#0099FF",  "#CCFA00", "mediumpurple2", "forestgreen", "darkorange" , "brown", "tan" , "deeppink3",  "lightpink" ,"black", "mediumspringgreen" ,"magenta3", "#A8B200","red", "#AA606F" , "#7C7298" , "maroon",  "#F9F333", "#8B7388" , "#95569D" , "#7A8C71" , "#CB982F" , "darkgreen" ,"#CC8FAD" , "cyan")
patterns <-("Typical" ,"Compound" ,"Divided" ,"Ectopic_cambia" ,"Phloem_wedges" ,"Fissured" ,"Lobed" ,"Compound+Ectopic_cambia" ,"Compound+Phloem_wedges+Ectopic_cambia" ,"Compound+Phloem_wedges" ,"Divided+Ectopic_cambia" ,"Phloem_wedges+Ectopic_cambia" ,"Lobed+Phloem_wedges" ,"Compound+Fissured+Ectopic_cambia" ,"Nothing" ,"Typical+Ectopic_Cambia" ,"Typical+Phloem_wedges" ,"Typical+Phloem_wedges+Ectopic_Cambia" ,"Typical+Fissured" ,"Typical+Fissured+Ectopic_Cambia" ,"Typical+Lobed" ,"Typical+Lobed+Phloem_wedges+Ectopic_Cambia" ,"Typical+Lobed+Phloem_wedges" ,"Typical+Lobed+Phloem_wedges+Ectopic_Cambia" ,"Typical+Divided" ,"Typical+Compound" ,"Typical+Compound+Ectopic_Cambia" ,"Divided+Fissured" ,"Compound+Fissured" ,"Compound+Lobed" ,"Lobed+Ectopic_Cambia")
states <- c("2111" ,"1311" ,"1211" ,"1112" ,"1121" ,"1131" ,"1141" ,"1312" ,"1322" ,"1321" ,"1212" ,"1122" ,"1151" ,"1332" ,"1111" ,"2112" ,"2121" ,"2122" ,"2131" ,"2132" ,"2141" ,"2142" ,"2151" ,"2152" ,"2211" ,"2311" ,"2312" ,"1231" ,"1331" ,"1341" ,"1142")
names(col_vec) <- states

# process the maps
processed_maps <- processStochMaps(tree, simmap = dat, states = states)

plotStochMaps(tree = tree,
              maps = processed_maps,
            #  color_by = "MAP",
               colors = col_vec,
              tree_layout = "circular",
              tip_labels = TRUE,line_width = .6,
            tip_labels_italics = TRUE,
            tip_labels_size = 1.6, tip_labels_offset = .5)


# Reassociate colors with descriptive pattern labels instead of codes
# Assume states and patterns are in the same order
names(col_vec) <- patterns

# Create dummy data — one row per pattern
legend_df <- data.frame(
  pattern = factor(patterns, levels = patterns),
  x = 1,
  y = 1
)

# Create dummy plot just to extract the legend
dummy_plot <- ggplot(legend_df, aes(x = x, y = y, fill = pattern)) +
  geom_point(shape = 21, size = 4) +
  scale_fill_manual(values = col_vec) +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 10),
    legend.text = element_text(size = 8)
  ) +
  guides(fill = guide_legend(title = "Anatomical Pattern", ncol = 3))

# Extract only the legend
legend_only <- cowplot::get_legend(dummy_plot)

# Display just the legend
grid.newpage()
grid.draw(legend_only)
