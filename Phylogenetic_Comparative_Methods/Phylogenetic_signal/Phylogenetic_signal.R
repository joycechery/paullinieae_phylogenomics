# Modified by Israel L. Cunha Neto from: https://blog.phytools.org/2018/02/how-to-fit-tree-transformation-for.html
#and https://github.com/wabarr/DFA-phylosim/blob/master/TestDiscreteTraitPhylogeneticSignal.R

library(geiger)

setwd("")

# Load data
tree <- read.nexus("CI-combined_chronograms.tre")

prunned_tree<-drop.tip(tree, tip = c("Cupania_guatemalensis","Talisia_nervosa",
                             "Guindilia_trinervis", "Paullinia_sp._a",
                             "Serjania_sp._a", "Serjania_sp._b"))

data <- read.csv("Dataset_PhylogeneticSignal.csv", row.names = 1)
datum <- as.data.frame(cbind(rownames(data), data[,1]), row.names = FALSE)
datum <- datum[complete.cases(datum), ]

# reformat data
species <- datum$V1
datas <- datum$V2
names(datas) <- species

#check data
name.check(prunned_tree, data = datas)

#transition models
fitER <- fitDiscrete(prunned_tree,datas,model= "ER")
fitARD <- fitDiscrete(prunned_tree,datas,model= "ARD")

#Compare the AIC scores from the model outputs below, picking the model with the lowest AIC
fitER$opt$aic
fitARD$opt$aic

# Phylogenetic signal
fit.lambda <-fitDiscrete(prunned_tree,datas,model= "ARD", transform ="lambda") #Fit discrete model with lambda transformation

# Tree with lambda = 0
tree0 <- rescale(prunned_tree, model = "lambda", 0) 

# Fit ARD model to rescaled tree
fit.lambda.0 <- fitDiscrete(tree0, datas, model = "ARD") # Fit null model with lambda = 0 (no signal)

# Likelihood Ratio Test
lnL1 <- fit.lambda$opt$lnL
lnL0 <- fit.lambda.0$opt$lnL
LR <- -2 * (lnL0 - lnL1)
pval <- pchisq(LR, df = 1, lower.tail = FALSE)

cat("Pagel's lambda =", round(fit.lambda$opt$lambda, 3), "\n")
cat("LogLik (lambda) =", lnL1, "| LogLik (lambda=0) =", lnL0, "\n")
cat("LR =", round(LR, 3), "| p =", signif(pval, 3), "\n")

### End of code ###