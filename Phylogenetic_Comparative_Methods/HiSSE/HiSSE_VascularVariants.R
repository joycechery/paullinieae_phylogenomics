# Written by Rebeca Hernández-Gutiérrez (2023) and edited by Israel L. Cunha Neto (2025)
# Hidden State Speciation and Extinction (HiSSE; Beaulieu & O'Meara 2016)

#install.packages("diversitree")
#install.packages("hisse")
#install.packages("deSolve")
#install.packages("GenSA")
#install.packages(c("subplex", "nloptr"))

#Loading libraries 
lapply(X=c("ape","geiger", "deSolve", "GenSA","subplex","nloptr"),FUN = library, character.only=TRUE)
suppressWarnings(library(diversitree))

#Trying HiSSE from https://github.com/thej022214/hisse 
#remotes::install_github("thej022214/hisse")
library(hisse)

# Trait-dependent diversification in Paullinieae

#Set directory
setwd("/Users/israelneto/Library/CloudStorage/Box-Box/Onyenedum_Lab/1 - Israel_L_Cunha_Neto/Projects/Paullinieae/NYU/Diversification_analyses/README/HiSSE/One_series")

# reading tree
tree <- read.nexus("CI-combined_chronograms.tre")

#Note: polymorphic state and missing data (?) are not allowed; they would have to be removed here if existing.

tree <-drop.tip(tree, tip = c("Cupania_guatemalensis","Talisia_nervosa",
                               "Guindilia_trinervis", "Paullinia_sp._a",
                               "Serjania_sp._a", "Serjania_sp._b"))

# reading data file
data <-read.csv("Dataset_HiSSE.csv") 

# reformat dataset
trait<-data.frame(data$Taxon, data$Vascular_Variant, row.names = data$Taxon)
trait

# Check if tree tip names and data names match
name.check(tree, data = trait)

# Begin HiSSE analyses

############---MODEL SELECTION---############
#Our anatomical database has 218 species; Paullinieae + outgroups have 707 spp.;
#Therefore, taxon sampling (218/707) is equal to 0.3074753

# The proportion of each state is BASED ON THE TOTAL ESTIMATED NUMBER OF SPECIES:

# State 0: 93/477= 0.194, rounding: 0.19 #there are 93 species with state 0 (typical growth) in the tree;
#there are 477 species with state 0 across Paullinieae + outgroups 

# State 1: 125/230= 0.5434, rounding: 0.54 #We found 230 species of Paullinieae with vascular variants

f<-c(0.19,0.54) # This is the sampling fraction that we will use for all the models;

#Begin model calculations 

# Model 1. Dull-null: turnover and extinction fraction equal for the two states
# no hidden character.
turnover <- c(1,1) #equal turnover rate
extinction.fraction<-c(1,1) #equal extinction fraction

# generating a transition matrix (q)
trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits=0)
print(trans.rates.bisse)
#    (0) (1)
#(0)  NA   2
#(1)   1  NA

# running Model 1
Null <- hisse(phy=tree, data=trait, f=f, turnover=turnover, 
                eps=extinction.fraction, hidden.states=FALSE,
                trans.rate=trans.rates.bisse)

# Model 2. BiSSE: turnover different; extinction fraction equal for the two states.
turnover <- c(1,2) #different turnover rate
extinction.fraction<-c(1,1) #equal extinction fraction
trans.rates.bisse <- TransMatMakerHiSSE(hidden.traits=0) # different transition rates

# running Model 2
BiSSElike <- hisse(phy=tree, data=trait, f=f, turnover=turnover,
                     eps=extinction.fraction, hidden.states=FALSE,
                     trans.rate=trans.rates.bisse)

# Model 3. HiSSE: turnover different for all the states (observed and hidden).
turnover <-c(1,2,3,4)
extinction.fraction<- c(1,1,1,1)
trans.rate.hisse<-TransMatMakerHiSSE(hidden.traits = 1)
print(trans.rate.hisse)

#(0A) (1A) (0B) (1B)
#(0A)   NA    2    5   NA
#(1A)    1   NA   NA    5
#(0B)    5   NA   NA    4
#(1B)   NA    5    3   NA

# running Model 3
HiSSE_trait <- hisse(phy=tree, data=trait, f=f, turnover=turnover,
                       eps=extinction.fraction, hidden.states=TRUE,
                       trans.rate=trans.rate.hisse)

# Model 4. CID2: turnover different for the hidden states and equal for the observed states,
# that is, character-independent diversification.
turnover<- c(1,1,2,2)
extinction.fraction<-c(1,1,1,1)
trans.rate<-TransMatMakerHiSSE(hidden.traits = 1, make.null = TRUE) #Three rates
trans.rate
#     (0A) (1A) (0B) (1B)
#(0A)   NA    2    3   NA
#(1A)    1   NA   NA    3
#(0B)    3   NA   NA    2
#(1B)   NA    3    1   NA

CID2 <- hisse(phy=tree, data=trait, f=f, turnover=turnover, eps=extinction.fraction, 
              hidden.states=TRUE,trans.rate=trans.rate)

# Model 5. CID4: turnover different for the hidden states but there are two additional
# states (C and D), equal turnover for the observed states, it has the same number of 
# parameters than the Model 3 (HiSSE) so this evaluates the character-independent 
# model as a null model to compare the HiSSE model.

turnover<-c(1,1,2,2,3,3,4,4)
extinction.fraction<-rep(1,8)
trans.rate<-TransMatMakerHiSSE(hidden.traits = 3, make.null = TRUE)
trans.rate
#     (0A) (1A) (0B) (1B) (0C) (1C) (0D) (1D)
#(0A)   NA    2    3   NA    3   NA    3   NA
#(1A)    1   NA   NA    3   NA    3   NA    3
#(0B)    3   NA   NA    2    3   NA    3   NA
#(1B)   NA    3    1   NA   NA    3   NA    3
#(0C)    3   NA    3   NA   NA    2    3   NA
#(1C)   NA    3   NA    3    1   NA   NA    3
#(0D)    3   NA    3   NA    3   NA   NA    2
#(1D)   NA    3   NA    3   NA    3    1   NA

CID4 <- hisse(phy=tree, data=trait, f=f, turnover=turnover, eps=extinction.fraction, 
              hidden.states=TRUE,trans.rate=trans.rate)

# End of model calculations

# Summarize the results of the models to compare the likelihoods and select models
library(tidyr)

models_summary <- as.data.frame(rbind(Null=c(Null$loglik, Null$AICc),
                                       BiSSE_like=c(BiSSElike$loglik, BiSSElike$AICc),
                                       HiSSE=c(HiSSE_trait$loglik, HiSSE_trait$AICc),
                                       CID2=c(CID2$loglik, CID2$AICc),
                                       CID4=c(CID4$loglik, CID4$AICc)))
models_summary                          
colnames(models_summary) <- c("loglik", "AICc")                             

# Order models in the dataframe
models_summary <- models_summary[order(models_summary$AICc), ]
models_summary
#loglik     AICc
#CID4       -754.5647 1525.818
#HiSSE      -757.5129 1536.089
#CID2       -767.6555 1547.709
#BiSSE_like -781.1782 1572.639
#Null       -803.0950 1614.378

# Obtaining Akaike values
aic <- models_summary$AICc

# Selecting the best model (i.e., minimum value of Akaike)
best <- min(aic)

# Obtaining the difference between Akaike values
delta <- aic - best

# Obtaining Akaike weights by first calculating the relative likelihood of the models using delta.
sumDelta <- sum(exp(-0.5 * delta))

# The Akaike weight for a model is this value divided by the sum of these values across all models.
# http://brianomeara.info/aic.html

w_models_variants <- round((exp(-0.5 * delta)/sumDelta), 3) %>%
  as.matrix()
colnames(w_models_variants) <- c("AICw")
w_models_variants

# Saving into a matrix
models_vascular_variants <- cbind(models_summary, w_models_variants) %>%
  as.data.frame() %>%
  write.csv("Hisse_models.csv")

# Save results as rdata
save(Null,HiSSE_trait,CID2,CID4,BiSSElike,
     file="./Paullinieae_hisse.rdata")

############---MODEL AVERAGING---############
# Model averaging will check if there are significant differences among 
#the models in diversification rates.

# Reading data file
data<-read.csv("Dataset_HiSSE.csv") 
data<-data.frame(data$Taxon, data$Vascular_Variant, row.names = data$Taxon)

f<-c(0.19,0.54) #The same fraction you used to calculate the models 

#Now, for each model we are interested in exploring further (e.g., model averaging),
#we estimate the likeliest states for internal nodes and tips of the phylogeny, 
#using the marginal reconstruction algorithm;
#For vascular variants, we will run model averaging including only 
#HiSSE and CID4 (AICw>0)

# Model 1: Null
#recon_null <- MarginReconHiSSE(phy = tree, data=data, f=f,
#                             pars = Null$solution, AIC = Null$AIC)
# Model 2: BiSSElike
#recon_bisse <- MarginReconHiSSE(phy = tree, data=data, f=f,
#                              pars = BiSSElike$solution, AIC = BiSSElike$AIC)

# Model 3: HiSSE_trait
recon_HiSSE_trait <- MarginReconHiSSE(phy = tree, data=data, f=f,
                             pars = HiSSE_trait$solution, AIC = HiSSE_trait$AIC,
                             hidden.states = 1)
# Model 4: CID2
#recon_CID2 <- MarginReconHiSSE(phy = tree, data=data, f=f,
#                             pars = CID2$solution, AIC = CID2$AIC,
#                             hidden.states = 2)

# Model 5: CID4
recon_CID4 <- MarginReconHiSSE(phy = tree, data=data, f=f,
                             pars = CID4$solution, AIC = CID4$AIC,
                             hidden.states = 4)

#save(recon_null,recon_bisse,recon_HiSSE_trait,recon_CID2,recon_CID4,
#     file = "Recon_objects.rdata")

save(recon_HiSSE_trait,recon_CID4, file = "Recon_objects.rdata")

# Model averaging including those with AICw > 0
models <- list(recon_CID4,recon_HiSSE_trait)

model_av <- GetModelAveRates(x= models,
                             type = "tips", AIC.weights = NULL,
                             bound.par.matrix = cbind(c(0,0,0,0,0), c(10,1,10,1,1)))

boxplot(net.div~state, data = model_av, notch=TRUE, col=c("burlywood4","cornflowerblue"),
        xlab="Vascular Variants", ylab = "Net diversification rate")

mean_difference<-t.test(x= model_av$state, y= model_av$net.div, var.equal = FALSE )
mean_difference

#Welch Two Sample t-test

#data:  model_av.1$state and model_av.1$net.div
#t = 13.195, df = 217.91, p-value < 2.2e-16
#alternative hypothesis: true difference in means is not equal to 0
#95 percent confidence interval:
#  0.3772533 0.5097364
#sample estimates:
#  mean of x mean of y 
#0.5733945 0.1298996  

####### End of script ######
