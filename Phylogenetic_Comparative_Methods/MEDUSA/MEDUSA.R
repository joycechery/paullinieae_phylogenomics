# Written by Rebeca Hernández-Gutiérrez (2023) and edited by Israel L. Cunha Neto (2025)
# Diversification rate estimation using BAMM (http://bamm-project.org/index.html)

library(ape)
library(BAMMtools)

setwd("")

# Preparing input tree

# Reading tree
tree <- read.tree("CI-combined_chronograms.tre")

tree <-drop.tip(tree, tip = c("Talisia_nervosa", 
                             "Guindilia_trinervis")) 
tree

write.tree(phy = tree, file = "paullinieae_V1.tree")

# Now, we start BAMM!

library(BAMMtools)
  
# Priors before running BAMM
setBAMMpriors(read.tree("paullinieae_V1.tree"))
# automatically writes a .txt file called "my_priors.txt"; the priors should be copied and
# pasted in the control file for BAMM to run correctly.

#Run the control file using the command "bamm -c BAMM_control_Paullinieae.txt" on Terminal 

# Once you have a BAMM completed, follow the instructions below to analyze BAMM output and plotting.

# Check convergence of chains
mcmcout <- read.csv("Paullinieae_mcmc_out.txt", header=T)
plot(mcmcout$logLik ~ mcmcout$generation)

# Discard burn-in
burnstart <- floor(0.1 * nrow(mcmcout))
postburn <- mcmcout[burnstart:nrow(mcmcout), ]

# Check Effective Sample Size, should be more than 200
library(coda)
effectiveSize(postburn$N_shifts) #number of rate shifts sampled at each MCMC generation; tells how well the MCMC chain has sampled the posterior distribution of shift numbers
effectiveSize(postburn$logLik) #even though its Bayesian it runs likelihood estimates (in every generation)
#this output is the frequency (if > 200, convergence reached) 

# Now, we should obtain the number of rate shifts
# read tree (newick format)
tree <- read.tree("paullinieae_V1.tree")

# read event data file, one output of BAMM
edata <- getEventData(tree, eventdata = "Paullinieae_event_data.txt", burnin=0.1)
edata

edata$tipLambda #check the absolute values of the speciation (lambda) rate at the tips.

# Obtaining the probability of rate shifts in the data
shift_probs <- summary(edata)
shift_probs # List of possible number of shifts, each one is a model and each model has a posterior probability (prob).
#shifts     prob
#1       1 0.096500
#2       2 0.389000
#3       3 0.282000
#4       4 0.132000
#5       5 0.062500
#6       6 0.023200
#7       7 0.008890
#8       8 0.003440
#9       9 0.001000
#10     10 0.000889
#11     11 0.000333

# This indicates that the data most probably have 11 shifts in the diversification rate.

# Bayesian credible sets of shift configurations
# Identify the 95% credible set of distinct shift configurations (i.e., different combinations of where shifts occur).
css <- credibleShiftSet(edata, expectedNumberOfShifts=1, threshold=5, set.limit = 0.95) 

# Here is the number of distinct shift configurations in the data:
css$number.distinct
# [1] 7

# Here are the shift configurations
plot(css)

summary(css) #95 % credible set of rate shift configurations sampled with BAMM .
# There are 7 distinct shift configurations in the credible set . 
#The model with highest posterior probability (P = 0.52) support two diversification shifts.  

#rank     probability cumulative  Core_shifts
#1 0.51727586  0.5172759          2
#2 0.18975669  0.7070326          2
#3 0.14676147  0.8537940          1
#4 0.04288412  0.8966781          3
#5 0.02255305  0.9192312          2
#6 0.01910899  0.9383402          1
#7 0.01755361  0.9558938          3

# This last results indicate that the most common number of rate shifts was 2!

# Now, we can look at the configuration with the maximum a posteriori probability (MAP).
best <- getBestShiftConfiguration(edata, expectedNumberOfShifts=1) 
class(best)
#"bammdata"

#Plotting

#Plotting speciation rate
par(mar=c(1, 1, 1, 8))
plot(x=best, tau = 0.01, method="phylogram", labels = TRUE, legend=TRUE, spex ="s",
     lwd=0.5, cex=0.1, pal="RdYlBu")
addBAMMshifts(best, method="phylogram", par.reset=FALSE, cex=0.4, col= "red",
              bg="red") 
#shift circles are always the same, even though the plotted rate is speciation, 
#extinction or net diversification. The circles mean "a shift occurred,
#in the net diversification rate".

#Plotting extinction rate
plot(x=best, tau = 0.01, method="phylogram", labels = TRUE, legend=TRUE, spex ="e",
     lwd=0.5, cex=0.1, pal="RdYlBu")
addBAMMshifts(best, method="phylogram", par.reset=FALSE, cex=0.4, col= "red",
              bg="red") 

#Plotting net diversification rate
plot(x=best, tau = 0.01, method="phylogram", labels = TRUE, legend=TRUE, spex ="netdiv",
     lwd=0.5, cex=0.1, pal="RdYlBu")
addBAMMshifts(best, method="phylogram", par.reset=FALSE, cex=0.4, col= "red",
              bg="red") 

####### Rate Through Time #######

# Define start time as root age
st <- max(branching.times(tree))

# Define plot area
par(oma = c(4, 4, 2, 1), mar = c(5, 5, 4, 2)) 

# Speciation
plotRateThroughTime(edata,
                    ratetype = "speciation",
                    intervalCol = adjustcolor("red", 0.3),
                    avgCol = "red",
                    start.time = st,
                    ylim = c(0, 0.5),
                    cex.axis = 1.5,
                    cex.lab = 1.5,
                    yticks = 5
                    )
mtext("Speciation", side = 3, line = 0.5, font = 2, cex = 1.8)

# Extinction
plotRateThroughTime(edata,
                    ratetype = "extinction",
                    intervalCol = adjustcolor("blue", 0.3),
                    avgCol = "blue",
                    start.time = st,
                    ylim = c(0, 0.5),
                    cex.axis = 1.5,
                    cex.lab = 1.5,
                    yticks = 5
                    )
mtext("Extinction", side = 3, line = 0.5, font = 2, cex = 1.8)

# Net Diversification
plotRateThroughTime(edata,
                    ratetype = "netdiv",
                    intervalCol = adjustcolor("darkgreen", 0.3),
                    avgCol = "darkgreen",
                    start.time = st,
                    ylim = c(0, 0.5),
                    cex.axis = 1.5,
                    cex.lab = 1.5,
                    yticks = 5
                    )

mtext("Net Diversification", side = 3, line = 0.5, font = 2, cex = 1.8)

### End of code ###
