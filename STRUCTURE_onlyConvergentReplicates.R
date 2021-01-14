############################################################
#   INRA - ANR project Pseudorasbora
#
#       Only convergent replicates of STRUCTURE (will persist)
#
############################################################
# by Thomas Brazier
# brazier.thomas@gmail.com
# MSc internship

# Supervisor: Scott McCairns
#             INRA UMR ESE

# Purpose of this script is to remove all non-convergent replicates of STRUCTURE in a given output of STRUCTURE software before subsequent analyses

# A data set is first copied with the '_OC' suffix (e.g. '2b_Native' becomes '2b_Native_OC')
# Then the followig scrip is applied
# Analyses of DeltaK, STRUCTURE selector and CLUMPP can be run on the remaining replicates that are all convergent



#==========================================================
# LOADING ENVIRONMENT
#==========================================================

# clear global environment: remove all variables
rm(list=ls(all=TRUE))

#----------------------
# Loading packages
# SYSTEM
library(rstudioapi)
library(devtools)
# GENETICS
library(ade4)
library(adegenet)
# GRAPHIC LIBRARIES
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
# STATS
library(MCMCglmm)
library(coda)
library(mcmc)

#----------------------
# Loading variables & objects

# Get the directory of the file & set working directory
wd=dirname(rstudioapi::getSourceEditorContext()$path)
setwd(wd)
# Define data directory
datadir=paste(wd,"/Data",sep="")
# Define DIYABC directory
DIYABCdir=paste(gsub("/R","",wd),"/DIYABC",sep="")
# Define STRUCTURE directory
STRdir=paste(gsub("/R","",wd),"/STRUCTURE",sep="")
# Define the directory where to save figures
figuresdir=paste(wd,"/Figures",sep="")
# Define the directory where to save tables in .csv format
tablesdir=paste(wd,"/Tables",sep="")

#==========================================================
# LOADING SAMPLING CHAINS
#==========================================================
# We needed to import the trace of Alpha and posterior estimates of Ln P(D) (parameters to estimate and the sampling chain)
# One file per value of K
Kmin=1
path="2b_Native_OC" # The name of the data set to explore
Kmax=19
# path="2b1_Native_OC" # The name of the data set to explore
# Kmax=16
# path="3a1_Invasive_OC" # The name of the data set to explore
# Kmax=10

# Construct directory path
path=paste(STRdir,"/",path,"/Results/",sep="")
for (k in Kmin:Kmax) {
  print(k)
  assign(paste("sampling.chain.K",k,sep=""), read.table(paste(path,"sampling_chain_K",k,sep=""),header=T,sep=" "))
  # assign(paste("sampling.chain.K",k,sep=""), read.table(paste(STRdir,"/21_Native_1MillionIterations/Results/sampling_chain_K",k,sep=""),header=T,sep=" "))
}
nrep=20 # Number of replicates of the sampling chain

#==========================================================
# MAKE A LIST OF CONVERGENT CHAINS FOR SELECTION
#==========================================================
# We need a list of convergent replicates, so we can remove all non-convergent replicates for subsequent analyses (Harvester, CLUMPP)
# We need a list of regex to remove: 'outputKirunj' and 'Kirunj_f', the pattern 'Kirunj' is common to both

# We first used the Geweke diagnostic to select convergent replicates for Alpha and LnLikelihood

#========================================================================#
# Geweke Diagnostic of MC Chain Stability
#------------------------------------------------------------------------#
# if the whole chain is stationary, the means of the values early and late in the sequence should be similar
# convergence diagnostic 'Z' is the difference between the 2 means divided by the asymptotic standard error of their difference
# values of 'Z' near the extreme tails of the N(0,1) indicates lack of convergence
# can also estimate p-value of 'Z' from the normal distribution
# yields the probability that the divided chain means are different: p-value < 0.05 indicates a lack of convergence
# ?geweke.diag

removal_Alpha=c() # Intermediate list

# As a first approach, we used a stringent removal (removal of all sampling chains non-convergent for at least one parameter)
# for each value of K to test, removal_Alpha contains all replicates to remove
for (k in Kmin:Kmax) {
  cat("==============================================\nK number",k,"\n")
  sampling.chain=get(paste("sampling.chain.K",k,sep=""))
  for (rep in 1:nrep) {
    cat("----------------------------------------------\nReplicate number",rep,"\n")
    # A p-value can be computed for this Z-score
    # print(2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]), frac1=0.1, frac2=0.5)$z)))
    # count the number of convergent chains in the replicate (i.e. p-value > 0.05)
    if (sum(2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]), frac1=0.1, frac2=0.5)$z))>0.05)<(k/2)) {  # If there is the same number of convergent alphas as number of K
      removal_Alpha=c(removal_Alpha, paste("K",k,"run",rep,sep=""))
    }
  }
}
print(removal_Alpha)

##################
# Second step: automated detection of non-convergent sampling chains
removal_LnLik=c() # Intermediate list
# for each value of K to test
for (k in Kmin:Kmax) {
  cat("==============================================\nK number",k,"\n")
  sampling.chain=get(paste("sampling.chain.K",k,sep=""))
  for (rep in 1:nrep) {
    cat("----------------------------------------------\nReplicate number",rep,"\n")
    # A p-value can be computed for this Z-score
    # print(2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]), frac1=0.1, frac2=0.5)$z)))
    # count the number of convergent chains in the replicate (i.e. p-value > 0.05)
    if (2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),ncol(sampling.chain)-1]), frac1=0.1, frac2=0.5)$z))<0.05) {
      removal_LnLik=c(removal_LnLik, paste("K",k,"run",rep,sep=""))
    }
  }
}
print(removal_LnLik)

#############
# List of replicates to keep in subsequent analyses
# Only those found in removal_Alpha & removal_LnLik
removal=c() # An empty list where to store regex to remove
# removal=unique(c(removal_Alpha, removal_LnLik))
# removal=removal_Alpha # Only aplha convergence taken in account
removal=removal_LnLik # Only Ln Likelihood convergence taken in account
print(removal)
length(removal) # 200 replicates to remove

# Verify that there is at least 3 convergent replicates per value of K (3 needed for subsequent Evanno's method)
# Otherwise, print a warning
n.conv=c() # Store the number of convergent replicates
for (i in Kmin:Kmax) {
  # How many Ki in the list 'removal'
  if (nrep-length(grep(paste("K",i,"run", sep=""),removal)) > 2) {
    n.conv=c(n.conv,nrep-length(grep(paste("K",i,"run", sep=""),removal)))
    cat(nrep-length(grep(paste("K",i,"run", sep=""),removal)), "convergent replicates for K=", i,"\n")
  } else {
    warning(paste("Not enough convergent replicates for K =", i,"\n", sep=" "))
  }
}
# We now have a list of replicates to remove before subsequent analyses

# Remove files with R in the 'path' directory
# Remove both output of STRUCTURE and output of the console
for (r in removal) {
  file.remove(paste(path,r,"_f",sep=""))
  file.remove(paste(path,"output",r,".txt",sep=""))
}

#############
# Problem with this method: not the same number of replicates for each K, whereas the variance is important in the DeltaK method and plateau method
# Statistically, should'nt we reduce the number of replicates to the same number for each K, the minimum number of convergent replicates for a given K

# in our case, it means to keep only 6 replicates per K (min number of convergent replciates for a given K is 6)
n.remain.rep=min(n.conv)
# For each K
# Get a list of remaining replicates in the 'path' directory
remain=list.files(path)
for (i in Kmin:Kmax) {
  # How many Ki remaining in the directory
  cat((length(grep(paste("^K",i,"run", sep=""),remain))), "convergent replicates for K=", i,"\n")
  
  # Remove replicates in excess, randomly chosen
  file.remove(paste(path,remain[sample(grep(paste("^K",i,"run", sep=""),remain), (length(grep(paste("^K",i,"run", sep=""),remain)) - n.remain.rep), replace=FALSE)],sep=""))  
}

# Rename the conserved replicates in sequential order, for CLUMPP usage...
remain=list.files(path)
for (i in Kmin:Kmax) {
  # How many Ki remaining in the directory
  cat((length(grep(paste("^K",i,"run", sep=""),remain))), "convergent replicates for K=", i,"\n")
  
  # Remove replicates in excess, randomly chosen
  file.rename(paste(path,remain[grep(paste("^K",i,"run", sep=""),remain)],sep=""),
              paste(path,"K",i,"run",c(1:n.remain.rep),"_f",sep=""))  
}


#==========================================================
# DISCUSSION
#==========================================================
# Not enough replicates retained with stringent filtering (only 26 replicates for k=1:9)
# We adopted a less stringent filtering of Alphas --> keeping partially convergent replicates (at least k/2 convergent sampling chains)
# 46 replicates retained for K=1:9, with at least 3 replicates for each number of K




#==========================================================
# END
#==========================================================
