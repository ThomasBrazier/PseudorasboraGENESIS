############################################################################
#   INRA - ANR project Pseudorasbora
#
#       STRUCTURE preparation
#
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com
# MSc internship

# Supervisor: Scott McCairns
#             INRA UMR ESE


#==========================================================
# LOADING ENVIRONMENT
#==========================================================

# clear global environment: remove all variables
rm(list=ls(all=TRUE))

#----------------------
# Loading packages
library(ade4)
library(adegenet)
library(hierfstat)
library(rstudioapi)
library(MCMCglmm)
library(coda)
library(RColorBrewer)


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

# Load Pparva object
load(paste(datadir,file="/Pparva.clean.Rda",sep=""))


# Export genind to STRUCTURE data frame
# with function provided by Lindsay V. Clark (2015)
source("Sources/genind2structure.R") # A function to convert a genind object to an object formated for STRUCTURE

# All populations
labelPops=data.frame(unique(Pparva@pop),
                     1:length(unique(Pparva@pop)))
colnames(labelPops)=c("Location_name","Location_index")
write.table(labelPops,"Data/STRUCTURE/labelPops.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# Convert genind object Pparva to a STRUCTURE format
genind2structure(Pparva,"Data/STRUCTURE/Pparva.structure.txt",pops=TRUE)
# # Reformating data file
# Pparva.structure=read.table("Data/STRUCTURE/Pparva.structure.txt",sep="\t")
# ## Remove first line
# Pparva.structure=Pparva.structure[-1,]
# ## Add population info
# Pparva.structure=cbind(Pparva.structure[,1],rep(NA,nrow(Pparva.structure)),Pparva.structure[,2:ncol(Pparva.structure)])
# for (i in 1:nrow(Pparva.structure)) {
#   # get the name of the pop from indiv. name
#   Pparva.structure[i,2]=gsub("_[0-9]*","",Pparva.structure[i,1])
#   # replace the pop. name by its index in labelPops
#   Pparva.structure[i,2]=labelPops[which(labelPops[,1]==Pparva.structure[i,2]),2]
# }
# ## Change label of ind. to integer (required by STRUCTURE)
# ## Change name of individual by its index in the list of names
# index=as.character(unique(Pparva.structure[,1]))
# write.table(index,"Data/STRUCTURE/ind.Index.structure.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
# Pparva.structure[,1]=as.character(Pparva.structure[,1])
# for (i in 1:length(Pparva.structure[,1])) {
#   Pparva.structure[i,1]=match(Pparva.structure[i,1],index)
# }
# rm(index)
# ## Add same location to row 2 of individuals
# index1=seq(1,nrow(Pparva.structure)-1,2)
# index2=seq(2,nrow(Pparva.structure),2)
# Pparva.structure[index2,2]=Pparva.structure[index1,2]
# Pparva.structure[index2,2]==Pparva.structure[index1,2]
# rm(index1)
# rm(index2)
# 
# # Rewrite the complete data set for STRUCTURE
# write.table(Pparva.structure,"Data/STRUCTURE/Pparva.structure.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)






##############################################################################################################################################################################################################
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
#										CONVERGENCE ASSESSMENT OF STRUCTURE IN R
#!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!#
##############################################################################################################################################################################################################



##############################################################################################################################################################################################################
#
#										CONVERGENCE ASSESSMENT FOR NATIVE POPULATIONS
#                   & 100,000 ITERATIONS AFTER A BURNIN OF 100,000
#
##############################################################################################################################################################################################################

# STRUCTURE was preliminary run with 100,000 ITERATIONS AFTER A BURNIN OF 100,000 in order to identify a close range of most probable K
# at a reasonable computational time
# However, we didn't expected to reach convergence because of the reduced iterations and burnin

# We needed to import the trace of Alpha and posterior estimates of Ln P(D) (parameters to estimate and the sampling chain)
# One file per value of K

for (k in 1:9) {
  assign(paste("sampling.chain.K",k,sep=""), read.table(paste(STRdir,"/2b_Native/Results/sampling_chain_K",k,sep=""),header=T,sep=" "))
  # assign(paste("sampling.chain.K",k,sep=""), read.table(paste(STRdir,"/21_Native_1MillionIterations/Results/sampling_chain_K",k,sep=""),header=T,sep=" "))
}

nrep=20 # Number of replicates of the sampling chain


##########################################################################
#!# Checking Convergence Visually
##########################################################################

#========================================================================#
# Alpha -- Sampling chain after burnin
#------------------------------------------------------------------------#
# Trace plot of multiple replicates

nrep=10 # Reduce the number of replicates to increase clarity of trace plots
cols=c(brewer.pal(n = 12, name = "Paired"),brewer.pal(n = 8, name = "Set3")) # Set a vector of 20 colors for 20 replicates

#####################################################
# Assessing convergence for K=2
k=2 # K value to test
for (i in 1:k) {
  a=paste("Alpha",i,"_K",k,sep="") # Legend of Y axis
  sampling.chain=get(paste("sampling.chain.K",k,sep=""))
  test.list<-list()
  for (r in 1:nrep){ test.list[[r]]=as.mcmc(sampling.chain[sampling.chain$Run==r & !is.na(sampling.chain$Ln.Like),(3+i)]) }
  rm(r)
  # just the sampling chain, after burnin
  traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
  rm(a)
  rm(test.list)
}
# Most sampling chains didn't converge for parameter Alpha in K=2 and 20 replicates
# High variance within and between chains


#####################################################
# Assessing convergence for K=3
k=3 # K value to test
for (i in 1:k) {
  a=paste("Alpha",i,"_K",k,sep="") # Legend of Y axis
  sampling.chain=get(paste("sampling.chain.K",k,sep=""))
  test.list<-list()
  for (r in 1:nrep){ test.list[[r]]=as.mcmc(sampling.chain[sampling.chain$Run==r & !is.na(sampling.chain$Ln.Like),(3+i)]) }
  rm(r)
  # just the sampling chain, after burnin
  traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
  rm(a)
  rm(test.list)
}
# Most sampling chains didn't converge for parameter Alpha in K=3 and 20 replicates
# High variance within and between chains


#####################################################
# Assessing convergence for K=5
k=5 # K value to test
for (i in 1:k) {
  a=paste("Alpha",i,"_K",k,sep="") # Legend of Y axis
  sampling.chain=get(paste("sampling.chain.K",k,sep=""))
  test.list<-list()
  for (r in 1:nrep){ test.list[[r]]=as.mcmc(sampling.chain[sampling.chain$Run==r & !is.na(sampling.chain$Ln.Like),(3+i)]) }
  rm(r)
  # just the sampling chain, after burnin
  traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
  rm(a)
  rm(test.list)
}
# Some sampling chains seemed close to convergence for parameter Alpha in K=5 and 20 replicates
# But there were still high variance and different means inferred.
# Parallel chains didn't converged on the same value.

#####################################################
# Assessing convergence for K=6
k=6 # K value to test
for (i in 1:k) {
  a=paste("Alpha",i,"_K",k,sep="") # Legend of Y axis
  sampling.chain=get(paste("sampling.chain.K",k,sep=""))
  test.list<-list()
  for (r in 1:nrep){ test.list[[r]]=as.mcmc(sampling.chain[sampling.chain$Run==r & !is.na(sampling.chain$Ln.Like),(3+i)]) }
  rm(r)
  # just the sampling chain, after burnin
  traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
  rm(a)
  rm(test.list)
}
# Some sampling chains seemed close to convergence for parameter Alpha in K=6 and 20 replicates
# But there were still high variance and different means inferred.
# Parallel chains didn't converged on the same value.

#####################################################
# Assessing convergence for K=7
k=7 # K value to test
for (i in 1:k) {
  a=paste("Alpha",i,"_K",k,sep="") # Legend of Y axis
  sampling.chain=get(paste("sampling.chain.K",k,sep=""))
  test.list<-list()
  for (r in 1:nrep){ test.list[[r]]=as.mcmc(sampling.chain[sampling.chain$Run==r & !is.na(sampling.chain$Ln.Like),(3+i)]) }
  rm(r)
  # just the sampling chain, after burnin
  traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
  rm(a)
  rm(test.list)
}
# Some sampling chains seemed close to convergence for parameter Alpha in K=7 and 20 replicates
# But there were still high variance and different means inferred.
# Parallel chains didn't converged on the same value.

#!!!!!!!!!!!!!!!!#
# Alpha parameter globally converged in some samplign chains, while most others chains didn't showed signs of convergence
# High variance whithin and between chains, and no clear graphical assessment of convergence
# A pattern was distinctive: there was no convergence at all for K=2-3, but convergence seemed to appear with higher K values
# There were stationarity of most sampling chains for K=6-7, despite some unsteadiness and outliers chains

# We need to go further with more quantitative diagnostics...


#========================================================================#
# Ln Likelihood - Sampling chain after burnin
#------------------------------------------------------------------------#
png("Figures/STRUCTURE/2b_Native/Traceplot_Likelihood_convergence.png",width=3000, height = 2000, pointsize = 60)
par(mfrow=c(2,3))
cols=c(brewer.pal(n = 12, name = "Paired"),brewer.pal(n = 8, name = "Set3")) # Set a vector of 20 colors for 20 replicated sampling chains
#------------------------------------------------------------------------#
# Assessing convergence for K=2
a="Ln Likelihood (K2)" # Legend of Y axis
test.list<-list()
for (r in 1:nrep){ test.list[[r]]<-as.mcmc(sampling.chain.K2[sampling.chain.K2$Run==r & !is.na(sampling.chain.K2$Ln.Like),]$Ln.Lik) }
rm(r)
# just the sampling chain, after burnin
traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a,cex=4)
rm(a)
rm(test.list)

#------------------------------------------------------------------------#
# Assessing convergence for K=3
a="Ln Likelihood (K3)" # Legend of Y axis
test.list<-list()
for (r in 1:nrep){ test.list[[r]]<-as.mcmc(sampling.chain.K3[sampling.chain.K3$Run==r & !is.na(sampling.chain.K3$Ln.Like),]$Ln.Lik) }
rm(r)
# just the sampling chain, after burnin
traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
rm(a)
rm(test.list)

#------------------------------------------------------------------------#
# Assessing convergence for K=4
a="Ln Likelihood (K4)" # Legend of Y axis
test.list<-list()
for (r in 1:nrep){ test.list[[r]]<-as.mcmc(sampling.chain.K4[sampling.chain.K4$Run==r & !is.na(sampling.chain.K4$Ln.Like),]$Ln.Lik) }
rm(r)
# just the sampling chain, after burnin
traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
rm(a)
rm(test.list)

#------------------------------------------------------------------------#
# Assessing convergence for K=5
a="Ln Likelihood (K5)" # Legend of Y axis
test.list<-list()
for (r in 1:nrep){ test.list[[r]]<-as.mcmc(sampling.chain.K5[sampling.chain.K5$Run==r & !is.na(sampling.chain.K5$Ln.Like),]$Ln.Lik) }
rm(r)
# just the sampling chain, after burnin
traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
rm(a)
rm(test.list)

#------------------------------------------------------------------------#
# Assessing convergence for K=6
a="Ln Likelihood (K6)" # Legend of Y axis
test.list<-list()
for (r in 1:nrep){ test.list[[r]]<-as.mcmc(sampling.chain.K6[sampling.chain.K6$Run==r & !is.na(sampling.chain.K6$Ln.Like),]$Ln.Lik) }
rm(r)
# just the sampling chain, after burnin
traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
rm(a)
rm(test.list)

#------------------------------------------------------------------------#
# Assessing convergence for K=7
a="Ln Likelihood (K7)" # Legend of Y axis
test.list<-list()
for (r in 1:nrep){ test.list[[r]]<-as.mcmc(sampling.chain.K7[sampling.chain.K7$Run==r & !is.na(sampling.chain.K7$Ln.Like),]$Ln.Lik) }
rm(r)
# just the sampling chain, after burnin
traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
rm(a)
rm(test.list)
dev.off()
par(mfrow=c(1,1))

#  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  #
# Likelihood sampling chains seemed all convergent, graphically
# But only K=3 (range -53400;-53100) and K=7 (range -48200;-47600) showed convergence of all replicates to a single mode
# Others K were multimodal, with different means between replicates
# 2 modes for K=2 (range -57200;-56200), 3 modes for K=4-5 (range -52200;-50800 & -50700;-49200 respectively)
# and one mode with one outlier for K=6 (range -49200;-48400)


#------------------------------------------------------------------------#
# r
#------------------------------------------------------------------------#
# Assessing convergence for K=6
png("Figures/STRUCTURE/2b_Native/Traceplot_r_convergence.png",width=3000, height = 2000, pointsize = 60)
par(mfrow=c(2,3))
cols=c(brewer.pal(n = 12, name = "Paired"),brewer.pal(n = 8, name = "Set3")) # Set a vector of colors
#-----------------------------------
a="r (K2)" # Legend of Y axis
test.list<-list()
for (r in 1:nrep){ test.list[[r]]<-as.mcmc(sampling.chain.K2[sampling.chain.K2$Run==r & !is.na(sampling.chain.K2$Ln.Like),]$r) }
rm(r)
# just the sampling chain, after burnin
traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
rm(a)
rm(test.list)
#-----------------------------------
a="r (K3)" # Legend of Y axis
test.list<-list()
for (r in 1:nrep){ test.list[[r]]<-as.mcmc(sampling.chain.K3[sampling.chain.K3$Run==r & !is.na(sampling.chain.K3$Ln.Like),]$r) }
rm(r)
# just the sampling chain, after burnin
traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
rm(a)
rm(test.list)
#-----------------------------------
a="r (K4)" # Legend of Y axis
test.list<-list()
for (r in 1:nrep){ test.list[[r]]<-as.mcmc(sampling.chain.K4[sampling.chain.K4$Run==r & !is.na(sampling.chain.K4$Ln.Like),]$r) }
rm(r)
# just the sampling chain, after burnin
traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
rm(a)
rm(test.list)
#-----------------------------------
a="r (K5)" # Legend of Y axis
test.list<-list()
for (r in 1:nrep){ test.list[[r]]<-as.mcmc(sampling.chain.K5[sampling.chain.K5$Run==r & !is.na(sampling.chain.K5$Ln.Like),]$r) }
rm(r)
# just the sampling chain, after burnin
traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
rm(a)
rm(test.list)
#-----------------------------------
a="r (K6)" # Legend of Y axis
test.list<-list()
for (r in 1:nrep){ test.list[[r]]<-as.mcmc(sampling.chain.K6[sampling.chain.K2$Run==r & !is.na(sampling.chain.K6$Ln.Like),]$r) }
rm(r)
# just the sampling chain, after burnin
traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
rm(a)
rm(test.list)
#-----------------------------------
a="r (K7)" # Legend of Y axis
test.list<-list()
for (r in 1:nrep){ test.list[[r]]<-as.mcmc(sampling.chain.K7[sampling.chain.K7$Run==r & !is.na(sampling.chain.K7$Ln.Like),]$r) }
rm(r)
# just the sampling chain, after burnin
traceplot(test.list,lty=1, col=cols,xlab="Iterations (x100)", ylab=a)
rm(a)
rm(test.list)
dev.off()
par(mfrow=c(1,1))

#  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  #
# The r parameter didn't seemed convergent, with a high variance within and between chains


##########################################################################
# Testing Within Chain Convergence & Stability
##########################################################################

#========================================================================#
# Basic diagnostics of MCMC sampling chain & convergence assessment
#------------------------------------------------------------------------#

# Basic statistics computed on the sampling chains
nrep=20
k=3
sampling.chain=get(paste("sampling.chain.K",k,sep=""))

# Alpha chains
for (rep in 1:nrep) { # With rep the index of the replicate
  cat("----------------------------------------------\nReplicate number",rep,"\n")
  print(summary(as.mcmc(sampling.chain.K6[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),c(4:9)])))
}
# Ln Likelihood chain
for (rep in 1:nrep) { # With rep the index of the replicate
  cat("----------------------------------------------\nReplicate number",rep,"\n")
  print(summary(as.mcmc(sampling.chain.K6[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),ncol(sampling.chain)-1])))
}# r parameter chain
for (rep in 1:nrep) { # With rep the index of the replicate
  cat("----------------------------------------------\nReplicate number",rep,"\n")
  print(summary(as.mcmc(sampling.chain.K6[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),(3+2*k+1)])))
}


#------------------------------------------------------------------------#
# Trace & density
#------------------------------------------------------------------------#

rep=1 # Index of replicate to plot and summarize
########### Alphas
par(mfrow=c(3,2))
a.draws=as.mcmc(sampling.chain.K6[which(sampling.chain.K6$Run==rep & !(is.na(sampling.chain.K6$Ln.Like))),4])
summary(a.draws)
plot(a.draws)
a.draws=as.mcmc(sampling.chain.K6[which(sampling.chain.K6$Run==rep & !(is.na(sampling.chain.K6$Ln.Like))),5])
summary(a.draws)
plot(a.draws)
a.draws=as.mcmc(sampling.chain.K6[which(sampling.chain.K6$Run==rep & !(is.na(sampling.chain.K6$Ln.Like))),6])
summary(a.draws)
plot(a.draws)
a.draws=as.mcmc(sampling.chain.K6[which(sampling.chain.K6$Run==rep & !(is.na(sampling.chain.K6$Ln.Like))),7])
summary(a.draws)
plot(a.draws)
a.draws=as.mcmc(sampling.chain.K6[which(sampling.chain.K6$Run==rep & !(is.na(sampling.chain.K6$Ln.Like))),8])
summary(a.draws)
plot(a.draws)
a.draws=as.mcmc(sampling.chain.K6[which(sampling.chain.K6$Run==rep & !(is.na(sampling.chain.K6$Ln.Like))),9])
summary(a.draws)
plot(a.draws)
par(mfrow=c(1,1))

########### Ln Likelihood
LnLk.draws=as.mcmc(sampling.chain.K6[which(sampling.chain.K6$Run==rep & !(is.na(sampling.chain.K6$Ln.Like))),ncol(sampling.chain)-1])
summary(LnLk.draws)
plot(LnLk.draws)

########## Parameter r
r.draws=as.mcmc(sampling.chain.K6[which(sampling.chain.K6$Run=rep & !(is.na(sampling.chain.K6$Ln.Like))),(3+2*k+1)])
summary(r.draws)
plot(r.draws)

rm(rep)

#------------------------------------------------------------------------#
# Autocorrelation of Alpha in MCMC
#------------------------------------------------------------------------#



#------------------------------------------------------------------------#
# Gelman-Rubin diagnostic for Alpha
#------------------------------------------------------------------------#




#========================================================================#
### TESTS OF CONVERGENCE
#========================================================================#
# Subsequently, we ran a battery of diagnostic tests of convergence on the sampling chains of Alpha, r and Ln Likelihood of K

#========================================================================#
# Geweke Diagnostic of MC Chain Stability
#------------------------------------------------------------------------#
# if the whole chain is stationary, the means of the values early and late in the sequence should be similar
# convergence diagnostic 'Z' is the difference between the 2 means divided by the asymptotic standard error of their difference
# values of 'Z' near the extreme tails of the N(0,1) indicates lack of convergence
# can also estimate p-value of 'Z' from the normal distribution
# yields the probability that the divided chain means are different
?geweke.diag


#------------------------------------------------------------------------#
# Alpha
#------------------------------------------------------------------------#
# We investigated K=2-3-5-6-7 for 1 replicate and all alpha values at a time
# Remind that it were 100,000 burnin + 100,000 iterations
# We had to discard the burnin by not considering iterations with NA values in Ln.likelihood

# First step: assess convergence
nrep=20
k=5
sampling.chain=get(paste("sampling.chain.K",k,sep=""))
convergent=c()
for (rep in 1:nrep) {
  cat("----------------------------------------------\nReplicate number",rep,"\n")
  print(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]), frac1=0.1, frac2=0.5))
  # A p-value can be computed for this Z-score
  print(2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]), frac1=0.1, frac2=0.5)$z)))
  # count the number of convergent chains in the replicate (i.e. p-value > 0.05)
  print(sum(2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]), frac1=0.1, frac2=0.5)$z))>0.05))
  convergent[rep]=sum(2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]), frac1=0.1, frac2=0.5)$z))>0.05)
  # As well as a diagnostic plot
  geweke.plot(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]), frac1=0.1, frac2=0.5, nbins=100)
}
print(convergent)
sum(convergent==k) # number of totally convergent replicates
sum(convergent>=(k/2) & convergent!=k) # number of partially convergent replicates (at least half of parameters)
rm(convergent)


##########################
#   IMPORTANT INFORMATION
# We need a list of convergent replicates, so we can remove all non-convergent replicates for subsequent analyses (Harvester, CLUMPP)



#------------------------
# K=2
# 10 replicates were fully convergent for Alpha, and 8 more replicates were partially convergent (1/2)

#------------------------
# K=3
# 6 replicates were fully convergent for Alpha, and 8 more replicates were partially convergent (2/3)

#------------------------
# K=5
# Only 3 replicates were fully convergent, whereas only 3 more were partially convergent (>= 3/5)

#------------------------
# K=6
# But it was partially convergent for 9 replicates
# and totally convergent for 4 replicates
# Hence, the convergence is not great on this K value for Alphas

#------------------------
# K=7
# 4 replicates were fully convergent for Alpha, and 6 more replicates were partially convergent (>= 4/6)


#------------------------------------------------------------------------#
# Ln Likelihood
#------------------------------------------------------------------------#

# First step: assess convergence

nrep=20
k=7
sampling.chain=get(paste("sampling.chain.K",k,sep=""))
convergent=c()
for (rep in 1:nrep) {
  cat("----------------------------------------------\nReplicate number",rep,"\n")
  print(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),ncol(sampling.chain)-1]), frac1=0.1, frac2=0.5))
  # A p-value can be computed for this Z-score
  print(2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),ncol(sampling.chain)-1]), frac1=0.1, frac2=0.5)$z)))
  # count the number of convergent chains in the replicate (i.e. p-value > 0.05)
  print(sum(2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),ncol(sampling.chain)-1]), frac1=0.1, frac2=0.5)$z))>0.05))
  convergent[rep]=sum(2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),ncol(sampling.chain)-1]), frac1=0.1, frac2=0.5)$z))>0.05)
  # As well as a diagnostic plot
  geweke.plot(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),ncol(sampling.chain)-1]), frac1=0.1, frac2=0.5)
}
print(convergent)
sum(convergent==1) # number of convergent replicates
rm(convergent)

#------------------------
# K=2
# There was only 9 replicates that achieved convergence for likelihood on K=2

#------------------------
# K=3
# There was only 8 replicates that achieved convergence for likelihood on K=3

#------------------------
# K=5
# There was 13 replicates that achieved convergence for likelihood on K=5

#------------------------
# K=6
# There was only 8 replicates that achieved convergence for likelihood on K=6


#------------------------
# K=7
# There was 11 replicates that achieved convergence for likelihood on K=7

# It is difficult to choose a K from these results





#------------------------------------------------------------------------#
# Parameter r
#------------------------------------------------------------------------#
(3+2*k+1)

nrep=20
k=7
sampling.chain=get(paste("sampling.chain.K",k,sep=""))
convergent=c()
for (rep in 1:nrep) {
  cat("----------------------------------------------\nReplicate number",rep,"\n")
  print(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),(3+2*k+1)]), frac1=0.1, frac2=0.5))
  # A p-value can be computed for this Z-score
  print(2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),(3+2*k+1)]), frac1=0.1, frac2=0.5)$z)))
  # count the number of convergent chains in the replicate (i.e. p-value > 0.05)
  print(sum(2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),(3+2*k+1)]), frac1=0.1, frac2=0.5)$z))>0.05))
  convergent[rep]=sum(2*pnorm(-abs(geweke.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),(3+2*k+1)]), frac1=0.1, frac2=0.5)$z))>0.05)
  # As well as a diagnostic plot
  geweke.plot(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),(3+2*k+1)]), frac1=0.1, frac2=0.5)
}
print(convergent)
sum(convergent==1) # number of convergent replicates
rm(convergent)

#------------------------
# K=2
# There was only 9 replicates that achieved convergence for r on K=2

#------------------------
# K=3
# There was only 6 replicates that achieved convergence for r on K=3

#------------------------
# K=5
# There was 5 replicates that achieved convergence for r on K=5

#------------------------
# K=6
# There was only 7 replicates that achieved convergence for r on K=6

#------------------------
# K=7
# There was 9 replicates that achieved convergence for r on K=7


#========================================================================#
# Heidelberger and Welch's Convergence Diagnostic
#------------------------------------------------------------------------#
# probability of rejecting hypothesis that Markov Chain is a stable/stationary distribution
# if Halfwidth test fails, chain should be extended
?heidel.diag

#------------------------------------------------------------------------#
# Alpha
#------------------------------------------------------------------------#
nrep=20
k=2
sampling.chain=get(paste("sampling.chain.K",k,sep=""))
for (rep in 1:nrep) {
  cat("----------------------------------------------\nReplicate number",rep,"\n")
  print(heidel.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))])))
}

# Convergence of Alpha was unclear, with high variability between replicates

#------------------------------------------------------------------------#
# Ln Likelihood
#------------------------------------------------------------------------#
nrep=20
k=7
sampling.chain=get(paste("sampling.chain.K",k,sep=""))
for (rep in 1:nrep) {
  cat("----------------------------------------------\nReplicate number",rep,"\n")
  print(heidel.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),ncol(sampling.chain)-1])))
}

#------------------------
# K=2
# 18 replicates passed the stationarity test

#------------------------
# K=3
# 14 replicates passed the stationarity test

#------------------------
# K=5
# 18 replicates passed the stationarity test

#------------------------
# K=6
# 14 replicates passed the stationarity test

#------------------------
# K=7
# 17 replicates passed the stationarity test



#------------------------------------------------------------------------#
# Parameter r
#------------------------------------------------------------------------#
nrep=20
k=2
sampling.chain=get(paste("sampling.chain.K",k,sep=""))
for (rep in 1:nrep) {
  cat("----------------------------------------------\nReplicate number",rep,"\n")
  print(heidel.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep & !(is.na(sampling.chain$Ln.Like))),(3+2*k+1)])))
}



#========================================================================#
# Raftery and Lewis' Diagnostic
#------------------------------------------------------------------------#
# run length control diagnostic
# intended for use on a short pilot run of a Markov chain
# to estimate number of iterations required to estimate the quantile 'q' to within an accuracy of +/- 'r' with probability 'p'
# only tests marginal convergence on each parameter
#!# high dependence factors (i.e. >5) are worrisome
# may indicate influential starting values, high correlations between coefficients, or poor mixing
?raftery.diag

#------------------------------------------------------------------------#
# Alpha
#------------------------------------------------------------------------#
nrep=20
k=2
sampling.chain=get(paste("sampling.chain.K",k,sep=""))
for (rep in 1:nrep) {
  cat("----------------------------------------------\nReplicate number",rep,"\n")
  #!!!!!!!! q=0.05 (default 0.025) and precision at 5% (0.05) to get results: instead sample size was not enough
  print(raftery.diag(as.mcmc(sampling.chain[which(sampling.chain$Run==rep),c(4:(3+k))]), q=0.025, r=0.05, s=0.95))
}

# Dependence factor (I) was well beyond 5, so we had to worry about autocorrelation in the sampling chain 
# BUT PROBLEM WITH UNDERSTANDING OF THIS DIAGNOSTIC
# Thinning has lareaday been done (thinning of 100): what influence for this pilot run?
# Analyses in standby until more comprehension...

#------------------------------------------------------------------------#
# Ln Likelihood
#------------------------------------------------------------------------#

#------------------------------------------------------------------------#
# Parameter r
#------------------------------------------------------------------------#


#------------------------------------------------------------------------#
# Summarizing Run Length Inference
#------------------------------------------------------------------------#
r=1
chain.length=data.frame(k=g, Rep=r,
                         Alpha.N=raftery.diag(as.mcmc(subset(trace.files, k==g & Rep==r)$Alpha), q=0.025, r=0.005, s=0.95)[[2]][,"N"],
                         Alpha.DepFactor=raftery.diag(as.mcmc(subset(trace.files, k==g & Rep==r)$Alpha), q=0.025, r=0.005, s=0.95)[[2]][,"I"])
row.names(chain.length)<-NULL
chain.length

for (r in 2:5)
{ temp.df<-data.frame(k=g, Rep=r,
                      Alpha.N=raftery.diag(as.mcmc(subset(trace.files, k==g & Rep==r)$Alpha), q=0.025, r=0.005, s=0.95)[[2]][,"N"],
                      Alpha.DepFactor=raftery.diag(as.mcmc(subset(trace.files, k==g & Rep==r)$Alpha), q=0.025, r=0.005, s=0.95)[[2]][,"I"])
row.names(temp.df)<-NULL
chain.length<-rbind(chain.length, temp.df)
rm(temp.df)
}
chain.length
rm(g)
rm(r)

for (g in 2:7)
{ for (r in 1:5)
{ temp.df<-data.frame(k=g, Rep=r,
                      Alpha.N=raftery.diag(as.mcmc(subset(trace.files, k==g & Rep==r)$Alpha), q=0.025, r=0.005, s=0.95)[[2]][,"N"],
                      Alpha.DepFactor=raftery.diag(as.mcmc(subset(trace.files, k==g & Rep==r)$Alpha), q=0.025, r=0.005, s=0.95)[[2]][,"I"])
row.names(temp.df)<-NULL
chain.length<-rbind(chain.length, temp.df)
rm(temp.df)
}
  rm(r)
}
rm(g)
chain.length
summary(chain.length$Alpha.N)

boxplot(Alpha.N~k, data=chain.length, xlab="k", ylab="Min. Chain Length Required for Convergence", log="y", ylim=c(20000,500000), cex.lab=1.5)

##########################################################################
#!# Testing Convergence Among Chains
##########################################################################


#========================================================================#
# Gelman diagnostic: 'potential scale reduction factor'
#------------------------------------------------------------------------#
?gelman.diag
#------------------------------------------------------------------------#
# Alphas
#------------------------------------------------------------------------#
# testing convergence among chains
#!# scaling factor should be less than 1.2
# values substantially greater than 1 indicate lack of convergence
# No autoburnin, as burnin had already been removed
nrep=20
k=2
sampling.chain=get(paste("sampling.chain.K",k,sep=""))
gelman.diag(mcmc.list(list(
  as.mcmc(sampling.chain[which(sampling.chain$Run==1 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==2 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==3 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==4 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==5 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==6 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==7 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==8 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==9 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==10 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==11 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==12 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==13 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==14 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==15 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==16 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==17 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==18 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==19 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==20 & !(is.na(sampling.chain$Ln.Like))),c(4:(3+k))])
)),autoburnin=FALSE)

#------------------------
# K=2
# Alpha 1 CI (1.56-1.89) and Alpha 2 (1.70-2.09)
# Approximate convergence is diagnosed when the upper limit is close to 1
# Thus, for K=2, convergence is not achieved
#------------------------
# K=3
#         Point est. Upper C.I.
# Alpha1       1.92       2.40
# Alpha2       1.97       2.46
# Alpha3       1.44       1.71
#------------------------
# K=5
#         Point est. Upper C.I.
# Alpha1       1.16       1.29
# Alpha2       1.32       1.54
# Alpha3       1.32       1.54
# Alpha4       1.49       1.79
# Alpha5       1.31       1.52
#------------------------
# K=6
#         Point est. Upper C.I.
# Alpha1       1.21       1.35
# Alpha2       1.30       1.53
# Alpha3       1.37       1.61
# Alpha4       1.37       1.61
# Alpha5       1.41       1.70
# Alpha6       1.36       1.62
#------------------------
# K=7
#         Point est. Upper C.I.
# Alpha1       1.25       1.42
# Alpha2       1.23       1.39
# Alpha3       1.23       1.39
# Alpha4       1.23       1.39
# Alpha5       1.27       1.45
# Alpha6       1.28       1.49
# Alpha7       1.19       1.33

# Globally, Gelman and Rubin diagnostic didn't supported the approximate convergence for Alpha,
# despite the fact that mean upper limits were around 1.5 for K=5-6-7
# Alpha parameter was better estimated with higher K values (K=5-6-7) than the K=3 given by Delta(K) method

#------------------------------------------------------------------------#
# Ln Likelihood
#------------------------------------------------------------------------#
k=7
sampling.chain=get(paste("sampling.chain.K",k,sep=""))
gelman.diag(mcmc.list(list(
  as.mcmc(sampling.chain[which(sampling.chain$Run==1 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==2 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==3 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==4 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==5 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==6 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==7 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==8 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==9 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==10 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==11 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==12 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==13 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==14 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==15 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==16 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==17 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==18 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==19 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)]),
  as.mcmc(sampling.chain[which(sampling.chain$Run==20 & !(is.na(sampling.chain$Ln.Like))),(ncol(sampling.chain)-1)])
)),autoburnin=FALSE)

#------------------------
# K=2
# Point est. Upper C.I.
# 8.66       11.4
#------------------------
# K=3
# 1.04       1.06
#------------------------
# K=5
# 5.84       7.64
#------------------------
# K=6
# 2.22       2.78
#------------------------
# K=7
# 2.21       2.79

# According to Gelman and Rubin's diagnostic, Ln Likelihood were convergent for K=3
# and close but not convergent for K=6-7 (upper limit 2.8)




# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# Using what we've learned fro preliminary runs, we re-launched the analysis with new parameters and config
# Number of burnin: 500,000
# Number of iterations: 200,000


##############################################################################################################################################################################################################
#
#										CONVERGENCE ASSESSMENT FOR NATIVE POPULATIONS
#                   & 200,000 ITERATIONS AFTER A BURNIN OF 500,000
#
##############################################################################################################################################################################################################



#=================================
# THE END
#=================================
