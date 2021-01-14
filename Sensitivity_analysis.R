############################################################################
#   INRA - ANR project Pseudorasbora
#
#       Sensitivity analysis of the SNPs data set for ABC
#
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com
# MSc internship

# Supervisor: Scott McCairns
#             INRA UMR ESE


# Huge large-scale data set of SNP are computationnaly demanding, as well as Bayesian approaches.
# Thus, a reduction of the original data set is often recommended and performed before Bayesian analysis.
# Hence, we needed first to assess the sensitivity of ABC approaches to our data set and look at which extent SuSt of ABC were affected by
# subsetting inside the data set.
# Different paradigms of subsetting:
# - Random sampling
# - Quality sampling, of the most informative markers (threshold of missing data proportion)


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
library(ggplot2)


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

# Load Pparva object, global data set before cleaning and subsetting
load(paste(datadir,file="/Pparva.3000.Rda",sep=""))
load(paste(datadir,file="/Pparva.Rda",sep=""))

# Manipulation of a hierfstat object is easier than a genind
Pparva.hier=genind2hierfstat(Pparva)

###########################################################
#
#       Sensitivity analysis for each independent SuSt
#
###########################################################

# Sensitivity analysis is the study of how the uncertainty in the output of a mathematical model or system (numerical or otherwise) can be apportioned to different sources of uncertainty in its inputs
# Here, we tested sensitivity to missing data and data reduction (subsetting)

#==========================================================
# RANDOM SET OF MARKERS
#==========================================================
# Basic idea was here to study the evolution of SuSt (summary statistics)
# as a function of a number of markers randomly sampled inside the global data set

# Usual SuSt for ABC are (Jeffries et al. 2016):
# - mean genetic diversity across all polymorphic loci 
# - variance of genetic diversity across all polymorphic loci
# - mean genetic diversity across all loci                                                ---> DONE with Hst
# - variance of genetic diversity across all loci
# - mean/variance of Fst (for Fst > 0 for 2 loci & for all loci)
# - mean/variance of Nei's distance (for Nei's distance > 0 for 2 loci & for all loci)
# - max-likelihood estimates of admixture proportions

# Consult the user manual of DIYABC for more detailed information as above: 12 SuSt
# Single sample statistics :
# - proportion of loci with null gene diversity (= proportion of monomorphic loci)        ---> DONE with Hs
# - mean gene diversity across polymorphic loci (Nei, 1987)                               ---> WAITING
# - variance of gene diversity across polymorphic loci                                    ----> WAITING
# - mean gene diversity across all loci                                                   ---> DONE with Hst
# Two sample statistics :
# - proportion of loci with null FST distance between the two samples (Weir and Cockerham, 1984)
# - mean across loci of non null FST distances between the two samples
# - variance across loci of non null FST distances between the two samples
# - mean across loci of FST distances between the two samples
# - proportion of loci with null Nei's distance between the two samples (Nei, 1972)
# - mean across loci of non null Nei's distances between the two samples
# - variance across loci of non null Nei's distances between the two samples
# - mean across loci of Nei's distances between the two samples


#=====================================================================================================
# STUDY OF A SINGLE SU. ST.
#-----------------------------------------------------------------------------------------------------


#-----------------------------------------------------------------------------------------------------
# MEAN GENETIC DIVERSITY ACROSS ALL LOCI: Hst
#-----------------------------------------------------------------------------------------------------
# We chose to investigate first: mean genetic diversity across all loci, the Hst stat. of hierfstat's function basic.stats() (Nei, 1987)

mean(basic.stats(Pparva.hier)$overall[3])
# Confidence interval (95%) of the "true" Hst on the global data set was computed with 1,000 bootstraps (resampling of individuals)
# BOOTSTRAP MEAN Hst
boots=numeric(1000)
for (i in 1:1000){
  print(i)
  boots[i]=basic.stats(Pparva.hier[,c(1,sample(2:(ncol(Pparva.hier)),replace=T))])$overall[3]
}
save(boots,file="Data/bootstrap.mean.Hst")
load(file="Data/bootstrap.mean.Hst")
hist(boots)
limBoot=quantile(boots,c(.025,.975))
(IC.Hst=limBoot)
(mean.Hst=mean(boots))
(sd.Hst=sd(boots))

# Build a matrix with n iterations for each value of N (nb of randomly sampled SNPs), with iterations in columns and SuSt in rows
niter=100
samp.n=seq(100,(ncol(Pparva.hier)-101),100) # Number of markers to sample
SuSt.Hst=matrix(NA,nrow=niter,ncol=length(samp.n)) # Index of each row correspond to one randomised replicate of the SuSt
# Index of the column give the number of markers sampled: (N -101)/100 markers sampled (-1 for first column of pop & -100 for setting lowest number of markers to 100)
# Indices are a sequence from 100 to N markers by step of 100 markers


# From 1 sampled SNP only to all SNPs sampled
c=0
for (N in samp.n) {
  c=c+1
  # The SuSt is computed from a random subset of N markers for each iteration
  for (i in 1:niter) {
    cat("Iteration:", i,"of",N,"markers sampled (",c,"/",length(samp.n),")\n")
    # Subset N marker randomly at each iteration (subset N columns in genind@tab)
    subset=Pparva.hier[,c(1,sample(2:ncol(Pparva.hier),N,replace=FALSE))]
    # Compute the SuSt, and add it to the matrix 'SuSt'
    SuSt.Hst[i,c]=basic.stats(subset)$overall[3]
  }
}
rm(c)
save(SuSt.Hst,file="Data/SuSt_Hst.Rda")
load("Data/SuSt_Hst.Rda")

# !!! #  Computing time was very long and time increased exponentially with number of markers,
# !!! #  so we needed to implement parallel programmming to achieve a large number of random bootstraps

# Plotting the mean +- sd of the stat for each number of randomly sampled markers
# Add a grey ribbon for the confidence interval (95%) of the "true" Hst on the global data set, computed with 1,000 bootstraps
# df=data.frame(N=samp.n, Mean=apply(SuSt.Hst, 2, mean), sd=apply(SuSt.Hst, 2, sd))

# Use CI instead of the SD of the resampled mean?
# suppose a symetric distribution (one value of CI only, the range between the mean and the limit at 0.975/0.025)
df=data.frame(N=samp.n, Mean=apply(SuSt.Hst, 2, mean), sd=apply(SuSt.Hst, 2,function(x) quantile(x,0.975))-apply(SuSt.Hst, 2, mean)) # Plot the CI instead of sd of the resampled means

plot.Hst=ggplot(data=df, aes(x=N, y=Mean)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2) +
    geom_hline(yintercept=c(IC.Hst[1], IC.Hst[2]), linetype="dashed", size=1) +
    xlab("Number of markers randomly sampled") + ylab("Mean genetic diversity across all loci (± sd)") +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
          plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
          axis.title.x = element_text(color="black", size=14),
          axis.title.y = element_text(color="black", size=14),
          axis.text=element_text(size=14, colour="black"),
          legend.key = element_rect(fill = "white", size = 1),
          legend.text=element_text(size=14),
          legend.title=element_text(size=14))
plot.Hst
ggsave(paste(figuresdir,"/Sensitivity/plot.Hst.png",sep=""),
       device="png",dpi=320,units="cm",width=20,height=16)

###########
# PROPORTION OF REPLICATES INSIDE THE CONFIDENCE INTERVAL OF THE TRUE VALUE OF THE SU. STATISTIC
# TYPE I ERROR: rejection of a true null hypothesis (false positive), i.e. seing a difference between the original data set and a subset

# For each number of marker, compute the proportion of values inside the CI of the true value
# A matrix of 2 columns with 1/ number of markers and 2/ proportion of exact values
true.Hst=matrix(NA,ncol=2,nrow=length(samp.n))
true.Hst[,1]=samp.n
true.Hst[,2]=apply(SuSt.Hst, 2, function(x) (sum(x>IC.Hst[1] & x<IC.Hst[2])/length(x)))
colnames(true.Hst)=c("N.Markers","Prop.Exact.Values")
true.Hst
# Plotting the proportion of values inside the CI of the true value
# Add a horizontal line at 95%
plot.prop.Hst=ggplot(data=as.data.frame(true.Hst), aes(x=N.Markers, y=Prop.Exact.Values)) +
  geom_point(colour="Black",size=1)+
  geom_line() +
  geom_hline(yintercept=0.95, linetype="dashed", size=1) +
  xlab("Number of markers randomly sampled") + ylab("Proportion of values inside the CI (95%) of the true value\nfor mean genetic diversity across all loci") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
plot.prop.Hst
ggsave(paste(figuresdir,"/Sensitivity/plot.trueValues.Hst.png",sep=""),
       device="png",dpi=320,units="cm",width=20,height=16)


#-----------------------------------------------------------------------------------------------------
# MEAN GENETIC DIVERSITY ACROSS POLYMORPHIC LOCI: Hst_poly
#-----------------------------------------------------------------------------------------------------
# We chose to investigate first: mean genetic diversity across polymorphic loci only, the Hst stat. of hierfstat's function basic.stats() (Nei, 1987)

# Confidence interval (95%) of the "true" Hst on the global data set was computed with 1,000 bootstraps
# BOOTSTRAP MEAN Hst
boots=numeric(100)
for (i in 1:100){
  print(i)
  temp=c()
  temp=basic.stats(Pparva.hier[,c(1,sample(2:(ncol(Pparva.hier)),replace=T))])$Hs
  boots[i]=mean(apply(temp,2,function(x) mean(x[x>0],na.rm=TRUE)))
  rm(temp)
}
save(boots,file="Data/bootstrap.mean.Hst_poly")
load(file="Data/bootstrap.mean.Hst_poly")
hist(boots)
limBoot=quantile(boots,c(.025,.975))
(IC.Hst=limBoot)
(mean.Hst=mean(boots))
(sd.Hst=sd(boots))

# Build a matrix with n iterations for each value of N (nb of randomly sampled SNPs), with iterations in columns and SuSt in rows
niter=100
samp.n=seq(100,(ncol(Pparva.hier)-101),100) # Number of markers to sample
SuSt.Hst.poly=matrix(NA,nrow=niter,ncol=length(samp.n)) # Index of each row correspond to one randomised replicate of the SuSt
# Index of the column give the number of markers sampled: (N -101)/100 markers sampled (-1 for first column of pop & -100 for setting lowest number of markers to 100)
# Indices are a sequence from 100 to N markers by step of 100 markers


# From 1 sampled SNP only to all SNPs sampled
c=0
for (N in samp.n) {
  c=c+1
  # The SuSt is computed from a random subset of N markers for each iteration
  for (i in 1:niter) {
    cat("Iteration:", i,"of",N,"markers sampled (",c,"/",length(samp.n),")\n")
    temp=c()
    # Subset N marker randomly at each iteration (subset N columns in genind@tab)
    subset=Pparva.hier[,c(1,sample(2:ncol(Pparva.hier),N,replace=FALSE))]
    # Compute the SuSt on all loci
    temp=basic.stats(subset)$Hs
    # remove monomorphic loci (0 values)
    SuSt.Hst.poly[i,c]=mean(apply(temp,2,function(x) mean(x[x>0],na.rm=TRUE)))
    rm(temp)
  }
}
rm(c)
save(SuSt.Hst.poly,file="Data/SuSt_Hst_poly.Rda")
load("Data/SuSt_Hst_poly.Rda")

# !!! #  Computing time was very long and time increased exponentially with number of markers,
# !!! #  so we needed to implement parallel programmming to achieve a large number of random bootstraps

# Plotting the mean +- sd of the stat for each number of randomly sampled markers
# Add a grey ribbon for the confidence interval (95%) of the "true" Hst on the global data set, computed with 1,000 bootstraps
df=data.frame(N=samp.n, Mean=apply(SuSt.Hst.poly, 2, mean), sd=apply(SuSt.Hst.poly, 2, sd))
plot.Hst=ggplot(data=df, aes(x=N, y=Mean)) +
  geom_point(colour="Black",size=1)+
  geom_line() +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2) +
  geom_hline(yintercept=c(IC.Hst[1], IC.Hst[2]), linetype="dashed", size=1) +
  xlab("Number of markers randomly sampled") + ylab("Mean genetic diversity across polymorphic loci (± sd)") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
plot.Hst
ggsave(paste(figuresdir,"/Sensitivity/plot.Hst.poly.png",sep=""),
       device="png",dpi=320,units="cm",width=20,height=16)

###########
# PROPORTION OF REPLICATES INSIDE THE CONFIDENCE INTERVAL OF THE TRUE VALUE OF THE SU. STATISTIC
# TYPE I ERROR: rejection of a true null hypothesis (false positive), i.e. seing a difference between the original data set and a subset

# For each number of marker, compute the proportion of values inside the CI of the true value
# A matrix of 2 columns with 1/ number of markers and 2/ proportion of exact values
true.Hst.poly=matrix(NA,ncol=2,nrow=length(samp.n))
true.Hst.poly[,1]=samp.n
true.Hst.poly[,2]=apply(SuSt.Hst.poly, 2, function(x) (sum(x>IC.Hst[1] & x<IC.Hst[2])/length(x)))
colnames(true.Hst.poly)=c("N.Markers","Prop.Exact.Values")
true.Hst.poly
# Plotting the proportion of values inside the CI of the true value
# Add a horizontal line at 95%
plot.prop.Hst=ggplot(data=as.data.frame(true.Hst.poly), aes(x=N.Markers, y=Prop.Exact.Values)) +
  geom_point(colour="Black",size=1)+
  geom_line() +
  geom_hline(yintercept=0.95, linetype="dashed", size=1) +
  xlab("Number of markers randomly sampled") + ylab("Proportion of values inside the CI (95%) of the true value\nfor mean genetic diversity across polymorphic loci") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
plot.prop.Hst
ggsave(paste(figuresdir,"/Sensitivity/plot.trueValues.Hst_poly.png",sep=""),
       device="png",dpi=320,units="cm",width=20,height=16)



#-----------------------------------------------------------------------------------------------------
# PROPORTION OF LOCI WITH NULL GENE DIVERSITY (I.E. MONOMORPHIC LOCI)
#-----------------------------------------------------------------------------------------------------
# Then we chose to investigate proportion of loci with null gene diversity
# We computed the mean proportion

# Confidence interval (95%) of the "true" proportion on the global data set was computed with 1,000 bootstraps
# BOOTSTRAP PROPORTION OF LOCI WITH NULL GENE DIVERSITY
(gene.div=basic.stats(Pparva.hier)$Hs) # bootstrap on populations
boots=numeric(1000)
for (i in 1:1000){
  print(i)
  # We first computed the bootstrap by resampling Hs values
  #boots[i]=mean(sample(apply(gene.div, 2, mean, na.rm=TRUE),replace=TRUE))
  # BUT it was a mistake and we computed next by resampling markers randomly
  # This way, it corresponded to the mean value at N=3999
  boots[i]=mean(basic.stats(Pparva.hier[,c(1,sample(2:(ncol(Pparva.hier)),replace=T))])$Hs,na.rm=TRUE)
}
save(boots,file="Data/bootstrap.mean.Hs")
load(file="Data/bootstrap.mean.Hs")
hist(boots)
(IC.Hs=quantile(boots,c(.025,.975)))
(mean.Hs=mean(boots))
(sd.Hs=sd(boots))

# Build a matrix with n iterations for each value of N (nb of randomly sampled SNPs), with iterations in columns and SuSt in rows
niter=100
samp.n=seq(100,(ncol(Pparva.hier)-101),100) # Number of markers to sample
SuSt.Hs=matrix(NA,nrow=niter,ncol=length(samp.n)) # Index of each row correspond to one randomised replicate of the SuSt
# Index of the column give the number of markers sampled: (N -101)/100 markers sampled (-1 for first column of pop & -100 for setting lowest number of markers to 100)
# Indices are a sequence from 100 to N markers by step of 100 markers


# From 1 sampled SNP only to all SNPs sampled
c=0
for (N in samp.n) {
  c=c+1
  # The SuSt is computed from a random subset of N markers for each iteration
  for (i in 1:niter) {
    cat("Iteration:", i,"of",N,"markers sampled (",c,"/",length(samp.n),")\n")
    # Subset N marker randomly at each iteration (subset N columns in genind@tab)
    subset=Pparva.hier[,c(1,sample(2:ncol(Pparva.hier),N,replace=FALSE))]
    # Compute the SuSt, and add it to the matrix 'SuSt'
    SuSt.Hs[i,c]=mean(basic.stats(subset)$Hs,na.rm=TRUE)
  }
}
rm(c)
save(SuSt.Hs,file="Data/SuSt_GeneDiversityHs.Rda")
load("Data/SuSt_GeneDiversityHs.Rda")

# !!! #  Computing time was very long and time increased exponentially with number of markers,
# !!! #  so we needed to implement parallel programmming to achieve a large number of random bootstraps

# Plotting the mean +- sd of the stat for each number of randomly sampled markers
# Add a grey ribbon for the confidence interval (95%) of the "true" Hst on the global data set, computed with 1,000 bootstraps
df=data.frame(N=samp.n, Mean=apply(SuSt.Hs, 2, mean), sd=apply(SuSt.Hs, 2, sd))
plot.Hs=ggplot(data=df, aes(x=N, y=Mean)) +
  geom_point(colour="Black",size=1)+
  geom_line() +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2) +
  geom_hline(yintercept=c(IC.Hs[1], IC.Hs[2]), linetype="dashed", size=1) +
  xlab("Number of markers randomly sampled") + ylab("Proportion of monomorphic loci (mean ± sd)") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
plot.Hs
ggsave(paste(figuresdir,"/Sensitivity/plot.GeneDiversityHs.png",sep=""),
       device="png",dpi=320,units="cm",width=20,height=16)

###########
# PROPORTION OF REPLICATES INSIDE THE CONFIDENCE INTERVAL OF THE TRUE VALUE OF THE SU. STATISTIC
# TYPE I ERROR: rejection of a true null hypothesis (false positive), i.e. seing a difference between the original data set and a subset

# For each number of marker, compute the proportion of values inside the CI of the true value
# A matrix of 2 columns with 1/ number of markers and 2/ proportion of exact values
true.Hs=matrix(NA,ncol=2,nrow=length(samp.n))
true.Hs[,1]=samp.n
true.Hs[,2]=apply(SuSt.Hs, 2, function(x) (sum(x>IC.Hs[1] & x<IC.Hs[2])/length(x)))
colnames(true.Hs)=c("N.Markers","Prop.Exact.Values")
true.Hs
# Plotting the proportion of values inside the CI of the true value
# Add a horizontal line at 95%
plot.prop.Hs=ggplot(data=as.data.frame(true.Hs), aes(x=N.Markers, y=Prop.Exact.Values)) +
  geom_point(colour="Black",size=1)+
  geom_line() +
  geom_hline(yintercept=0.95, linetype="dashed", size=1) +
  xlab("Number of markers randomly sampled") + ylab("Proportion of values inside the CI of the true value\nfor the proportion of monomorphic loci") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
plot.prop.Hs
ggsave(paste(figuresdir,"/Sensitivity/plot.trueValues.GeneDiversityHs.png",sep=""),
       device="png",dpi=320,units="cm",width=20,height=16)



#==========================================================
# SUBSETTING ONLY THE MOST INFORMATIVE MARKERS (SUBSETTING AFTER SORTING BY PROPORTION OF MISSING DATA)
#==========================================================


#=====================================================================================================
# STUDY OF A SINGLE SU. ST.
#-----------------------------------------------------------------------------------------------------

# MAKE AN INDEX OF LOCI SORTED PER PROPORTION OF MISSING DATA
missing.locus=c() # Proportion of missing data per locus
for (i in 2:ncol(Pparva.hier)) {
  missing.locus[i]=sum(is.na(Pparva.hier[,i]))/length(Pparva.hier[,i])
}
# Sort loci per increasing proportion of missing data and create an ordered index
(order.idx=order(missing.locus))
missing.locus[order(missing.locus)]
# Hence, by sampling the N first loci of this index, we'll be assured to sample only the N most informative markers (in terms of missing data)



#-----------------------------------------------------------------------------------------------------
# MEAN GENETIC DIVERSITY ACROSS ALL LOCI: Hst
#-----------------------------------------------------------------------------------------------------
# We chose to investigate first: mean genetic diversity across all loci, the Hst stat. of hierfstat's function basic.stats() (Nei, 1987)

# Confidence interval (95%) of the "true" Hst on the global data set was computed with 1,000 bootstraps
# BOOTSTRAP MEAN Hst
boots=numeric(1000)
for (i in 1:1000){
  print(i)
  boots[i]=basic.stats(Pparva.hier[,c(1,sample(2:(ncol(Pparva.hier)),replace=T))])$overall[3]
}
save(boots,file="Data/bootstrap.mean.Hst.ordered")
load(file="Data/bootstrap.mean.Hst.ordered")
hist(boots)
(IC.Hst.ordered=quantile(boots,c(.025,.975)))
(mean.Hst.ordered=mean(boots))
(sd.Hst.ordered=sd(boots))

# We suspected a bias in the global data set that could be imputed to missing data (see subsequent figures)
# So we also computed the:
# Confidence interval (95%) of the "true" Hst on the reduced data set (N=2112) with 1,000 bootstraps
# BOOTSTRAP MEAN Hst
boots=numeric(1000)
for (i in 1:1000){
  print(i)
  # Random resampling within the 2112 best markers
  boots[i]=basic.stats(Pparva.hier[,c(1,sample(order.idx[1:2112],replace=T))])$overall[3] # random sampling within the 2112 best loci, those taht have been sampled
}
save(boots,file="Data/bootstrap.mean.Hst.reduced")
load(file="Data/bootstrap.mean.Hst.reduced")
hist(boots)
(IC.Hst.reduced=quantile(boots,c(.025,.975)))
(mean.Hst.reduced=mean(boots))
(sd.Hst.reduced=sd(boots))


# Build a vector for each value of N (nb of randomly sampled SNPs), with iterations in columns and SuSt in rows
samp.n=c(seq(100,3900,100),3999) # Number of markers to sample
SuSt.Hst.ordered=rep(NA,length(samp.n)) # Index of each row correspond to one randomised replicate of the SuSt
SuSt.Hst.upper=rep(NA,length(samp.n))
SuSt.Hst.lower=rep(NA,length(samp.n))
# Index of the column give the number of markers sampled: (N -101)/100 markers sampled (-1 for first column of pop & -100 for setting lowest number of markers to 100)
# Indices are a sequence from 100 to N markers by step of 100 markers


# From 1 sampled SNP only to all SNPs sampled
c=0
for (N in samp.n) {
  c=c+1
  # The SuSt is computed from a random subset of N markers
  cat("\n", N,"markers sampled\n")
  # Subset N most informative markers at each iteration (subset N columns in genind@tab given the index 'order.idx')
  subset=Pparva.hier[,c(1,order.idx[1:N])]
  # Compute the SuSt, and add it to the matrix 'SuSt'
  SuSt.Hst.ordered[c]=basic.stats(subset)$overall[3]
  # Compute the CI of the SuSt by resampling loci (variance of loci chosen)
  boots=numeric(1000)
  for (i in 1:1000){
    cat("Bootstrap #", i, "... ")
    # Random resampling of loci within the N best markers
    boots[i]=basic.stats(subset[,c(1,sample(2:(ncol(subset)),replace=T))])$overall[3] # random sampling within the N best loci, those that have been sampled
  }
  SuSt.Hst.upper[c] = quantile(boots,c(.975))
  SuSt.Hst.lower[c] = quantile(boots,c(.025))
}
SuSt.Hst.ordered = data.frame(SuSt.Hst.ordered, SuSt.Hst.upper, SuSt.Hst.lower)
save(SuSt.Hst.ordered,file="Data/SuSt_Hst_ordered.boot.Rda")
load("Data/SuSt_Hst_ordered.boot.Rda")

# Find the value where missing data overcome 0.55, correspond to the inflection point of the SuSt around locus n°3500
length(missing.locus)-length(order.idx[which(missing.locus[order(missing.locus)]>0.55)]) # Locus 3511
idx=length(missing.locus)-length(order.idx[which(missing.locus[order(missing.locus)]>0.55)])

# Plotting the mean of the stat for each number of sampled markers by order
# Add a grey ribbon for the confidence interval (95%) of the "true" Hst on the global data set, computed with 1,000 bootstraps
df=data.frame(N=samp.n, Mean=SuSt.Hst.ordered$SuSt.Hst.ordered, lower = SuSt.Hst.lower, upper = SuSt.Hst.upper)
plot.Hst.ordered=ggplot(data=df, aes(x=N, y=Mean)) +
  geom_point(colour="Black",size=1)+
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
  geom_hline(yintercept=c(IC.Hst.ordered[1], IC.Hst.ordered[2]), linetype="dashed", size=1) +
  geom_hline(yintercept=c(IC.Hst.reduced[1], IC.Hst.reduced[2]), linetype="dashed", col="DarkGrey", size=1) +
  geom_vline(xintercept=idx,col="Red") +
  xlab("Number of selected markers") + ylab("Mean genetic diversity across all loci ") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
plot.Hst.ordered
# Error bars are the CI by resampling loci (variance of loci sampling)
ggsave(paste(figuresdir,"/Sensitivity/plot.Hst.ordered.png",sep=""),
       device="png",dpi=320,units="cm",width=20,height=16)




#-----------------------------------------------------------------------------------------------------
# PROPORTION OF LOCI WITH NULL GENE DIVERSITY (I.E. MONOMORPHIC LOCI) --> Hs
#-----------------------------------------------------------------------------------------------------
# Then we chose to investigate proportion of loci with null gene diversity (Hs)
# We computed the mean proportion


# Confidence interval (95%) of the "true" Hs on the global data set was computed with 1,000 bootstraps
# BOOTSTRAP MEAN Hs
boots=numeric(1000)
for (i in 1:1000){
  print(i)
  boots[i]=mean(basic.stats(Pparva.hier[,c(1,sample(2:(ncol(Pparva.hier)),replace=T))])$Hs,na.rm=TRUE)
}
save(boots,file="Data/bootstrap.mean.Hs.ordered")
load(file="Data/bootstrap.mean.Hs.ordered")
hist(boots)
(IC.Hs.ordered=quantile(boots,c(.025,.975)))
(mean.Hs.ordered=mean(boots))
(sd.Hs.ordered=sd(boots))

# We suspected a bias in the global data set that could be imputed to missing data (see subsequent figures)
# So we also computed the:
# Confidence interval (95%) of the "true" Hst on the reduced data set (N=2112) with 1,000 bootstraps
# BOOTSTRAP MEAN Hst

# We suspected a bias in the global data set that could be imputed to missing data (see subsequent figures)
# So we also computed the:
# Confidence interval (95%) of the "true" Hst on the reduced data set (N=2112) with 1,000 bootstraps
# BOOTSTRAP MEAN Hst
boots=numeric(1000)
for (i in 1:1000){
  print(i)
  # Random resampling within the 2112 best markers
  boots[i]=mean(basic.stats(Pparva.hier[,c(1,sample(order.idx[1:2112],replace=T))])$Hs,na.rm=TRUE) # random sampling within the 2112 best loci, those that have been sampled for analysis
}
boots[i]=mean(basic.stats(Pparva.hier[,c(1,sample(order.idx[1:2112],replace=T))])$Hs,na.rm=TRUE) # random sampling within the 2112 best loci, those that have been sampled for analysis

save(boots,file="Data/bootstrap.mean.Hs.reduced")
load(file="Data/bootstrap.mean.Hs.reduced")
hist(boots)
(IC.Hs.reduced=quantile(boots,c(.025,.975)))
(mean.Hs.reduced=mean(boots))
(sd.Hs.reduced=sd(boots))


# Build a vector for each value of N (nb of randomly sampled SNPs), with iterations in columns and SuSt in rows
samp.n=c(seq(100,3900,100),3999) # Number of markers to sample
SuSt.Hst.ordered=rep(NA,length(samp.n)) # Index of each row correspond to one randomised replicate of the SuSt
# Index of the column give the number of markers sampled: (N -101)/100 markers sampled (-1 for first column of pop & -100 for setting lowest number of markers to 100)
# Indices are a sequence from 100 to N markers by step of 100 markers


# From 1 sampled SNP only to all SNPs sampled
c=0
for (N in samp.n) {
  c=c+1
  # The SuSt is computed from a random subset of N markers
  cat(N,"markers sampled\n")
  # Subset N most informative markers at each iteration (subset N columns in genind@tab given the index 'order.idx')
  subset=Pparva.hier[,c(1,order.idx[1:N])]
  # Compute the SuSt, and add it to the matrix 'SuSt'
  SuSt.Hs.ordered[c]=mean(basic.stats(subset)$Hs,na.rm=TRUE)
}
save(SuSt.Hs.ordered,file="Data/SuSt.Hs.ordered.Rda")
load("Data/SuSt.Hs.ordered.Rda")

# Find the value where missing data overcome 0.55, correspond to the inflection point of the SuSt around locus n°3500
length(missing.locus)-length(order.idx[which(missing.locus[order(missing.locus)]>0.55)])
idx=length(missing.locus)-length(order.idx[which(missing.locus[order(missing.locus)]>0.55)])

# Plotting the mean of the stat for each number of sampled markers by order
# Add a grey ribbon for the confidence interval (95%) of the "true" Hst on the global data set, computed with 1,000 bootstraps
df=data.frame(N=samp.n, Mean=SuSt.Hs.ordered)
plot.Hs.ordered=ggplot(data=df, aes(x=N, y=Mean)) +
  geom_point(colour="Black",size=1)+
  geom_line() +
  geom_hline(yintercept=c(IC.Hs.ordered[1], IC.Hs.ordered[2]), linetype="dashed", size=1) +
  geom_hline(yintercept=c(IC.Hs.reduced[1], IC.Hs.reduced[2]), linetype="dashed", col="DarkGrey", size=1) +
  geom_vline(xintercept=idx,col="Red") +
  xlab("Number of selected markers") + ylab("Proportion of monomorphic loci") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
        axis.title.x = element_text(color="black", size=14),
        axis.title.y = element_text(color="black", size=14),
        axis.text=element_text(size=14, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=14),
        legend.title=element_text(size=14))
plot.Hs.ordered
ggsave(paste(figuresdir,"/Sensitivity/plot.GeneDiversityHs.ordered.png",sep=""),
       device="png",dpi=320,units="cm",width=20,height=16)




###################
# Randomly selected and sorted by increasing proportion of missing data on the same plot for comparison


###########################################################
#
#       Sensitivity analysis for SuSt Euclidian distances
#
###########################################################

# Our previous sensitivity studies did not care for population structure and assessed only independent SuSt at a global scale (all populations mixed together in a single mean value)
# Hence we now consider population structure by computing SuSt per population

# Because this rapidly produced a lot of values hard to interpret graphically, populations tructure had to be resumed
# As the sensitivity analysis is designed for ABC, we used the same resuming approach, i.e. the Euclidian distance between SuSt

#==========================================================
# RANDOM SET OF MARKERS
#==========================================================





#==========================================================
# SUBSETTING MOST INFORMATIVE MARKERS (LESS OF MISSING DATA)
#==========================================================


#==========================================================
# Save the environment 
#==========================================================
save.image(file="Sensitivity.RData")

#==========================================================
# THE END