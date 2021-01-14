############################################################################
#   INRA - ANR project Pseudorasbora
#
#       SUMMARY
#
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com
# MSc internship

# Supervisor: Scott McCairns
#             INRA UMR ESE


# DESCRIPTION
# This script performs basic descriptive analySes and computes populations summary statistics
# But all "main" analyses are performed in separate scripts with the name of the analysis as file name

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
library(data.table)
# GENETICS
library(ade4)
library(adegenet)
library(hierfstat)
library(ape)
library(pegas)
library(genepop)
library(diveRsity)
library(poppr) # for function 'private_alleles'
# devtools::install_github("thierrygosselin/radiator")
library(radiator) # Genomic converter among different genomic data formats
#devtools::install_github("rystanley/genepopedit",dependencies=TRUE)
library(genepopedit)
#devtools::install_github('wrengels/HWxtest', subdir='pkg')
library(HWxtest)

# CARTOGRAPHY
library(rgdal)
library(raster)
library(gdistance)
library(maps)
library(mapdata)
library(mapplots)
library(rnaturalearth) # For finer cartography of river network
# Install gstudio, for spatial genetic analysis
# install.packages( c("RgoogleMaps","geosphere","proto","sampling","seqinr","spacetime","spdep"), dependencies=TRUE )
# devtools::install_github("dyerlab/gstudio")
# library(gstudio)

# GRAPHIC LIBRARIES
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(ggplotify)
# STATS
library(MCMCglmm)
library(coda)
library(abc)

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


# Import sources of functions
# Only functions developped by myself will be printed in appendix,
# some other functions were offered to scientific community by their owner (e.g. Github)
source("Sources/map_nomenclature.R") # Function to draw scalebar and north arrow
source("Sources/structureLauncher.R") # A function to generate a bash or a qsub files for computation servers


#---------------------------------
# Import data and objects from Scott McCairns, INRA ESE Epix
#---------------------------------
load(paste(datadir,"/Pparva.Rda",sep=""))
load(paste(datadir,"/Pparva.3000.Rda",sep=""))

#---------------------------------
# Import coordinates of sampled locations
Pparva.sites=read.table(paste(datadir,"/Geographic/Pparva.sites.csv",sep=""),sep=";",h=T)

# write.table(genind2df(Pparva),paste(datadir,"/Pparva.alleles.txt",sep=""),sep=" ",quote=FALSE,row.names=F,col.names=F)
labelPops=read.table("Data/STRUCTURE/labelPops.txt",header = TRUE)


# Lists of sampled sites names with their associated coordinates
native.names=c("Jap", "S1", "S2", "S3", "S4", "S6","S9",
               "S10", "S11", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "Tib")
native.coords=Pparva.sites[which(Pparva.sites$Pop %in% native.names),]
invasive.names=c("Aus", "Bel", "Bul1", "Bul2", "Hun", "Ira", "Ita", "Pol", "Spa", "Tur", "UK")
invasive.coords=Pparva.sites[which(Pparva.sites$Pop %in% invasive.names),]

#----------------------------------
# GENIND FOR NATIVE & INVASIVE
# Create native and invasise genind
Native=Pparva[Pparva$pop %in% native.names]
Invasive=Pparva[Pparva$pop %in% invasive.names]

#==========================================================
# EXPLORE & DESCRIBE DATA SET
#==========================================================
# Populations:
# 18 populations in Chinese native area
# 3 additional sites in Tibet, Japan and Taiwan
# 13 populations in European invasive area
# 34 populations in data set
length(levels(Pparva@pop))

# final.df is a data frame object that contains 3999 SNPs loci for 741 indivuals
# all.df is a data frame object that contains 3999 SNPs loci for 849 indivuals
# Pparva is a genind object that contains genotypes, with individuals and populations names

#==========================================================
# SAMPLED SITES
#----------------------------------------------------------

# Lists of sampled sites names with their associated coordinates
native.names=c("Jap", "S1", "S2", "S3", "S4", "S6","S9",
               "S10", "S11", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "Tib")
native.coords=Pparva.sites[which(Pparva.sites$Pop %in% native.names),]
invasive.names=c("Aus", "Bel", "Bul1", "Bul2", "Hun", "Ira", "Ita", "Pol", "Spa", "Tur", "UK")
invasive.coords=Pparva.sites[which(Pparva.sites$Pop %in% invasive.names),]

#----------------------------------
# GENIND FOR NATIVE & INVASIVE
# Create native and invasise genind
Native=Pparva[Pparva$pop %in% native.names]
Invasive=Pparva[Pparva$pop %in% invasive.names]

# populations size
(pop.size=table(Pparva@pop))

#-----------------------------------
# INDIVIDUALS

# Description of data set of individuals
head(Pparva)
summary(Pparva)

# We need another data format, also from adegenet
# Instead of individuals, the new unit is populations
Pparva.pop=genind2genpop(Pparva)
# To which we add population coordinates
Pparva.pop@other$xy=Pparva.sites[,c(2,3)]

# Description of data set of populations
summary(Pparva.pop)

#------------------------------
# MISSING DATA
#------------------------------
load(paste(datadir,"/Pparva.Rda",sep=""))
Pparva.hier=genind2hierfstat(Pparva)

# Missing data per indiviudal
## First, compute a vector of missing data % per individual
missing.ind=c()
for (i in 1:nrow(Pparva.hier)) {
  missing.ind[i]=sum(is.na(Pparva.hier[i,]))/length(Pparva.hier[i,])
}
sum(missing.ind<0.1) # 112 indiv (15%) with less than 10% of missing data
# Arbitrary threshold of 45% that is the percentage of missing data on the global data set
sum(missing.ind<0.45) # 412 indiv (56%) with less than 45% of missing data (percentage of missing data in the global data set)

# Missing data per population
missing.pop=data.frame(unique(Pparva@pop),rep(NA,length(unique(Pparva@pop))))
for (i in unique(Pparva@pop)) {
  missing.pop[which(missing.pop[,1]==i),2]=sum(is.na(genind2df(Pparva[Pparva@pop==i,])))/(sum(Pparva@pop==i)*3999)
}
colnames(missing.pop)=c("Population","Missing data")
missing.pop
# Choose an arbitrary threshold of 0.7 for missing data in a population
missing.pop[missing.pop[,2]>0.7,]
# Some populations have very high amounts of missing data, over 70%: Morocco (89%), Slovenia (85%), S7 (77%), S5 (73%), Tai (76%)
# Zoom on details of the Moroccan population
Mor.data=genind2df(Pparva[Pparva@pop=="Mor",])

# Missing data per allele
## First, compute a vector of missing data % per allele
missing.all=c()
for (i in 1:ncol(Pparva@tab)) {
  missing.all[i]=sum(is.na(Pparva@tab[,i]))/length(Pparva@tab[,i])
}
sum(missing.all<0.1) # There is no allele with less than 10% of missing data
sum(missing.all<0.45) # 4682 of 9146 alleles (52%) have less than 45% of missing data

# Missing data per locus
## First, compute a vector of missing data % per locus
missing.locus=c()
for (i in 1:ncol(Pparva.hier)) {
  missing.locus[i]=sum(is.na(Pparva.hier[,i]))/length(Pparva.hier[,i])
}
sum(missing.locus<0.1) # There is only 1 locus with less than 10% of missing data
sum(missing.locus<0.45) # 2112 of 3999 locus (53%) have less than 45% of missing data
sum(missing.locus>0.5) # 

#------------------------------
# STRATEGY FOR MISSING DATA
#------------------------------
# Do we filter individuals or loci?

png("Figures/missing.ind.png",res=150,width=800,height=600)
hist(missing.ind,breaks=20)
# Inflection point around 0.5-0.6 of missing data per individual
# We removed individuals with more than 60% of missing data
sum(missing.ind<0.6)
sum(missing.ind<0.45)
dev.off()

png("Figures/missing.pop.png",res=150,width=800,height=600)
missing.pop
hist(missing.pop[,2],breaks=20)
abline(v=0.7,lwd=2,col="red")
# % of missing data in pop. dropped after 0.6, and then there was a peak after 0.7
# Cutting pop. above 0.7 of missing data removed these outliers
dev.off()


hist(missing.all,breaks=20)
abline(v=0.45,lwd=2,col="red")

png("Figures/missing.locus.png",res=150,width=800,height=600)
hist(missing.locus,breaks=20)
abline(v=0.45,lwd=2,col="red")
dev.off()
# By retaining only loci with less than 45% of missing data, we kept only half of original markers
# that was the most informative markers.
# We assumed that markers with more than 45% of missing data do not carried more information than markers with less than 45% of missing,
# but this threshold is arbitrary


# N=468 (64%) if we removed all individuals with more than 60% of missing data
# & populations with more than 70% of missing data (5 pops) leading to too few individuals per pop. to retain
(Pparva.clean=Pparva[!(Pparva@pop %in% c("Mor","Tai","Slo","S5","S7")) & missing.ind<0.6])
table(Pparva.clean@pop)
#-------------------------------------
# Keep only SNPs with less than 45% of missing data
# Recompute missing data per locus on the new data set
missing.locus2=c()
for (i in 1:ncol(Pparva.hier)) {
  missing.locus2[i]=sum(is.na(Pparva.hier[,i]))/length(Pparva.hier[,i])
}
sum(missing.locus2<0.1) # There is only 1 locus with less than 10% of missing data
sum(missing.locus2<0.45) # 2112 of 3999 locus (53%) have less than 45% of missing data

(Pparva.clean=Pparva.clean[loc=(missing.locus2<0.45)])
# Percentage of missing data is reduced to 22% instead of 45% by eliminating the half of markers that were the less informative

#-----------------------------------
# Finally, we needed to retain more loci, to gain genetic information for DIY ABC
# We chose to retain the 3,000 best loci instead of 2,112 previously
missing.locus3=c()
for (i in 1:ncol(Pparva.hier)) {
  missing.locus3[i]=sum(is.na(Pparva.hier[,i]))/length(Pparva.hier[,i])
}
sum(missing.locus3<0.1) # There is only 1 locus with less than 10% of missing data
sum(missing.locus3<0.45) # 2112 of 3999 locus (53%) have less than 45% of missing data

(Pparva.3000=Pparva.clean[loc=order(missing.locus3,decreasing=FALSE)[1:3000]])

# From now, Pparva.clean became the new Pparva
# Pparva=Pparva.clean
# save(Pparva,file=file(paste(datadir,"/Pparva.clean.Rda",sep="")))
# After running STRUCTURE on Pparva.clean (2,112 loci), we reconsidered the number of loci to retain, in order to improve genetic information for DIY ABC
# Hence, we took 3,000 loci
# Pparva=Pparva.3000
# save(Pparva,file=file(paste(datadir,"/Pparva.3000.Rda",sep="")))
load(paste(datadir,"/Pparva.3000.Rda",sep=""))
# Recompute genpop
Pparva.pop=genind2genpop(Pparva)


#==========================================================
# Sensitivity analysis to missing data
#==========================================================
#-----------------------------------------------------------------------------------------------------
# MEAN GENETIC DIVERSITY ACROSS ALL LOCI: Hst
#-----------------------------------------------------------------------------------------------------
# Manipulation of a hierfstat object is easier than a genind
load(paste(datadir,"/Pparva.Rda",sep=""))
Pparva.hier=genind2hierfstat(Pparva)
#-----------------------------------------------------------------------
# Markers randomly sampled
#-----------------------------------------------------------------------
# We chose to investigate first: mean genetic diversity across all loci, the Hst stat. of hierfstat's function basic.stats() (Nei, 1987)
load(file="Data/bootstrap.mean.Hst")
limBoot=quantile(boots,c(.025,.975))
(IC.Hst=limBoot)
(mean.Hst=mean(boots))
(sd.Hst=sd(boots))
load("Data/SuSt_Hst.Rda")
samp.n=seq(100,(ncol(Pparva.hier)-101),100) # Number of markers to sample

df=data.frame(N=samp.n, Mean=apply(SuSt.Hst, 2, mean), sd=apply(SuSt.Hst, 2,function(x) quantile(x,0.975))-apply(SuSt.Hst, 2, mean)) # Plot the CI instead of sd of the resampled means
plot.Hst=ggplot(data=df, aes(x=N, y=Mean)) +
  geom_point(colour="Black",size=1)+
  geom_line() +
  geom_errorbar(aes(ymin=Mean-sd, ymax=Mean+sd), width=.2) +
  geom_hline(yintercept=c(IC.Hst[1], IC.Hst[2]), linetype="dashed", size=1) +
  xlab("Number of markers randomly sampled") + ylab("Mean genetic diversity across all loci (± C.I.)") +
  xlim(0, 4000) +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=18, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=18,hjust = 0.5),
        axis.title.x = element_text(color="black", size=18),
        axis.title.y = element_text(color="black", size=18),
        axis.text=element_text(size=18, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18))
plot.Hst

#-----------------------------------------------------------------------
# Markers ordered by proportion of missing data
#-----------------------------------------------------------------------
# MAKE AN INDEX OF LOCI SORTED PER PROPORTION OF MISSING DATA
missing.locus=c() # Proportion of missing data per locus
for (i in 2:ncol(Pparva.hier)) {
  missing.locus[i]=sum(is.na(Pparva.hier[,i]))/length(Pparva.hier[,i])
}
# Sort loci per increasing proportion of missing data and create an ordered index
(order.idx=order(missing.locus))
missing.locus[order(missing.locus)]

load(file="Data/bootstrap.mean.Hst.ordered")
(IC.Hst.ordered=quantile(boots,c(.025,.975)))
(mean.Hst.ordered=mean(boots))
(sd.Hst.ordered=sd(boots))
load(file="Data/bootstrap.mean.Hst.reduced")
(IC.Hst.reduced=quantile(boots,c(.025,.975)))
(mean.Hst.reduced=mean(boots))
(sd.Hst.reduced=sd(boots))
load("Data/SuSt_Hst_ordered.boot.Rda")
samp.n=c(seq(100,3900,100),3999) # Number of markers to sample

# Find the value where missing data overcome 0.55, correspond to the inflection point of the SuSt around locus n°3500
length(missing.locus)-length(order.idx[which(missing.locus[order(missing.locus)]>0.55)]) # Locus 3511
idx=length(missing.locus)-length(order.idx[which(missing.locus[order(missing.locus)]>0.55)])

# Plotting the mean of the stat for each number of sampled markers by order
# Add a grey ribbon for the confidence interval (95%) of the "true" Hst on the global data set, computed with 1,000 bootstraps
df=data.frame(N=samp.n, Mean=SuSt.Hst.ordered$SuSt.Hst.ordered, lower = SuSt.Hst.ordered$SuSt.Hst.lower,
              upper = SuSt.Hst.ordered$SuSt.Hst.upper)
plot.Hst.ordered=ggplot(data=df, aes(x=N, y=Mean)) +
  geom_point(colour="Black",size=1)+
  geom_line() +
  geom_errorbar(aes(ymin=lower, ymax=upper), width=.2) +
  geom_hline(yintercept=c(IC.Hst.ordered[1], IC.Hst.ordered[2]), linetype="dashed", size=1) +
  geom_hline(yintercept=c(IC.Hst.reduced[1], IC.Hst.reduced[2]), linetype="dashed", col="DarkGrey", size=1) +
  # geom_vline(xintercept=idx,col="Red", size = 0.8) +
  xlab("Number of selected markers") + ylab("Mean genetic diversity across all loci (± C.I.)") +
  theme(axis.line = element_line(colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=18, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=18,hjust = 0.5),
        axis.title.x = element_text(color="black", size=18),
        axis.title.y = element_text(color="black", size=18),
        axis.text=element_text(size=18, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=18),
        legend.title=element_text(size=18))
plot.Hst.ordered
# Error bars are the CI by resampling loci (variance of loci sampling)

ggarrange(plot.Hst,plot.Hst.ordered,widths=1:1,heights=1:1,labels="auto")
ggsave(paste(figuresdir,"/Sensitivity/PlotArrange.Hst.ordered.png",sep=""),
       device="png",dpi=320,units="cm",width=40,height=16)


load(paste(datadir,"/Pparva.3000.Rda",sep=""))


#==========================================================
# Part 1. Describe genetic diversity in sampled sites
#==========================================================

#==========================================================
# Basic population genetics statistics
#----------------------------------------------------------
# Use the dataset with 3,000 loci

# Compute a table with 6 predictors (columns) on 28 populations (rows):
# N — number of sampled individuals
# AR —allelic richness
# HE — expected heterozygosity
# HO — observed heterozygosity
# Fst
# Fis

# Barker et al. 2017 described its RAD seq data set like that:
# N, number of individuals analysed; P, proportion of variable loci; number of private alleles (Priv.); NAR, mean allelic richness corrected
# for uneven sample size; NPAR, mean private allelic richness corrected for uneven sample size; HO, observed heterozygosity for
# polymorphic loci; p, nucleotide diversity. Standard errors (SE) for NAR, NPAR, HO and p are indicated in parentheses.



# Hohenlohe et al. 2011 described He and Hobs for SNP (p.3):
# At each putative
# locus, we calculated expected heterozygosity as Hexp =
#   1)Pni(ni ) 1) ⁄ n(n ) 1), where ni is the count of allele i in
# the sample and n = Pni, and observed heterozygosity
# Hobs as the proportion of individuals that appear
# heterozygous at the locus. We also calculated FIS = 1 )
# (Hobs ⁄Hexp), which provides a measure of the deviation
# of genotype frequencies (i.e. observed heterozygosity)
# from Hardy–Weinberg proportions.
load(paste(datadir,"/Pparva.3000.Rda",sep=""))
# convert the genind to a hierfstat object
Pparva.hier=genind2hierfstat(Pparva)
# load(paste(datadir,"/Pparva.clean.Rda",sep=""))

basic.genestats=as.data.frame(matrix(nrow = 29, ncol = 14))
colnames(basic.genestats)=c("Population","Range","N","Allelic_Richness","Private","Monomorphic_Loci","HE","HE.sem","HO","HO.sem","Fst","Fst.sem","Fis","Fis.sem")
# Get population names
basic.genestats$Population=unique(Pparva$pop)

#---------------------------------------------------------
# Range — Asia or Europe
basic.genestats$Range=c(rep("Europe",7),"Asia","Europe",rep("Asia",16),"Europe","Asia",rep("Europe",2))

#---------------------------------------------------------
# N — number of sampled individuals
basic.genestats$N=table(Pparva.hier$pop)

#---------------------------------------------------------
# AR - mean Allelic Richness per population, corrected for sample size
# Estimation by rarefaction
# Hierfstat
basic.genestats$Allelic_Richness = round(colMeans(allelic.richness(Pparva.hier)$Ar, na.rm = TRUE), digits = 3)

# for (i in 1:29) {
#   print(i)
#   sub = subset(Pparva.hier, pop==as.character(basic.genestats$Population[i]))
#   basic.genestats$Allelic_Richness[i]=round(mean(allelic.richness(sub)$Ar,na.rm=TRUE),digits=3)
# }


#---------------------------------------------------------
# NUMBER OF PRIVATE ALLELES PER POPULATION
# basic.genestats$Private=apply(private_alleles(Pparva),1,sum)
# basic.genestats$Private = apply(private_alleles(Pparva,form = alleles ~ ., report = "table",
                                              # level = "population"),1,sum)
# private_alleles(Pparva, form = allele ~ ., level = "population")
private = as.data.frame(private_alleles(Pparva, form = allele ~ .,report = "table", level = "population"))
private
# private_alleles count the number of occurences of the private allele, hence it can be different of 0/1
private[private > 0] = 1
basic.genestats$Private = apply(private, 1, sum)
sum(private)

#---------------------------------------------------------
# PROPORTION OF LOCI WITH NULL GENE DIVERSITY (I.E. MONOMORPHIC LOCI) --> Hs
basic.genestats$Monomorphic_Loci=round(apply(basic.stats(Pparva.hier)$Hs,2,function(x) sum(x==0,na.rm=TRUE)/length(x)),digits=3)       



#---------------------------------------------------------
# HE — mean expected heterozygosity across polymorphic loci +- s.e.m.
#basic.genestats$HE=round(Hs(Pparva),digits=3)

# Bootstrap of HE - random resampling of individuals within populations
# for (i in 1:29) {
#   cat("pop",i,"\n")
#   boots=numeric(10)
#   for (j in 1:10){
#     cat("boot",j,"\n")
#     boots[j]=Hs(Pparva[sample(which(Pparva$pop==basic.genestats$Population[i]),replace=TRUE),])
#   }
#   basic.genestats$HE[i]=round(mean(boots),digits=3)
#   basic.genestats$HE.sem[i]=round(sd(boots),digits=3)
# }


# Bootstrap of HE - random resampling of individuals within populations
n.pop=seppop(Pparva)
for (i in 1:29) {
  cat("--------------------------\npop",i,"\n--------------------------\n")
  boots=numeric(100)
  for (j in 1:100){
    cat("pop",i,"boot",j,"\n")
    temp=c()
    subs=n.pop[[i]]
    temp=summary(subs[sample(1:nrow(n.pop[[i]]@tab),replace=TRUE),])$He
    boots[j]=mean(temp[temp>0],na.rm=TRUE)
    # boots[j]=mean(summary(Pparva[which(Pparva$pop==basic.genestats$Population[i]),sample(1:4216,replace=TRUE)])$He)
    rm(temp)
  }
  basic.genestats$HE[i]=round(mean(boots),digits=3)
  basic.genestats$HE.sem[i]=round(sd(boots),digits=3)
  print(round(mean(boots),digits=3))
  print(round(sd(boots),digits=3))
}

barplot(basic.genestats$HE)


#---------------------------------------------------------
# HO — mean observed individual heterozygosity across polymorphic loci +- s.e.m.
n.pop=seppop(Pparva) 
# (basic.genestats$HO=round(do.call("c", lapply(n.pop, function(x) mean(summary(x)$Hobs,na.rm=TRUE))),digits=3)) # Compute mean Hobs for each pop

# Bootstrap of Ho - random resampling individuals within populations
for (i in 1:29) {
  cat("--------------------------\npop",i,"\n--------------------------\n")
  boots=numeric(100)
  for (j in 1:100){
    cat("pop",i,"boot",j,"\n")
    temp=c()
    subs=n.pop[[i]]
    temp=summary(subs[sample(1:nrow(n.pop[[i]]@tab),replace=TRUE),])$Hobs
    boots[j]=mean(temp[temp>0],na.rm=TRUE)
    # boots[j]=mean(summary(Pparva[which(Pparva$pop==basic.genestats$Population[i]),sample(1:4216,replace=TRUE)])$He)
    rm(temp)
  }
  basic.genestats$HO[i]=round(mean(boots),digits=3)
  basic.genestats$HO.sem[i]=round(sd(boots),digits=3)
}

barplot(basic.genestats$HO)


#---------------------------------------------------------
# Weir and Cockerham Fst per population  +- s.e.m.
# Compute Weir and Cockerham Fstat
# As the mean Fst per population in:
# - intercontinental range, to compare native and invasive populations in a single metapopulation
# - within continental range (Asia OR Europe) for estimating differentiation with geographical neighbours
(Pparva.pairwiseWC=pairwise.WCfst(Pparva.hier))
basic.genestats$Fst=round(apply(Pparva.pairwiseWC,2, function(x) mean(x, na.rm=TRUE)),digits=3) # Mean Fst
# Bootstrap of Fst s.e.m. with 1,000 bootstraps, by randomly resampling Fst from pairwise Fst list for a given population
sem=c()
for (j in 1:29) {
  print(j)
  boots=numeric(1000)
  for (i in 1:1000){
    boots[i]=mean(sample(Pparva.pairwiseWC,replace=TRUE),na.rm=TRUE)
  }
  sem[j]=sd(boots)# Fst s.e.m
}
basic.genestats$Fst.sem=round(sem,digits=3)

# # Asia
# Pparva.hier.native=genind2hierfstat(native)
# (Pparva.pairwiseWC.native=pairwise.WCfst(Pparva.hier.native))
# 
# # Europe
# Pparva.hier.invasive=genind2hierfstat(invasive)
# (Pparva.pairwiseWC.invasive=pairwise.WCfst(Pparva.hier.invasive))

#----------------------------------------------------------
# Fis per population +- s.e.m.
# with 1000 bootstraps
boot=boot.ppfis(dat=Pparva.hier,nboot=1000,quant=c(0.025,0.975),diploid=TRUE)
basic.genestats$Fis=round(apply(boot$fis.ci,1, mean),digits=3) # Mean
basic.genestats$Fis.sem=round((as.numeric(unlist(boot$fis.ci[2]))-as.numeric(unlist(boot$fis.ci[1])))/3.92,digits=3) # Standard error


########
# Sort populations by range (first Asia, the Europe)
basic.genestats=basic.genestats[c(8,10,20,22:25,11:19,21,27,1:7,9,26,28,29),]

# And then export table:
save(basic.genestats,file="Data/basic.genestats.Rda")
write.csv(basic.genestats,file="Tables/basic.genestats.csv",row.names = FALSE)
load(basic.genestats,file="Data/basic.genestats.Rda")

#----------------------------------------------------------
# Format in publication mode
basic.genestats.formatted = basic.genestats[,1:6]
# Paste estimates with s.e.m. in a readable format
basic.genestats.formatted$He = paste(basic.genestats$HE, " (±", basic.genestats$HE.sem, ")", sep = "" )
basic.genestats.formatted$Ho = paste(basic.genestats$HO, " (±", basic.genestats$HO.sem, ")", sep = "" )
basic.genestats.formatted$Fst = paste(basic.genestats$Fst, " (±", basic.genestats$Fst.sem, ")", sep = "" )
basic.genestats.formatted$Fis = paste(basic.genestats$Fis, " (±", basic.genestats$Fis.sem, ")", sep = "" )
write.csv(basic.genestats.formatted,file = "Tables/basic.genestats.formatted.csv", row.names = FALSE)




#--------------------------------------------------------------
# !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
# HIERARCHICAL FST
# Fst among sampled sites, among demes (DAPC and STRUCTURE results), among continents and isolated populations (Jap, Ira)

(Pparva.wc=wc(Pparva.hier)) # Mean Fst and Fis


fstat(Pparva)

Fst(as.loci(Pparva))
(Gtest=gstat.randtest(Pparva,nsim=999))
plot(Gtest)

# Compute Weir & Cockerham Fst global pairwise
(Pparva.pairwiseWC=pairwise.WCfst(Pparva.hier))

# Compute Weir & Cockerham Fst native pairwise
(Pparva.pairwiseWC.native=pairwise.WCfst(Pparva.hier.native))

# Compute Weir & Cockerham Fst invasive pairwise
(Pparva.pairwiseWC.invasive=pairwise.WCfst(Pparva.hier.invasive))


#===========================================================
# DEME ANALYSIS
#===========================================================
# From now, we switched from sampling level to deme level
# Deme defined according to STRUCTURE results

#===========================================================
# Pairwise Fst among demes: 
# Pairwise Fst per demes (Western Europe, Easter Europe, Iran, Japan, North China, South China, North Central China, South Central China, Continental China, Tibet)
# Weir & Cockerham pairwise Fst, significant ones in bold

# Putative demes were:
#   in ASIA
# - North China: S13 S14 S15 S17
# - North Central China: S1 S2 S3 S16
# - Admixed coastal zone: S3 S4 S6
# - Admixed continental zone: S10 S11
# - South China: S9 S18 S19 S20
# - Japan
# - Tibet
#   in EUROPE
# - Western Europe: AUS BEL HUN ITA POL SPA UK
# - Eastern Europe: BUL1 BUL2 TUR
# - Iran
load(paste(datadir,"/Pparva.3000.Rda",sep=""))

Pparva.hier=genind2hierfstat(Pparva)
Pparva.demes=Pparva.hier
# Remove the undesirable populations (outlier, admixed zone)
# Pparva.demes=Pparva.demes[-which(Pparva.demes$pop %in% c("Mor","Slo","Tai","S5","S7")),]
# Assign populations to demes in the genind object (genind object of 3,000 loci)
Pparva.demes$pop=as.character(Pparva.demes$pop)
Pparva.demes$pop[Pparva.demes$pop %in% c("Aus","Bel","Hun","Ita","Spa","Pol","UK")]="WesternEurope"
Pparva.demes$pop[Pparva.demes$pop %in% c("Bul1","Bul2","Tur")]="EasternEurope"
Pparva.demes$pop[Pparva.demes$pop %in% c("S1","S2","S3","S16")]="NorthCentralChina"
Pparva.demes$pop[Pparva.demes$pop %in% c("S13","S14","S15","S17")]="NorthChina"
Pparva.demes$pop[Pparva.demes$pop %in% c("S4","S6")]="AdmixedCoastalChina"
Pparva.demes$pop[Pparva.demes$pop %in% c("S9","S18","S19","S20")]="SouthChina"
Pparva.demes$pop[Pparva.demes$pop %in% c("S10","S11")]="AdmixedContinentalChina"

# Get pop names (as factors) in the right order
Pparva.demes$pop=as.factor(Pparva.demes$pop)
Pparva.demes$pop=factor(Pparva.demes$pop,levels=c("WesternEurope","EasternEurope","Ira","Jap","AdmixedContinentalChina","SouthChina",
                                                  "AdmixedCoastalChina","NorthChina","NorthCentralChina","Tib"))
# Due to a bug in boot.ppfst, populations need to be encoded as numeric and the data frame ordered by population number
v.pop=as.data.frame(Pparva.demes$pop) # save a data frame of population names for the significance matrix
v.pop$num=as.numeric(Pparva.demes$pop)
v.pop=unique(v.pop)
write.table(v.pop,file="Tables/Pairwise.Fst.labels.csv") # Save labels of pop names with corresponding numeric label
Pparva.demes$pop=as.numeric(Pparva.demes$pop)
setorder(Pparva.demes,pop)

# Compute pairwise Fst (Weir & Cockerham's Fst) between demes
(Pparva.pairwiseWC=pairwise.WCfst(Pparva.demes))
(Pparva.pairwiseWC=round(Pparva.pairwiseWC,digits=3))
# Keep only the upper triangular, lower triangular will be filled with pairwise Nei's distances
Pparva.pairwiseWC[lower.tri(Pparva.pairwiseWC, diag = TRUE)]=0
# And now compute the Nei's distances
(Pparva.pairwiseNei=pairwise.neifst(Pparva.demes))
(Pparva.pairwiseNei=round(Pparva.pairwiseNei,digits=3))
Pparva.pairwiseNei[upper.tri(Pparva.pairwiseNei, diag = TRUE)]=0
# Combine the two triangular matrices of WC Fst and Nei's Fst
Pparva.pairwiseWC[is.na(Pparva.pairwiseWC)]=0
Pparva.pairwiseNei[is.na(Pparva.pairwiseNei)]=0
(mat=Pparva.pairwiseWC + Pparva.pairwiseNei)
write.table(mat,file="Tables/Pairwise.Fst.csv")
rm(mat)

# Compute significance of the Fst by bootstrap
bootppfst=boot.ppfst(dat=Pparva.demes,nboot=1000,quant=c(0.005,0.995),dig=3) # CI at 1%
# make a matrix with upper limit in top triangular matrix and lower limit in bottom triangular matrix
# Fst were not significantly different from 0 when the lower CI was negative or equal to 0
(mat.ll=t(bootppfst$ll))
(mat.ul=bootppfst$ul)
mat.ll[is.na(mat.ll)]=0
mat.ul[is.na(mat.ul)]=0
(mat=mat.ll + mat.ul)
write.table(mat,file="Tables/Pairwise.Fst.significance.csv")
rm(mat)
# All pairwise Fst were significant (higher than expected just by chance). Putative demes seemed well defined

# All pairwise comparisons are signiﬁcantly differentiated (P < 0.01)


#----------------------------------------------------------
# Neighbour-Joining tree estimated on pairwise distances
# with 'ape' package
?nj() # require a distance matrix
# NJ tree with Weir & Cockerham distances
# and populations in order:
# "WesternEurope","EasternEurope","Ira","Jap","AdmixedContinentalChina","SouthChina",
#   "AdmixedCoastalChina","NorthChina","NorthCentralChina","Tib"
(Pparva.pairwiseWC=pairwise.WCfst(Pparva.demes))
colnames(Pparva.pairwiseWC)=c("WesternEurope","EasternEurope","Ira","Jap","AdmixedContinentalChina","SouthChina",
                              "AdmixedCoastalChina","NorthChina","NorthCentralChina","Tib")
rownames(Pparva.pairwiseWC)=c("WesternEurope","EasternEurope","Ira","Jap","AdmixedContinentalChina","SouthChina",
                              "AdmixedCoastalChina","NorthChina","NorthCentralChina","Tib")
nj.fst=nj(as.dist(Pparva.pairwiseWC))
plot(nj.fst, main="NJ tree for Weir & Cockerham's distances", type="phylogram")
plot(nj.fst, main="NJ tree for Weir & Cockerham's distances", type="unrooted")

# NJ tree with Nei's distances
(Pparva.pairwiseNei=pairwise.neifst(Pparva.demes))
colnames(Pparva.pairwiseNei)=c("WesternEurope","EasternEurope","Ira","Jap","AdmixedContinentalChina","SouthChina",
                               "AdmixedCoastalChina","NorthChina","NorthCentralChina","Tib")
rownames(Pparva.pairwiseNei)=c("WesternEurope","EasternEurope","Ira","Jap","AdmixedContinentalChina","SouthChina",
                               "AdmixedCoastalChina","NorthChina","NorthCentralChina","Tib")
nj.Nei=nj(as.dist(Pparva.pairwiseNei))
plot(nj.Nei, main="NJ tree for Nei's distances", type="phylogram")
plot(nj.Nei, main="NJ tree for Nei's distances", type="unrooted")

# Compute a bootstrapped tree of pairwise genetic distances





# Heterozygosities weighted by group sizes
(Pparva.pairwiseNei=pairwise.fst(Pparva,res.type="matrix"))
save(Pparva.pairwiseNei,file=file(paste(datadir,"/Pparva.pairwiseNei.Rda",sep="")))
load(paste(datadir,"/Pparva.pairwiseNei.Rda",sep=""))
# Summary statistics of pairwise Nei's distance
mean(Pparva.pairwiseNei) # Mean
sd(Pparva.pairwiseNei)/sqrt(nrow(Pparva.pairwiseNei)) # standard error of the mean
max(Pparva.pairwiseNei)
min(Pparva.pairwiseNei)
# Mean pairwise Fst (Nei's Fst) was of  (S.E.M = ), in range -


pairNei.tree=nj(Pparva.pairwiseNei)
# Draw a tree plot of genetic distances (Fst) between populations
plot(pairNei.tree, type="unr", tip.col=funky(17)[-17], font=2)
annot <- round(pairNei.tree$edge.length,2)
edgelabels(annot[annot>0], which(annot>0), frame="n")
add.scale.bar()
# Draw a boxplot of Fst per population
temp=Pparva.pairwiseNei
diag(temp)=NA
boxplot(temp, col=funky(17)[-17], xlab="Population", ylab="Fst")


# basic statistics (hierfstat) per loci and overall
(Pparva.stat=basic.stats(Pparva.hier))

# other basic statistics (allelic richness) not included in the basic.stats function (hierfstat)
# Allelic richness per loci: number of alleles corrected for sample size
#(Pparva.allrich=(as.data.frame(allelic.richness(Pparva.hier))))




#==========================================================
# Basic statistics per demes
#==========================================================
# Individuals pooled in defined demes
load(paste(datadir,"/Pparva.3000.Rda",sep=""))
Pparva.hier=genind2hierfstat(Pparva)
Pparva.demes=Pparva.hier
# Remove the undesirable populations (outlier, admixed zone)
# Pparva.demes=Pparva.demes[-which(Pparva.demes$pop %in% c("Mor","Slo","Tai","S5","S7")),]
# Assign populations to demes in the genind object (genind object of 3,000 loci)
Pparva.demes$pop=as.character(Pparva.demes$pop)
Pparva.demes$pop[Pparva.demes$pop %in% c("Aus","Bel","Hun","Ita","Spa","Pol","UK")]="WesternEurope"
Pparva.demes$pop[Pparva.demes$pop %in% c("Bul1","Bul2","Tur")]="EasternEurope"
Pparva.demes$pop[Pparva.demes$pop %in% c("S1","S2","S3","S16")]="NorthCentralChina"
Pparva.demes$pop[Pparva.demes$pop %in% c("S13","S14","S15","S17")]="NorthChina"
Pparva.demes$pop[Pparva.demes$pop %in% c("S4","S6")]="AdmixedCoastalChina"
Pparva.demes$pop[Pparva.demes$pop %in% c("S9","S18","S19","S20")]="SouthChina"
Pparva.demes$pop[Pparva.demes$pop %in% c("S10","S11")]="AdmixedContinentalChina"

# Get pop names (as factors) in the right order
Pparva.demes$pop=as.factor(Pparva.demes$pop)
Pparva.demes$pop=factor(Pparva.demes$pop,levels=c("WesternEurope","EasternEurope","Ira","Jap","AdmixedContinentalChina","SouthChina",
                                                  "AdmixedCoastalChina","NorthChina","NorthCentralChina","Tib"))

#==========================================================
# Private alleles per deme
#==========================================================
# Count of private alleles per deme
# list of private alleles per deme





#---------------------------------------------------------
# NUMBER OF PRIVATE ALLELES PER DEME
# Define demes in a strata
load(paste(datadir,"/Pparva.3000.Rda",sep=""))
# Convert populations to demes
Pparva_demes = Pparva
pop(Pparva_demes) = Pparva.demes$pop
pop(Pparva_demes)

# Use the 'poppr' package and private_alleles() on a Genind object
# Return a matrix with the count of private alleles per population and loci
# library(poppr)
# private = private_alleles(Pparva)
# Analyze private alleles based on the deme of interest
private_alleles(Pparva_demes, form = allele ~ ., level = "population")
private = as.data.frame(private_alleles(Pparva_demes, form = allele ~ .,report = "table", level = "population"))
private
# private_alleles count the number of occurences of the private allele, hence it can be different of 0/1
private[private > 0] = 1
demes.private = apply(private, 1, sum)
demes.private
sum(demes.private)
# 1,452 private alleles out of 3,000 loci
length(locNames(Pparva_demes))
length(unlist(alleles(Pparva_demes)))
# 3,000 loci and 5987 alleles


#---------------------------------------------------------
# IDENTIFY THE LIST OF PRIVATE ALLELES PER DEME
# in 'private', make the list of loci == 1 for each deme
# Loci names are in colNames
private
# Remove columns with no private allele (colSum == 0)
private = private[,which(colSums(private) > 0)]
list_private_loci = data.frame(genindID = colnames(private), locus = gsub("_[0-9]*.[0-9]*$", "", colnames(private)), deme = NA)
for (i in 1:nrow(list_private_loci)) {
  # print(i)
  list_private_loci$deme[i] = as.character(rownames(private)[which(private[,i] == 1)])
}
# Save results
write.table(list_private_loci, "Tables/List_private_alleles.csv", col.names = T, row.names = F, quote = F)

#---------------------------------------------------------
# WHICH ALLELE IS THE PRIVATE ALLELE?
# The alternative allele in the vcf at a given private locus must be the private allele
# Or more simply, based on teh translation table between genind locus name (01...04) and nucleotides (A, C, G, T)
# Wa can say that the private allele is the last part of the genindID (i.e. '.0x, with with = 1...4 corresponding to the nucleotide)
# Import vcf
# Pparva3000_vcf = read.vcf("data/Pparva3000.recode.vcf")
# Pparva3000_vcf
list_private_loci$private_allele = gsub("[0-9]*_[0-9]*.0", "", list_private_loci$genindID)
# Convert from numeric to nucleotide
list_private_loci$private_allele = as.character(c("A", "C", "G", "T")[as.numeric(list_private_loci$private_allele)])
# Save results with private allele information
write.table(list_private_loci, "Tables/List_private_alleles.csv", col.names = T, row.names = F, quote = F)





#---------------------------------------------------------
# IDENTIFY THE LIST OF PRIVATE SNPs PER DEME
# Pparva.demes is a hierfstat object with demes instead of sampling location
summary(Pparva.demes)
Pparva.demes$pop
# Estimate relative frequencies of alleles in demes
# In columns are demes and in rows are alleles
pop.freq(Pparva.demes)
# https://johnbhorne.wordpress.com/2017/07/12/identifying-private-snps-in-r/

# Relative frequencies
rel_freq = pop.freq(Pparva.demes)
length(rel_freq)
rel_freq[1]
# We use the lapply() function to “apply” a function to every item in the list. And the function we want to use needs to pull out all the zero frequencies from each population.
zeros = lapply(rel_freq, function(x) which(x==0))
zeros[1:2]
rel_freq$X6_26
# We are interested mostly in loci where n-1 out of n populations have fixed differences, so now lets iterate over the zeros list and look for loci with n-1 zeros.
n = 10 # 10 demes
private = lapply(zeros, function(x) which(length(x) == n-1))
private[2000:2005]
# In this list, loci with ten zeros register a 1 and those with less or more register 0
# So now we just need to know which loci in the ten_zeros list have a 1.
length(which(private == 1))
# 1404 private loci out of 3000 loci
private_SNPs = rel_freq[which(private == 1)]
private_SNPs[1:10]

# Another to achieve the quest of private alleles (however, not the same number of private SNPs with this method)
# Compare results of alleles in genind format (i.e. 1...4) with vcf format REF/ALT (A, C, G, T)
private_SNPs[1:10]
# Import vcf
Pparva3000_vcf = read.vcf("data/Pparva3000.recode.vcf")
Pparva3000_vcf

lapply(private_SNPs[1:10], function(x) row.names(x))
apply(Pparva3000_vcf[,c(2:3, 5, 10, 12, 13, 15, 16:18)], 2, function(x) unique(x))
# 1 means A
# 2 means C
# 3 means G
# 4 means T
# Hence we can identify which is the private allele directly in 'private_SNPs'
# The allele with the lowest frequency in the dataset
private_allele = NA
for (i in 1:length(private_SNPs)) {
  df = as.data.frame(private_SNPs[i])
  A_freq = sum(df[which(df[,1] == 1), 3], na.rm = TRUE)
  C_freq = sum(df[which(df[,1] == 2), 3], na.rm = TRUE)
  G_freq = sum(df[which(df[,1] == 3), 3], na.rm = TRUE)
  T_freq = sum(df[which(df[,1] == 4), 3], na.rm = TRUE)
  # Which is the lowest non-zero frequency?
  # Vector of frequencies with 0 replaced by NA
  freq = c(A_freq, C_freq, G_freq, T_freq)
  freq[which(freq == 0)] = NA
  # Translate 1...4 into nucleotide
  private_allele[i] = c("A", "C", "G", "T")[which.min(freq)]
}



#==========================================================
# Compute summary statistics with demes
#==========================================================

# convert the genind to a hierfstat object
Pparva.hier = genind2hierfstat(Pparva_demes)
length(levels(Pparva_demes$pop))

basic.genestats = as.data.frame(matrix(nrow = 10, ncol = 13))
colnames(basic.genestats) = c("Population", "N","Allelic_Richness","Private","Monomorphic_Loci","HE","HE.sem","HO","HO.sem","Fst","Fst.sem","Fis","Fis.sem")
# Get population names
basic.genestats$Population = as.character(levels(Pparva_demes$pop))

#---------------------------------------------------------
# Range — Asia or Europe
# basic.genestats$Range = c(rep("Europe",7),"Asia","Europe",rep("Asia",16),"Europe","Asia",rep("Europe",2))

#---------------------------------------------------------
# N — number of sampled individuals
basic.genestats$N = table(Pparva.hier$pop)

#---------------------------------------------------------
# AR - mean Allelic Richness per population, corrected for sample size
# Estimation by rarefaction
# Hierfstat
basic.genestats$Allelic_Richness = round(colMeans(allelic.richness(Pparva.hier)$Ar, na.rm = TRUE), digits = 3)

#---------------------------------------------------------
# NUMBER OF PRIVATE ALLELES PER POPULATION
# basic.genestats$Private=apply(private_alleles(Pparva),1,sum)
# basic.genestats$Private = apply(private_alleles(Pparva,form = alleles ~ ., report = "table",
# level = "population"),1,sum)
# private_alleles(Pparva, form = allele ~ ., level = "population")
private = as.data.frame(private_alleles(Pparva_demes, form = allele ~ .,report = "table", level = "population"))
private
# private_alleles count the number of occurences of the private allele, hence it can be different of 0/1
private[private > 0] = 1
basic.genestats$Private = apply(private, 1, sum)
sum(private)

#---------------------------------------------------------
# PROPORTION OF LOCI WITH NULL GENE DIVERSITY (I.E. MONOMORPHIC LOCI) --> Hs
basic.genestats$Monomorphic_Loci=round(apply(basic.stats(Pparva.hier)$Hs,2,function(x) sum(x==0,na.rm=TRUE)/length(x)),digits=3)       



#---------------------------------------------------------
# HE — mean expected heterozygosity across polymorphic loci +- s.e.m.
#basic.genestats$HE=round(Hs(Pparva),digits=3)

# Bootstrap of HE - random resampling of individuals within populations
n.pop = seppop(Pparva_demes)
for (i in 1:10) {
  cat("--------------------------\npop",i,"\n--------------------------\n")
  boots=numeric(100)
  for (j in 1:100){
    cat("pop",i,"boot",j,"\n")
    temp=c()
    subs=n.pop[[i]]
    temp=summary(subs[sample(1:nrow(n.pop[[i]]@tab),replace=TRUE),])$He
    boots[j]=mean(temp[temp>0],na.rm=TRUE)
    # boots[j]=mean(summary(Pparva[which(Pparva$pop==basic.genestats$Population[i]),sample(1:4216,replace=TRUE)])$He)
    rm(temp)
  }
  basic.genestats$HE[i]=round(mean(boots),digits=3)
  basic.genestats$HE.sem[i]=round(sd(boots),digits=3)
  print(round(mean(boots),digits=3))
  print(round(sd(boots),digits=3))
}

barplot(basic.genestats$HE)


#---------------------------------------------------------
# HO — mean observed individual heterozygosity across polymorphic loci +- s.e.m.
# Bootstrap of Ho - random resampling individuals within populations
for (i in 1:10) {
  cat("--------------------------\npop",i,"\n--------------------------\n")
  boots=numeric(100)
  for (j in 1:100){
    cat("pop",i,"boot",j,"\n")
    temp=c()
    subs=n.pop[[i]]
    temp=summary(subs[sample(1:nrow(n.pop[[i]]@tab),replace=TRUE),])$Hobs
    boots[j]=mean(temp[temp>0],na.rm=TRUE)
    # boots[j]=mean(summary(Pparva[which(Pparva$pop==basic.genestats$Population[i]),sample(1:4216,replace=TRUE)])$He)
    rm(temp)
  }
  basic.genestats$HO[i]=round(mean(boots),digits=3)
  basic.genestats$HO.sem[i]=round(sd(boots),digits=3)
}

barplot(basic.genestats$HO)


#---------------------------------------------------------
# Weir and Cockerham Fst per population  +- s.e.m.
# Compute Weir and Cockerham Fstat
# As the mean Fst per population in:
# - intercontinental range, to compare native and invasive populations in a single metapopulation
# - within continental range (Asia OR Europe) for estimating differentiation with geographical neighbours
(Pparva.pairwiseWC=pairwise.WCfst(Pparva.hier))
basic.genestats$Fst=round(apply(Pparva.pairwiseWC,2, function(x) mean(x, na.rm=TRUE)),digits=3) # Mean Fst
# Bootstrap of Fst s.e.m. with 1,000 bootstraps, by randomly resampling Fst from pairwise Fst list for a given population
sem=c()
for (j in 1:10) {
  print(j)
  boots=numeric(1000)
  for (i in 1:1000){
    boots[i]=mean(sample(Pparva.pairwiseWC,replace=TRUE),na.rm=TRUE)
  }
  sem[j]=sd(boots)# Fst s.e.m
}
basic.genestats$Fst.sem=round(sem,digits=3)

# # Asia
# Pparva.hier.native=genind2hierfstat(native)
# (Pparva.pairwiseWC.native=pairwise.WCfst(Pparva.hier.native))
# 
# # Europe
# Pparva.hier.invasive=genind2hierfstat(invasive)
# (Pparva.pairwiseWC.invasive=pairwise.WCfst(Pparva.hier.invasive))

#----------------------------------------------------------
# Fis per population +- s.e.m.
# with 1000 bootstraps
boot=boot.ppfis(dat=Pparva.hier,nboot=1000,quant=c(0.025,0.975),diploid=TRUE)
basic.genestats$Fis=round(apply(boot$fis.ci,1, mean),digits=3) # Mean
basic.genestats$Fis.sem=round((as.numeric(unlist(boot$fis.ci[2]))-as.numeric(unlist(boot$fis.ci[1])))/3.92,digits=3) # Standard error


########
# Sort populations by range (first Asia, the Europe)
# basic.genestats=basic.genestats[c(8,10,20,22:25,11:19,21,27,1:7,9,26,28,29),]

# And then export table:
write.csv(basic.genestats,file="Tables/basic.genestats.DEMES.csv",row.names = FALSE)

#----------------------------------------------------------
# Format in publication mode
basic.genestats.formatted = basic.genestats[,1:5]
# Paste estimates with s.e.m. in a readable format
basic.genestats.formatted$He = paste(basic.genestats$HE, " (±", basic.genestats$HE.sem, ")", sep = "" )
basic.genestats.formatted$Ho = paste(basic.genestats$HO, " (±", basic.genestats$HO.sem, ")", sep = "" )
basic.genestats.formatted$Fst = paste(basic.genestats$Fst, " (±", basic.genestats$Fst.sem, ")", sep = "" )
basic.genestats.formatted$Fis = paste(basic.genestats$Fis, " (±", basic.genestats$Fis.sem, ")", sep = "" )
write.csv(basic.genestats.formatted,file = "Tables/basic.genestats.formatted.DEMES.csv", row.names = FALSE)







#==========================================================
# Export filtered dataset from Genind to VCF
#==========================================================
# vcftools can filter a vcf per list of loci and list of individuals

# SITE ID FILTERING
# --snp <string>
#   Include SNP(s) with matching ID (e.g. a dbSNP rsID). This command can be used multiple times in order to include more than one SNP.
# --snps <filename>
# --exclude <filename>
#   Include or exclude a list of SNPs given in a file. The file should contain a list of SNP IDs (e.g. dbSNP rsIDs), with one ID per line. No header line is expected.

# INDIVIDUAL FILTERING OPTIONS
# These options are used to include or exclude certain individuals from any analysis being performed by the program.
# --indv <string>
# --remove-indv <string>
# Specify an individual to be kept or removed from the analysis. This option can be used multiple times to specify multiple individuals. If both options are specified, then the "--indv" option is executed before the "--remove-indv option".
# --keep <filename>
# --remove <filename>
# Provide files containing a list of individuals to either include or exclude in subsequent analysis. Each individual ID (as defined in the VCF headerline) should be included on a separate line. If both options are used, then the "--keep" option is executed before the "--remove" option.
# When multiple files are provided, the union of individuals from all keep files subtracted by the union of individuals from all remove files are kept. No header line is expected.

# In the vcf, loci identified by the three first columns #CHROM	POS	ID
# In the genind object, loci names are unique
locNames(Pparva)
# Names corresponding to the vcf (removed suffix)
gsub("_[0-9]*", "", locNames(Pparva))
list.loci = gsub("_[0-9]*", "", locNames(Pparva))
write(list.loci, "Data/list.loci")

# In the vcf, individuals are in columns
# In genind, use the function to get the accessor
indNames(Pparva)
list.ind = indNames(Pparva)
length(list.ind)
# 468 individuals to retain
# Beware of duplicated individuals
# That are individuals...
# List of individuals in the vcf
vcf.ind = scan(file = "Data/vcf.ind.txt", what = character(), sep = ",")
vcf.ind
# Make the list of loci conform to both vcf and genind
# Identify duplicated individuals and replace in list.ind by the name of the first duplicate (e.g. replace "S11_25" by "S11_25a")
# Individuals matching pattern "_[0-9]*a$"
duplicates = vcf.ind[grep("_[0-9]*a$", vcf.ind)]
duplicates
for (i in 1:length(list.ind)) {
  if (paste(list.ind[i], "a", sep = "") %in% duplicates) {
    list.ind[i] = paste(list.ind[i], "a", sep = "")
  }
}
list.ind
write(list.ind, "Data/list.ind")

# Filter loci and individuals with vcftools
system("vcftools --vcf Data/batch_1.vcf --out Data/Pparva3000 --snps Data/list.loci --keep Data/list.ind --recode")

# Produce a summary of the dataset in vcf



#==========================================================
# Basic statistics per locus
#----------------------------------------------------------

#----------------------------------------------------------
# Fis per locus
#----------------------------------------------------------


#----------------------------------------------------------
# Fst per locus
#----------------------------------------------------------



#----------------------------------------------------------
# Distribution of mean Fst per locus
#----------------------------------------------------------


#----------------------------------------------------------
# Identify loci showing a higher Fst than others,that could indicate adaptive divergence
#----------------------------------------------------------








#==========================================================
# Distribution of allele frequency
#----------------------------------------------------------

#----------------------------------------------------------
# Distribution of MAF: rare alleles (<0.005)
# Minor Allele frequencies



#==========================================================
# Linkage disequilibrium
#==========================================================




#==========================================================
#   Hardy Weinberg equilibrium
#==========================================================













#==========================================================
# THE END