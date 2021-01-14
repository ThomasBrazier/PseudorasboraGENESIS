############################################################################
#   INRA - ANR project Pseudorasbora
#
#       STRUCTURE input files
#
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com
# MSc internship

# Supervisor: Scott McCairns
#             INRA UMR ESE


# DESCRIPTION
# This script produces input files for STRUCTURE analyses
# It sampled populations of interest and assemble them in a correct format for STRUCTURE

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
# STRUCTURE input files
#==========================================================

# Pparva.structure is formatted to produce STRUCTURE data set
# See STRUCTURE.R for the complete procedure...
Pparva.structure=read.table("Data/STRUCTURE/Pparva.structure.txt",sep="\t")


#--------------------------------
# WRITE DATA SET FOR STRUCTURE
#--------------------------------
# Rewrite data file
# 468 individuals, so 936 lines in the file for verification
write.table(Pparva.structure,"Data/STRUCTURE/1a_Allpops/100a_Pparva_structure_allpops.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
# Add POPFLAG
# Asian populations are popflagged 1 whereas european populations (ancestry to predict) are popflagged 0
popflag=Pparva.structure[,2] # popflag contain population number
popflag[which(popflag %in% c(1:7,9,26,28,29))]=0
popflag[which(popflag != 0)]=1
Pparva.structure.popflag=cbind(Pparva.structure[,1:2],
                               popflag,
                               Pparva.structure[,3:ncol(Pparva.structure)])
write.table(Pparva.structure.popflag,"Data/STRUCTURE/1b_Allpops/100b_Pparva_structure_allpops.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

###### LABELS for Distruct
# Labels below figure were simply population number in data set
# Generate a file of labels atop figure (Population names)
# A first column of population's number (number in dataset) and a second column with the corresponding population name
labels_atop=data.frame(seq(1,length(levels(Pparva@pop))),levels(Pparva@pop))
write.table(labels_atop,"Data/STRUCTURE/labels_atop.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)


#--------------------------------
# Only populations of native range
# i.e. individuals in 'temp' whom populations (col. 2) were included in labelPops[c(8,11:28,31:32),2]
## Write data file
# 300 individuals, so 600 lines in the file for verification
write.table(as.data.frame(Pparva.structure[which(Pparva.structure[,2] %in% c(8,10:25,27)),]),"Data/STRUCTURE/2_Native/200_Pparva_structure_native.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(as.data.frame(Pparva.structure[which(Pparva.structure[,2] %in% c(8,10:25,27)),2]),"Data/STRUCTURE/2_Native/labelPops_native.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

# # Just a first batch run of 8 runs for K=1 to infer lambda
# structureLauncher(project="2_Native",path="/home/users/tbrazier/STRUCTURE/2_Native",niter=c(1,8),Kmin=1,Kmax=1,
#                   popsize=nrow(Pparva.structure[which(Pparva.structure[,2] %in% c(8,10:25,27)),])/2)
# # Then relaunch with the inferred value of lambda
# structureLauncher(project="2_Native",path="/home/users/tbrazier/STRUCTURE/2_Native",niter=c(1,8),Kmin=2,Kmax=16,
#                   popsize=nrow(Pparva.structure[which(Pparva.structure[,2] %in% c(8,10:25,27)),])/2)
# structureLauncher(project="2_Native",path="/home/users/tbrazier/STRUCTURE/2_Native",niter=c(9,15),Kmin=2,Kmax=16,
#                   popsize=nrow(Pparva.structure[which(Pparva.structure[,2] %in% c(8,10:25,27)),])/2)
# # 5 more replicates for K=2-16 and 100,000 iterations:
# structureLauncher(project="2_Native",path="/home/users/tbrazier/STRUCTURE/2_Native",niter=c(16,20),Kmin=2,Kmax=16,
#                   popsize=nrow(Pparva.structure[which(Pparva.structure[,2] %in% c(8,10:25,27)),])/2)
# # then 15 replicates for K=17-21 and 100,000 iterations:
# structureLauncher(project="2_Native",path="/home/users/tbrazier/STRUCTURE/2_Native",niter=c(1,15),Kmin=17,Kmax=21,
# popsize=nrow(Pparva.structure[which(Pparva.structure[,2] %in% c(8,10:25,27)),])/2)


# then 20 replicates for K=1-21 for 100,000 BURNIN + 100,000 iterations:
structureLauncher(project="2_Native",path="/home/users/tbrazier/STRUCTURE/2_Native",niter=c(1,20),Kmin=1,Kmax=21,
                  popsize=nrow(Pparva.structure[which(Pparva.structure[,2] %in% c(8,10:25,27)),])/2)

# HIERARCHICAL PROCEDURE
# STRUCTURE on a subset of populations in native area: explain genetic structure of the large Pan-asian admixed zone
write.table(as.data.frame(Pparva.structure[which(Pparva.structure[,2] %in% c(10:17,20,22:24,27)),]),"Data/STRUCTURE/2b1_Native/200b_Pparva_structure_native.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
structureLauncher(project="2b1_Native",path="/home/users/tbrazier/STRUCTURE/2b1_Native",niter=c(1,20),Kmin=1,Kmax=16,
                  popsize=nrow(Pparva.structure[which(Pparva.structure[,2] %in% c(10:17,20,22:24,27)),])/2)

#################7
# There was a lack of convergence for 100,000 iterations, but it was suspected that the K value was between 3, 6 and 7
# So, we decided to relaunch STRUCTURE for K=6 and 500,000 BURNIN + 500,000 iterations (previously, no sign of convergence after 100,000 iterations)
structureLauncher(project="21_Native_1MillionIterations",path="/home/users/tbrazier/STRUCTURE/21_Native_1MillionIterations",niter=c(1,16),Kmin=4,Kmax=8,
                  popsize=nrow(Pparva.structure[which(Pparva.structure[,2] %in% c(8,10:25,27)),])/2) # 8 cores on local server --> 2*8=16 replicates
structureLauncher(project="21_Native_1MillionIterations",path="/home/users/tbrazier/STRUCTURE/21_Native_1MillionIterations",niter=c(1,16),Kmin=1,Kmax=3,
                  popsize=nrow(Pparva.structure[which(Pparva.structure[,2] %in% c(8,10:25,27)),])/2) # 8 cores on local server --> 2*8=16 replicates


#-----------------------------------
# Only populations of invasive range
## Write data file
# 168 individuals, so 336 lines in the file for verification
write.table(as.data.frame(Pparva.structure[which(Pparva.structure[,2] %in% c(1:7,9,26,28:29)),]),"Data/STRUCTURE/3_Invasive/300_Pparva_structure_invasive.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
write.table(as.data.frame(Pparva.structure[which(Pparva.structure[,2] %in% c(1:7,9,26,28:29)),2]),"Data/STRUCTURE/3_Invasive/labelPops_invasive.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)

structureLauncher(project="3_Invasive",path="/home/users/tbrazier/STRUCTURE/3_Invasive",niter=c(1,20),Kmin=1,Kmax=14,
                  popsize=nrow(Pparva.structure[which(Pparva.structure[,2] %in% c(1:7,9,26,28:29)),2])/2)

# HIERARCHICAL PROCEDURE
# Western Europe  showed a pattern of strong admixture (spurious clusters suspected?) with 2 blended clusters
# STRUCTURE on a subset of populations in invasive area: explain genetic structure of the Western Europe admixed zone
# n=119
write.table(as.data.frame(Pparva.structure[which(Pparva.structure[,2] %in% c(1,2,5,7,9,26,29)),]),"Data/STRUCTURE/3a1_Invasive/300a1_Pparva_structure_invasive.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
structureLauncher(project="3a1_Invasive",path="/home/users/tbrazier/STRUCTURE/3a1_Invasive",niter=c(1,20),Kmin=1,Kmax=10,
                  popsize=nrow(Pparva.structure[which(Pparva.structure[,2] %in% c(1,2,5,7,9,26,29)),])/2)




#-----------------------------------
# WORLDWIDE ANALYSES
# We ran this intercontinental cluster analysis to see if some source populations were closest to invasive populations than their own neighbors in the native area
# It was also interesting to estimates admixture within populations in order to parameterize, DIY ABC

# 20 replicates for K=1-32 (all sampled sites + 3) for 100,000 BURNIN + 100,000 iterations:
structureLauncher(project="1a_AllPops",path="/home/users/tbrazier/STRUCTURE/1a_AllPops",niter=c(1,20),Kmin=1,Kmax=32,
                  popsize=nrow(Pparva.structure)/2)


# 20 replicates for K=1-21 (all sampled sites in Asia (reference) + 3) for 100,000 BURNIN + 100,000 iterations:
structureLauncher(project="1b_AllPops",path="/home/users/tbrazier/STRUCTURE/1b_AllPops",niter=c(1,20),Kmin=1,Kmax=21,
                  popsize=nrow(Pparva.structure)/2)

# It may be smarter to define populations in Asia as demes and not as sampled populations. 
# 20 replicates for K=1-13 (10 putative demes in Asia (reference) + 3) for 100,000 BURNIN + 100,000 iterations:
#-------------------------- NATIVE
# POP1 S4 S6
# POP2 Jap
# POP3 S10
# POP4 S13 S14 S15 S17
# POP5 S1 S2 S16
# POP6 S9 S18
#------------------------ INVASIVE
# POP7 Tib
# POP8 S19 S20
#------------------------ ADMIXTURE
# POP9 S3
# POP10 S11
Pparva.structure.demes=Pparva.structure.popflag
# Replace sampled sites in Asia by putative deme (putative source population): the deme number is the pop number of the first population in the vector
Pparva.structure.demes[which(Pparva.structure.demes$V2 %in% c("S4","S6")),3]="23"
Pparva.structure.demes[which(Pparva.structure.demes$V2 %in% c("S10")),3]="11"
Pparva.structure.demes[which(Pparva.structure.demes$V2 %in% c("S13", "S14" ,"S15", "S17")),3]="13"
Pparva.structure.demes[which(Pparva.structure.demes$V2 %in% c("S1", "S2", "S16")),3]="10"
Pparva.structure.demes[which(Pparva.structure.demes$V2 %in% c("S9", "S18")),3]="25"
Pparva.structure.demes[which(Pparva.structure.demes$V2 %in% c("S19", "S20")),3]="19"
# Write the data set for STRUCTURE
write.table(Pparva.structure.demes,"Data/STRUCTURE/1c_Allpops/100c_Pparva_structure_demes.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=FALSE)
structureLauncher(project="1c_AllPops",path="/home/users/tbrazier/STRUCTURE/1c_AllPops",niter=c(1,20),Kmin=1,Kmax=13,
                  popsize=nrow(Pparva.structure.demes)/2)
