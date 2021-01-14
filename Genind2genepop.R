############################################################################
#   INRA - ANR project Pseudorasbora
#
#       FILE CONVERSION - GENIND 2 GENEPOP
#
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com
# MSc internship

# Supervisor: Scott McCairns
#             INRA UMR ESE


# DESCRIPTION
# This script performs file conversion for the whole data analysis pipeline
# Currently:
# - genind 2 genepop


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

# convert the genind to a hierfstat object
Pparva.hier=genind2hierfstat(Pparva)



#==========================================================
# Format in genepop file
#==========================================================
# Genepop format is required for multiple R packages for genetic diversity (e.g. diveRsity)

# The Input File Format for Genepop(WWW) and LinkDos(WWW) :
#   
#   First line : Any characters. Use this line to store information about you data.
# Second line : the name of the first locus.
# Third line : the name of the second locus (if needed)
# ETC...
# 
# Alternatively, separate loci with commas (or a comma+space) on the same line.
# Line N+1 : the name of the Nth locus.
# Line N+2 : Type the word "POP" (or "Pop", or "pop").
# Line N+3 : And example is given below :
#   
#   ind#001fem, 0101 0202 0000 0410
genind.df=genind2df(Pparva)
# !!!!! Missing data should be indicated with 00
# Missing data (NA) was transformed in '0000'
genind.df[is.na(genind.df)]="0000"

# ............
# Open connection to new genepop files: 3 dataset (global, Asian and European) have been printed at once
file.create("Data/genepop.txt")
file.create("Data/genepop_native.txt")
file.create("Data/genepop_invasive.txt")
# print first a title line
write.table("Phylogeography of Pseudorasbora parva","Data/genepop.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
write.table("Phylogeography of Pseudorasbora parva: Native area","Data/genepop_native.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
write.table("Phylogeography of Pseudorasbora parva: Invasive area","Data/genepop_invasive.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)

# Locus names: They may be given one per line, or on the same line but separated by commas.
# write(paste(colnames(genind.df)[-1],collapse ="\n"),file="Data/genepop.txt",append=TRUE)
# write(paste(colnames(genind.df)[-1],collapse ="\n"),file="Data/genepop_native.txt",append=TRUE)
# write(paste(colnames(genind.df)[-1],collapse ="\n"),file="Data/genepop_invasive.txt",append=TRUE)

write(paste(colnames(genind.df)[-1],collapse =" , "),file="Data/genepop.txt",append=TRUE)
write(paste(colnames(genind.df)[-1],collapse =" , "),file="Data/genepop_native.txt",append=TRUE)
write(paste(colnames(genind.df)[-1],collapse =" , "),file="Data/genepop_invasive.txt",append=TRUE)

# Then, one individual per line, with the following nomenclature 'ind_name , locus1 locus2 ...'
# Each sample from a different geographical original is declared by a line with a POP statement.

# Print each population at a time
# Lists of population names with their associated coordinates
native.names=c("Jap", "S1", "S2", "S3", "S4", "S6","S9",
               "S10", "S11", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "Tib")
invasive.names=c("Aus", "Bel", "Bul1", "Bul2", "Hun", "Ira", "Ita", "Pol", "Spa", "Tur", "UK")

# Native populations were added first(Asia)
for (i in 1:length(native.names)) {
  write("POP",file="Data/genepop.txt",append=TRUE)
  write("POP",file="Data/genepop_native.txt",append=TRUE)
  print(native.names[i])
  subset=genind.df[which(genind.df$pop %in% native.names[i]),]
  for (s in 1:nrow(subset)) {
    write(paste(subset[s,1],paste(subset[s,-1],collapse=" "),sep=", "),file="Data/genepop.txt",append=TRUE)
    write(paste(subset[s,1],paste(subset[s,-1],collapse=" "),sep=", "),file="Data/genepop_native.txt",append=TRUE)
  }
}

# Then invasive populations were added (Europe)
for (i in 1:length(invasive.names)) {
  write("POP",file="Data/genepop.txt",append=TRUE)
  write("POP",file="Data/genepop_invasive.txt",append=TRUE)
  print(invasive.names[i])
  subset=genind.df[which(genind.df$pop %in% invasive.names[i]),]
  for (s in 1:nrow(subset)) {
    write(paste(subset[s,1],paste(subset[s,-1],collapse=" "),sep=", "),file="Data/genepop.txt",append=TRUE)
    write(paste(subset[s,1],paste(subset[s,-1],collapse=" "),sep=", "),file="Data/genepop_invasive.txt",append=TRUE)
  }
}

#---------------------------------
# GENEPOP files with demes
genind.df=genind2df(Pparva)
# !!!!! Missing data should be indicated with 00
# Missing data (NA) was transformed in '0000'
genind.df[is.na(genind.df)]="0000"
# Open connection to new genepop files: 2 dataset (Asian and European)
file.create("Data/genepop_native_demes.txt")
file.create("Data/genepop_invasive_demes.txt")
# print first a title line
write.table("Phylogeography of Pseudorasbora parva: demes of Native area","Data/genepop_native_demes.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
write.table("Phylogeography of Pseudorasbora parva: demes of Invasive area","Data/genepop_invasive_demes.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)

# Locus names: They may be given one per line, or on the same line but separated by commas.
write(paste(colnames(genind.df)[-1],collapse =" , "),file="Data/genepop_native_demes.txt",append=TRUE)
write(paste(colnames(genind.df)[-1],collapse =" , "),file="Data/genepop_invasive_demes.txt",append=TRUE)

# Then, one individual per line, with the following nomenclature 'ind_name , locus1 locus2 ...'
# Each sample from a different geographical original is declared by a line with a POP statement.

# Print each population at a time
# Lists of demes
genind.df.demes = genind.df
genind.df.demes[,1]=as.character(genind.df.demes[,1])
genind.df.demes[which(genind.df.demes$pop %in% c("S4","S6")),1]="AdmixedCoastalChina"
genind.df.demes[which(genind.df.demes$pop %in% c("S10","S11")),1]="AdmixedContinentalChina"
genind.df.demes[which(genind.df.demes$pop %in% c("S13","S14", "S15", "S17")),1]="NorthChina"
genind.df.demes[which(genind.df.demes$pop %in% c("S1", "S2", "S3", "S16")),1]="NorthCentralChina"
genind.df.demes[which(genind.df.demes$pop %in% c("S9", "S18", "S19", "S20")),1]="SouthChina"

genind.df.demes[which(genind.df.demes$pop %in% c("Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK")),1]="WesternEurope"
genind.df.demes[which(genind.df.demes$pop %in% c("Bul1", "Bul2", "Tur", "S20")),1]="EasternEurope"

native.demes=c("NorthChina", "NorthCentralChina", "AdmixedContinentalChina", "AdmixedCoastalChina",
               "SouthChina", "Jap", "Tib")
invasive.demes=c("WesternEurope", "EasternEurope", "Ira")

# Native demes
for (i in 1:length(native.demes)) {
  write("POP",file="Data/genepop_native_demes.txt",append=TRUE)
  print(native.demes[i])
  subset=genind.df.demes[which(genind.df.demes$pop %in% native.demes[i]),]
  for (s in 1:nrow(subset)) {
    write(paste(subset[s,1],paste(subset[s,-1],collapse=" "),sep=", "),file="Data/genepop_native_demes.txt",append=TRUE)
  }
}

# Then invasive demes
for (i in 1:length(invasive.demes)) {
  write("POP",file="Data/genepop_invasive_demes.txt",append=TRUE)
  print(invasive.demes[i])
  subset=genind.df.demes[which(genind.df.demes$pop %in% invasive.demes[i]),]
  for (s in 1:nrow(subset)) {
    write(paste(subset[s,1],paste(subset[s,-1],collapse=" "),sep=", "),file="Data/genepop_invasive_demes.txt",append=TRUE)
  }
}

## Initialising summary file
#  Successive results will be append to this file
write.table("# SUMMARY #\n\n","summary.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)







#==========================================================
# THE END