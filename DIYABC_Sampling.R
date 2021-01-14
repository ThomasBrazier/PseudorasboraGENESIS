############################################################################
#   INRA - ANR project Pseudorasbora
#
#       DIY ABC -- Sampling
#
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com
# MSc internship

# Supervisor: Scott McCairns
#             INRA UMR ESE


# DESCRIPTION
# In this script, we performed the sampling of populations for DIY ABC input dataset.

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
# STATS
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


#---------------------------------
# Format in DIY ABC file
# see page 91 of the DIY ABC user manual
load("Data/Pparva.3000.Rda") # Make sure that the dataset of 3,000 loci is used
          # genind.df=genind2df(Pparva)
          # # genind.df=genind2df(Pparva[loc=1:200])
          # # !!!!! Missing data should be indicated with 00
          # # Missing data (NA) was transformed in '0000'
          # genind.df[is.na(genind.df)]="0000"
          # 
          # # ............
          # # Open connection to new snp files: 3 dataset (global, Asian and European) have been printed at once
          # file.create("Data/ABC/ABC.txt")
          # file.create("Data/ABC/ABC_native.txt")
          # file.create("Data/ABC/ABC_invasive.txt")
          # # print first a title line with parameters NM=x NF (number of females per male, sex ratio) and required MAF to take in account (minimum allele frequency)
          # write.table("Phylogeography of Pseudorasbora parva <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
          # write.table("Phylogeography of Pseudorasbora parva in Native area <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_native.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
          # write.table("Phylogeography of Pseudorasbora parva in Invasive area <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_invasive.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
          # 
          # abc.df=cbind(names(Pparva$tab[,1]),rep(9,nrow(genind.df)),genind.df)
          # abc.df[,1]=as.character(abc.df[,1])
          # 
          # # one individual per line, with the following nomenclature:
          # # as many lines as there are genotyped individuals, with the code-name of the individual,
          # # a letter (M or F) indicating its sex ('9' for missing data),
          # # a code-name for its population and the values (0, 1 or 2) of the number of the (arbitrarily chosen) reference allele at each SNP locus. For instance in the case autosomal diploid SNP loci,
          # # we have 0 = homozygous genotype for the non reference allele, 1 = heterozygous genotype for the reference allele, 2 = homozygous genotype for the reference allele.
          # # And missing data is noted '9'
          # # It is worth noting that for autosomal haploid loci (denoted H), as well as for mitochondrial loci (denoted M) and Y-linked loci (denoted Y), the SNP genotypes will be 0 or 1.
          # 
          # # Transform the genotype data in the required format (0, 1, 2 or 9)
          # # Identify a reference allele (arbitrarily chosen) for each locus
          # for (i in 4:ncol(abc.df)) {
          #   # Reduce a locus to all unique alleles
          #   all.i=unique(c(substr(unique(abc.df[,i]),1,2),substr(unique(abc.df[,i]),3,4)))
          #   all.1=all.i[all.i!="00"][1] # allele 1
          #   all.2=all.i[all.i!="00"][2] # allele 2
          #   # Homozygote of allele 1 is coded '0'
          #   abc.df[abc.df[,i]==paste(all.1,all.1,sep=""),i]="0"
          #   # Heterozygote is coded '1'
          #   abc.df[abc.df[,i]==paste(all.1,all.2,sep="") | abc.df[,i]==paste(all.2,all.1,sep=""),i]="1"
          #   # Homozygote of allele 2 is coded '2'
          #   abc.df[abc.df[,i]==paste(all.2,all.2,sep=""),i]="2"
          #   
          #   rm(all.1)
          #   rm(all.2)
          #   rm(all.i)
          # }
          # # Missing data is expressed as '9'
          # abc.df[abc.df=="0000"]="9"
          # 
          # # remove loci with missing data
          # #abc.df=abc.df[,-c(632,908,1659,1602,2115)] # remove loci 629, 905, 1599, 1656, 2112 (+3 columns) --> manual removal too long
          # # We automated it...
          # # For each locus, compute a table per population: if 0 for two alleles in one pop, the locus was added to the removal.idx
          # removal.idx=c()
          # for (j in 4:ncol(abc.df)) {
          #   print(j)
          #   tab=table(abc.df[,c(3,j)])
          #   for (l in 1:nrow(tab)) {
          #     if (sum(tab[l,-ncol(tab)])==0) {
          #       removal.idx=c(removal.idx,j)
          #     }
          #   }
          # }
          # rm(tab)
          # # remove all incorrect loci (at least one population with only missing data) at once
          # abc.df=abc.df[,-unique(removal.idx)] # 195 loci removed, 2,805 loci retained
          # rm(removal.idx)
          # 
          # # The second line is starting with the three keywords IND SEX POP, separated by at least one space, followed by as many letters as SNP loci,
          # # the letter giving the location of the locus as above
          # # (< A > for autosomal diploid loci, < H > for autosomal haploid loci,
          # # < X > for X-linked (or haplo-diploid) loci,
          # # < Y > for Y-linked loci and < M > for mitochondrial loci). Letters are separated by a single space.
          # N.loci=ncol(abc.df)-3
          # # All loci (N=2112) were autosomal diploid loci < A >
          # write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC.txt",append=TRUE)
          # write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_native.txt",append=TRUE)
          # write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_invasive.txt",append=TRUE)
          # 
          # 
          # 
          # # Lists of population names
          # native.names=c("Jap", "S1", "S2", "S3", "S4", "S6","S9",
          #                "S10", "S11", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "Tib")
          # invasive.names=c("Aus", "Bel", "Bul1", "Bul2", "Hun", "Ira", "Ita", "Pol", "Spa", "Tur", "UK")
          # 
          # for (i in 1:nrow(abc.df)) {
          #   # Add to global data set
          #   write(paste(abc.df[i,],collapse=" "),file="Data/ABC/ABC.txt",append=TRUE)
          #   
          #   # Add to native (Asian)
          #   if (genind.df[i,1] %in% native.names) {
          #     write(paste(abc.df[i,],collapse=" "),file="Data/ABC/ABC_native.txt",append=TRUE)
          #   }
          #   
          #   # Add to invasive (Europe)
          #   if (genind.df[i,1] %in% invasive.names) {
          #     write(paste(abc.df[i,],collapse=" "),file="Data/ABC/ABC_invasive.txt",append=TRUE)
          #   }
          # }

# save(abc.df,file="Data/abc.df.Rda")

# RE-IMPORT ABC.DF
load("Data/abc.df.Rda")
N.loci=ncol(abc.df)-3


#-------------------------------------------------
# Custom data set for scenarios

# REMIND, FOR ALL SCENARIOS
# Computation time is precious, so
# KISS (Keep It Simple, Stupid !)

# test simple hypothesis:
# 1/ Origin of the Japanese population
# 2/ Source of 3 invasive demes in Europe
#       In 3 separated analyses 2.(a, b, c)
# 3/ Dispersal pattern in Western Europe


###############################
# 1/ Scenarios for Japan lineage phylogeography: Sampled locations at the extreme range didn't represented the pan-Asian deme
# though, they were discarded in this first step scenarios dedicated to inference of the history of the japanese deme

# KISS (https://en.wikipedia.org/wiki/KISS_principle), test simple scenario:
# We expect Japan to be either ancient divergence with mainland
# or recent introduction from East of China
# 3 pops to test North, South and Japan
abc.df.tmp=rbind(abc.df[which(abc.df$pop %in% c("Jap")),],
                 abc.df[which(abc.df$pop %in% c("S3", "S4", "S6", "S1", "S2", "S16")),]
)
# Cluster sampled locations into demes
abc.df.tmp[,3]=as.character(abc.df.tmp[,3])
abc.df.tmp[which(abc.df.tmp$pop %in% c("Jap")),3]="1"
abc.df.tmp[which(abc.df.tmp$pop %in% c("S4","S6")),3]="2"
abc.df.tmp[which(abc.df.tmp$pop %in% c("S1", "S2", "S3", "S16")),3]="3"

# Reorder the dataset, DIY ABC names pop in incoming order
abc.df.tmp = abc.df.tmp[order(abc.df.tmp$pop),]
abc.df.tmp$pop

# Print into file
file.create("Data/ABC/ABC_native_JapanOrigin_Simple.txt")
write.table("Phylogeography of Pseudorasbora parva <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_native_JapanOrigin_Simple.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_native_JapanOrigin_Simple.txt",append=TRUE)
for (i in 1:nrow(abc.df.tmp)) {
  write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_native_JapanOrigin_Simple.txt",append=TRUE)
}
rm(abc.df.tmp)



# Then, go to more complex ones...
# We had a prior that Japan were linked to the admixture zone near S4
# Although there exists an endemic lineage in Japan, this hypotheses need to be tested
# Besides, we must also test for the case of an unsampled source population

# We defined DEMES
# POP1 Jap
# POP2 S3 S4 S6 i.e. Admixed Zone
# POP3 S10 S11 i.e. continental China (admixed)
# POP4 S13 S14 S15 S17 i.e. North China
# POP5 S1 S2 S16 i.e. North Central China
# POP6 S9 S18 S19 S20 i.e. South China

# Tibet is not included, probability of a link with Japan is very low

# "Known" history of biogeography and species was used to draw the null hypothesis tree (i.e. the tree of populations history not to test)
# This is the same tree in all scenarios, except for the position of the Japanese DEME for which we search a source population

abc.df.tmp=rbind(abc.df[which(abc.df$pop %in% c("Jap")),],
                 abc.df[which(abc.df$pop %in% c("S3", "S4", "S6", "S1", "S2", "S16","S10","S11", "S13",
                                                "S14", "S15", "S17", "S9", "S18", "S19", "S20")),]
)
# Cluster sampled locations into demes
abc.df.tmp[,3]=as.character(abc.df.tmp[,3])
abc.df.tmp[which(abc.df.tmp$pop %in% c("Jap")),3]="1"
abc.df.tmp[which(abc.df.tmp$pop %in% c("S3", "S4","S6")),3]="2"
abc.df.tmp[which(abc.df.tmp$pop %in% c("S10","S11")),3]="3"
abc.df.tmp[which(abc.df.tmp$pop %in% c("S13","S14", "S15", "S17")),3]="4"
abc.df.tmp[which(abc.df.tmp$pop %in% c("S1", "S2", "S16")),3]="5"
abc.df.tmp[which(abc.df.tmp$pop %in% c("S9", "S18", "S19", "S20")),3]="6"

# Reorder the dataset, DIY ABC names pop in incoming order
abc.df.tmp = abc.df.tmp[order(abc.df.tmp$pop),]
abc.df.tmp$pop

# Print into file
file.create("Data/ABC/ABC_native_JapanOrigin.txt")
write.table("Phylogeography of Pseudorasbora parva <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_native_JapanOrigin.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_native_JapanOrigin.txt",append=TRUE)
for (i in 1:nrow(abc.df.tmp)) {
  write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_native_JapanOrigin.txt",append=TRUE)
}
rm(abc.df.tmp)


###############################
# 2.A/ Scenarios for routes of invasion in Europe: Sampled locations at the extreme range didn't represented the pan-Asian deme
# though, they were discarded in this first step scenarios dedicated to inference of routes of invasion

# S3 systematically failed in population assignment ('assignPOP')
# Besides, sample size is only 8, with high admixture rates in this particulaer population
# Hence, S3 was discarded for any subsequent analysis

# We defined DEMES
#-------------------------- NATIVE
# POP1 Jap
# POP2 S4 S6 i.e. Admixed Zone
# POP3 S10 S11 i.e. continental China (admixed)
# POP4 S13 S14 S15 S17 i.e. North China
# POP5 S1 S2 S16 i.e. North Central China
# POP6 S9 S18 S19 S20 i.e. South China

# Basically, it is the validated scenario of Japanese Origin, with the target population tested in each terminal branch
# --> 6 source pop. to test + unsampled pop.
#------------------------ INVASIVE
# POP7.A Western Europe AUS BEL HUN ITA POL SPA UK
# POP7.B Eastern Europe BUL1 BUL2 TUR
# POP7.C Iran IRA

############################
# Western European invasion: we tested Western Europe source in POP2 or POP3 or POP4 or POP5 or POP6 (or Unsampled)
# Western Europe is highly admixed, so there is a lot of plausible origins to discriminate
# Japan is not a plausible origin for Western Europe
# abc.df.tmp=rbind(abc.df[which(abc.df$pop %in% c("S4", "S6", "S1", "S2", "S16","S10","S11", "S13",
#                                                 "S14", "S15", "S17", "S9", "S18", "S19", "S20",
#                                                 "Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK")),]
# )
# # Cluster sampled locations into demes
# abc.df.tmp[,3]=as.character(abc.df.tmp[,3])
# abc.df.tmp[which(abc.df.tmp$pop %in% c("Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK")),3]="1" # Western Europe
# abc.df.tmp[which(abc.df.tmp$pop %in% c("S4","S6")),3]="4" # Admixed coastal China 
# abc.df.tmp[which(abc.df.tmp$pop %in% c("S10","S11")),3]="5" # Admixed continental China
# abc.df.tmp[which(abc.df.tmp$pop %in% c("S13","S14", "S15", "S17")),3]="2" # North China
# abc.df.tmp[which(abc.df.tmp$pop %in% c("S1", "S2", "S16")),3]="3" # North continental China
# abc.df.tmp[which(abc.df.tmp$pop %in% c("S9", "S18", "S19", "S20")),3]="6" # South China
# 
# # Reorder the dataset, DIY ABC names pop in incoming order
# abc.df.tmp = abc.df.tmp[order(abc.df.tmp$pop),]
# abc.df.tmp$pop
# 
# # Print into file
# file.create("Data/ABC/ABC_RoutesInvasionWesternEu.txt")
# write.table("Phylogeography of Pseudorasbora parva: Western Europe Origin <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_RoutesInvasionWesternEu.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
# write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_RoutesInvasionWesternEu.txt",append=TRUE)
# for (i in 1:nrow(abc.df.tmp)) {
#   write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_RoutesInvasionWesternEu.txt",append=TRUE)
# }
# rm(abc.df.tmp)



############################
# Western European invasion (SIMPLE): we tested Western Europe source in POP2 or POP3 or POP4 or POP5 or POP6 (or Unsampled)
# Western Europe is highly admixed, so there is a lot of plausible origins to discriminate
# Japan is not a plausible origin for Western Europe

# In a second thought, removing Ita, not a good replicate of the same pool of alleles 
# heritated from the original introduction, meaning that there were successive bottlenecks leading
# to those populations that we couldn't model due to complexity of the modle and low sampling(not enough loci)
# Hence, remove them because they had noise if they were not modelled
# abc.df.tmp=rbind(abc.df[which(abc.df$pop %in% c("S4", "S6", "S1", "S2", "S16","S10","S11",
#                                                 "Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK")),]
# )
abc.df.tmp=rbind(abc.df[which(abc.df$pop %in% c("S1", "S2", "S16","S10","S11",
                                                "S13", "S14", "S15", "S17",
                                                "Aus", "Bel", "Hun", "Pol", "Spa", "UK", "Ita")),]
)

# Scenario with 3 candidate native pop. (the one tested)
# Cluster sampled locations into demes
abc.df.tmp[,3]=as.character(abc.df.tmp[,3])
abc.df.tmp[which(abc.df.tmp$pop %in% c("Aus", "Bel", "Hun", "Pol", "Spa", "UK", "Ita")),3]="1" # Western Europe
abc.df.tmp[which(abc.df.tmp$pop %in% c("S1", "S2", "S16")),3]="2" # North central China, without S3 (unresolved deme assignment)
abc.df.tmp[which(abc.df.tmp$pop %in% c("S13", "S14", "S15", "S17")),3]="3" # North China
abc.df.tmp[which(abc.df.tmp$pop %in% c("S10","S11")),3]="4" # Admixed continental China
# abc.df.tmp[which(abc.df.tmp$pop %in% c("S4","S6")),3]="5" # Admixed coastal China
# Reorder the dataset, DIY ABC names pop in incoming order
abc.df.tmp = abc.df.tmp[order(abc.df.tmp$pop),]
abc.df.tmp$pop





# Scenario with 4 candidate native pop.
# Cluster sampled locations into demes
abc.df.tmp=rbind(abc.df[which(abc.df$pop %in% c("S1", "S2", "S16","S10","S11",
                                                "S13", "S14", "S15", "S17", "S4", "S6",
                                                "Aus", "Bel", "Hun", "Pol", "Spa", "UK", "Ita")),]
)
abc.df.tmp[,3]=as.character(abc.df.tmp[,3])

abc.df.tmp[which(abc.df.tmp$pop %in% c("Aus", "Bel", "Hun", "Pol", "Spa", "UK", "Ita")),3]="1" # Western Europe
abc.df.tmp[which(abc.df.tmp$pop %in% c("S1", "S2", "S16")),3]="2" # North central China, without S3 (unresolved deme assignment)
abc.df.tmp[which(abc.df.tmp$pop %in% c("S13", "S14", "S15", "S17")),3]="3" # North China
abc.df.tmp[which(abc.df.tmp$pop %in% c("S10","S11")),3]="4" # Admixed continental China
abc.df.tmp[which(abc.df.tmp$pop %in% c("S4","S6")),3]="5" # Admixed coastal China
# Reorder the dataset, DIY ABC names pop in incoming order
abc.df.tmp = abc.df.tmp[order(abc.df.tmp$pop),]
abc.df.tmp$pop

# Print into files
# file.create("Data/ABC/ABC_RoutesInvasionWesternEu_Simple.txt")
# write.table("Phylogeography of Pseudorasbora parva: Western Europe Origin (simple) <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_RoutesInvasionWesternEu_Simple.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
# write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_RoutesInvasionWesternEu_Simple.txt",append=TRUE)
# for (i in 1:nrow(abc.df.tmp)) {
#   write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_RoutesInvasionWesternEu_Simple.txt",append=TRUE)
# }
file.create("Data/ABC/ABC_Pparva_SourceWesternEurope.txt")
write.table("Phylogeography of Pseudorasbora parva: Western Europe Origin (simple) <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_Pparva_SourceWesternEurope.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_Pparva_SourceWesternEurope.txt",append=TRUE)
for (i in 1:nrow(abc.df.tmp)) {
  write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_Pparva_SourceWesternEurope.txt",append=TRUE)
}
rm(abc.df.tmp)



#---------------------------------
# Scenario Western Invasion without Italia
abc.df.tmp=rbind(abc.df[which(abc.df$pop %in% c("S1", "S2", "S16","S10","S11",
                                                "S13", "S14", "S15", "S17",
                                                "Aus", "Bel", "Hun", "Pol", "Spa", "UK")),]
)

# Scenario with 3 candidate native pop. (the one tested)
# Cluster sampled locations into demes
abc.df.tmp[,3]=as.character(abc.df.tmp[,3])
abc.df.tmp[which(abc.df.tmp$pop %in% c("Aus", "Bel", "Hun", "Pol", "Spa", "UK")),3]="1" # Western Europe
abc.df.tmp[which(abc.df.tmp$pop %in% c("S1", "S2", "S16")),3]="2" # North central China, without S3 (unresolved deme assignment)
abc.df.tmp[which(abc.df.tmp$pop %in% c("S13", "S14", "S15", "S17")),3]="3" # North China
abc.df.tmp[which(abc.df.tmp$pop %in% c("S10","S11")),3]="4" # Admixed continental China
# abc.df.tmp[which(abc.df.tmp$pop %in% c("S4","S6")),3]="5" # Admixed coastal China
# Reorder the dataset, DIY ABC names pop in incoming order
abc.df.tmp = abc.df.tmp[order(abc.df.tmp$pop),]
abc.df.tmp$pop
file.create("Data/ABC/ABC_Pparva_SourceWesternEurope_woIta.txt")
write.table("Phylogeography of Pseudorasbora parva: Western Europe Origin (simple, withour Italia) <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_Pparva_SourceWesternEurope.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_Pparva_SourceWesternEurope_woIta.txt",append=TRUE)
for (i in 1:nrow(abc.df.tmp)) {
  write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_Pparva_SourceWesternEurope_woIta.txt",append=TRUE)
}
rm(abc.df.tmp)

############################
# Eastern European invasion: we tested Eastern Europe source from 3 most probable origins
# Eastern Europe is highly admixed, so there is a lot of plausible origins to discriminate
# Japan is now a plausible origin for Eastern Europe
# We cannot make an economy of scenarios
# abc.df.tmp=rbind(abc.df[which(abc.df$pop %in% c("Jap", "S4", "S6", "S1", "S2", "S16","S10","S11", "S13",
#                                                 "S14", "S15", "S17", "S9", "S18", "S19", "S20",
#                                                 "Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK")),]
# )
# # Cluster sampled locations into demes
# abc.df.tmp[,3]=as.character(abc.df.tmp[,3])
# abc.df.tmp[which(abc.df.tmp$pop %in% c("Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK")),3]="1"
# abc.df.tmp[which(abc.df.tmp$pop %in% c("S4","S6")),3]="2"
# abc.df.tmp[which(abc.df.tmp$pop %in% c("S10","S11")),3]="3"
# abc.df.tmp[which(abc.df.tmp$pop %in% c("S13","S14", "S15", "S17")),3]="4"
# abc.df.tmp[which(abc.df.tmp$pop %in% c("S1", "S2", "S16")),3]="5"
# abc.df.tmp[which(abc.df.tmp$pop %in% c("S9", "S18", "S19", "S20")),3]="6"
# abc.df.tmp[which(abc.df.tmp$pop %in% c("Jap")),3]="7"
# 
# # Print into file
# file.create("Data/ABC/ABC_RoutesInvasionWesternEu.txt")
# write.table("Phylogeography of Pseudorasbora parva <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_RoutesInvasionWesternEu.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
# write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_RoutesInvasionWesternEu.txt",append=TRUE)
# for (i in 1:nrow(abc.df.tmp)) {
#   write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_RoutesInvasionWesternEu.txt",append=TRUE)
# }
# rm(abc.df.tmp)


############################
# Eastern European invasion (SIMPLE): we tested Eastern Europe source from 3 most probable source origins
abc.df.tmp=rbind(abc.df[which(abc.df$pop %in% c("S4", "S6", "S10", "S11", "S9", "S18", "S19", "S20",
                                                "Bul1", "Bul2", "Tur")),]
)
# Cluster sampled locations into demes
abc.df.tmp[,3]=as.character(abc.df.tmp[,3])
abc.df.tmp[which(abc.df.tmp$pop %in% c("Bul1", "Bul2", "Tur")),3]="1" # Eastern Europe
abc.df.tmp[which(abc.df.tmp$pop %in% c("S4","S6")),3]="2" # Admixed coastal China 
abc.df.tmp[which(abc.df.tmp$pop %in% c("S10","S11")),3]="3" # Admixed continental China
abc.df.tmp[which(abc.df.tmp$pop %in% c("S9", "S18", "S19", "S20")),3]="4" # South China
# Reorder the dataset, DIY ABC names pop in incoming order
abc.df.tmp = abc.df.tmp[order(abc.df.tmp$pop),]
abc.df.tmp$pop
# Print into files
# file.create("Data/ABC/ABC_RoutesInvasionEasternEu_Simple.txt")
# write.table("Phylogeography of Pseudorasbora parva: Eastern Europe Origin (simple) <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_RoutesInvasionEasternEu_Simple.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
# write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_RoutesInvasionEasternEu_Simple.txt",append=TRUE)
# for (i in 1:nrow(abc.df.tmp)) {
#   write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_RoutesInvasionEasternEu_Simple.txt",append=TRUE)
# }
file.create("Data/ABC/ABC_Pparva_SourceEasternEurope.txt")
write.table("Phylogeography of Pseudorasbora parva: Eastern Europe Origin (simple) <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_Pparva_SourceEasternEurope.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_Pparva_SourceEasternEurope.txt",append=TRUE)
for (i in 1:nrow(abc.df.tmp)) {
  write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_Pparva_SourceEasternEurope.txt",append=TRUE)
}
rm(abc.df.tmp)





############################
# Iranese invasion (SIMPLE): we tested Iran source from 3 most probable source origins
abc.df.tmp=rbind(abc.df[which(abc.df$pop %in% c("S4", "S6", "Jap", "S1", "S2", "S16",
                                                "Ira")),]
)
# Cluster sampled locations into demes
abc.df.tmp[,3]=as.character(abc.df.tmp[,3])
abc.df.tmp[which(abc.df.tmp$pop %in% c("Ira")),3]="1" # Iran
abc.df.tmp[which(abc.df.tmp$pop %in% c("S4","S6")),3]="2" # Admixed coastal China 
abc.df.tmp[which(abc.df.tmp$pop %in% c("Jap")),3]="3" # Japan
abc.df.tmp[which(abc.df.tmp$pop %in% c("S1", "S2", "S16")),3]="4" # North central China
# Reorder the dataset, DIY ABC names pop in incoming order
abc.df.tmp = abc.df.tmp[order(abc.df.tmp$pop),]
abc.df.tmp$pop
# Print into file
# file.create("Data/ABC/ABC_RoutesInvasionIran_Simple.txt")
# write.table("Phylogeography of Pseudorasbora parva: Iran Origin (simple) <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_RoutesInvasionIran_Simple.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
# write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_RoutesInvasionIran_Simple.txt",append=TRUE)
# for (i in 1:nrow(abc.df.tmp)) {
#   write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_RoutesInvasionIran_Simple.txt",append=TRUE)
# }
file.create("Data/ABC/ABC_Pparva_SourceIran.txt")
write.table("Phylogeography of Pseudorasbora parva: Iran Origin (simple) <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_Pparva_SourceIran.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_Pparva_SourceIran.txt",append=TRUE)
for (i in 1:nrow(abc.df.tmp)) {
  write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_Pparva_SourceIran.txt",append=TRUE)
}
rm(abc.df.tmp)

#==========================================================
# GLOBAL SCENARIOS OF INVASION
#==========================================================
# Worldwide pathways, including all invasive demes in the same scenario

# Scenarios were constructed based on 3 previous independent analysis (demes analysis)
abc.df.tmp=rbind(abc.df[which(abc.df$pop %in% c("S4", "S6", "S1", "S2", "S16","S10","S11",
                                                "S9", "S18", "S19", "S20", "Jap",
                                                "Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK",
                                                "Bul1", "Bul2", "Tur", "Ira")),]
)

# Cluster sampled locations into demes
abc.df.tmp[,3]=as.character(abc.df.tmp[,3])
abc.df.tmp[which(abc.df.tmp$pop %in% c("Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK")),3]="1" # Western Europe
abc.df.tmp[which(abc.df.tmp$pop %in% c("Bul1", "Bul2", "Tur")),3]="2" # Eastern Europe
abc.df.tmp[which(abc.df.tmp$pop %in% c("Ira")),3]="3" # Iran
abc.df.tmp[which(abc.df.tmp$pop %in% c("S1", "S2", "S16")),3]="4" # North central China
abc.df.tmp[which(abc.df.tmp$pop %in% c("S4","S6")),3]="5" # Admixed coastal China
abc.df.tmp[which(abc.df.tmp$pop %in% c("S10","S11")),3]="6" # Admixed continental China
abc.df.tmp[which(abc.df.tmp$pop %in% c("S9","S18", "S19", "S20")),3]="7" # South China
abc.df.tmp[which(abc.df.tmp$pop %in% c("Jap")),3]="8" # Japan
# Reorder the dataset, DIY ABC names pop in incoming order
abc.df.tmp = abc.df.tmp[order(abc.df.tmp$pop),]
abc.df.tmp$pop
# Print into file
# file.create("Data/ABC/ABC_RoutesInvasionIran_Simple.txt")
# write.table("Phylogeography of Pseudorasbora parva: Iran Origin (simple) <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_RoutesInvasionIran_Simple.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
# write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_RoutesInvasionIran_Simple.txt",append=TRUE)
# for (i in 1:nrow(abc.df.tmp)) {
#   write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_RoutesInvasionIran_Simple.txt",append=TRUE)
# }
file.create("Data/ABC/ABC_Pparva_Global.txt")
write.table("Phylogeography of Pseudorasbora parva: Global invasion pathways <NM=1.0NF> <MAF=hudson>","Data/ABC/ABC_Pparva_Global.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
write(paste("IND SEX POP", paste(rep("A",N.loci),collapse = " "), sep=" "),file="Data/ABC/ABC_Pparva_Global.txt",append=TRUE)
for (i in 1:nrow(abc.df.tmp)) {
  write(paste(abc.df.tmp[i,],collapse=" "),file="Data/ABC/ABC_Pparva_Global.txt",append=TRUE)
}
rm(abc.df.tmp)


#==========================================================
# THE END