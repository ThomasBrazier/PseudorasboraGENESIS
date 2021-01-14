############################################################################
#   INRA - ANR project Pseudorasbora
#
#       Analysis of MOlecular VAriance (AMOVA)
#
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com
# MSc internship

# Supervisor: Scott McCairns
#             INRA UMR ESE

# see https://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html for tutorials


#==========================================================
# LOADING ENVIRONMENT
#==========================================================

# clear global environment: remove all variables
rm(list=ls(all=TRUE))

#----------------------
# Loading packages
library(rstudioapi)
library(ade4)
library(adegenet)
library(hierfstat)
library(poppr)



#======================================
# Loading variables & objects
#======================================

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

# Load Pparva object (8.7 Mb object)
load(paste(datadir,file="/Pparva.clean.Rda",sep=""))


#======================================
# AMOVA Analysis
#======================================
# AMOVA with adegenet()
# Based on https://popgen.nescent.org/DifferentiationSNP.html
# We performed an AMOVA to assess the validity of our clustering according to STRUCTURE and PCA clustering
?varcomp.glob
loci = Pparva@tab  # A matrix of loci

sampling = Pparva@pop
sampling
# Reclassify according to structure results
population = as.character(sampling)
population[population %in% c("UK", "Bel", "Pol", "Aus", "Hun", "Ita", "Spa")] = "Western Europe"
population[population %in% c("Bul1", "Bul2", "Tur")] = "Eastern Europe"
population[population %in% c("Ira")] = "Iran"

population[population %in% c("Jap")] = "Japan"
population[population %in% c("Tib")] = "Tibet"
population[population %in% c("S13", "S14", "S15", "S17")] = "North China"
population[population %in% c("S1", "S2", "S16")] = "North Central China"
population[population %in% c("S4", "S6")] = "Admixed Coastal China"
population[population %in% c("S10", "S11")] = "Admixed Continental China"
population[population %in% c("S9", "S18", "S19", "S20")] = "South China"
population

# Reclassify in two continents
region = as.character(sampling)
region[region %in% levels(sampling)[c(1:7,9,26, 28, 29)]] = "Europe"
region[!(region == "Europe")] = "Asia"
region

# AMOVA (Hierarchical Fst)
varcomp = varcomp.glob(levels = data.frame(region, population, sampling), loci, diploid = TRUE) # Different levels of populations (sampling location, STRUCTURE assessed populations, native/inv)

# The function test.g() tests the effect of the population on genetic differentiation.
# Individuals are randomly permuted among states. The states influence genetic differentiation at a 5% level.
test.g(loci, level = population)

# With the function test.between(), the counties are permuted among states. The states influence significantly genetic structuring.
test.between(loci, test.lev = population, rand.unit = county, nperm = 100) 

# Pairwise Fst
genet.dist(Mydata1, method = "WC84")


#======================================
# AMOVA with poppr
#======================================
# Based on http://grunwaldlab.github.io/Population_Genetics_in_R/AMOVA.html

# Set populations hierarchy
pop = data.frame(region = as.factor(region), population = as.factor(population), sampling = as.factor(sampling))
strata(Pparva) = pop
Pparva@strata
strata(Pparva)

table(strata(Pparva, ~region))
table(strata(Pparva, ~region/population, combine = FALSE))

Pparva.amova = poppr.amova(Pparva, ~region/population/sampling)
Pparva.amova
Pparva.amova$componentsofcovariance

write.table(Pparva.amova$results, sep = ",", file = "Tables/PparvaAMOVAresults.csv")
write.table(Pparva.amova$componentsofcovariance, sep = ",", file = "Tables/PparvaAMOVA.csv")
write.table(Pparva.amova$statphi, sep = ",", file = "Tables/PparvaAMOVAPhi.csv")

#======================================
# AMOVA Significance testing
#======================================
set.seed(1999)
Pparva.signif = randtest(Pparva.amova, nrepet = 999)
Pparvacc.signif <- randtest(Pparva.amovacc, nrepet = 999)

plot(Pparva.signif)
Pparva.signif
plot(Pparvacc.signif)
Pparvacc.signif

# Randomized population structure
Pparva.new = Pparva
head(strata(Pparva)[, -1])

set.seed(9001)
head(strata(Pparva)[sample(nInd(Pparva)), -1])

set.seed(9001)
strata(Pparva.new) = strata(Pparva)[sample(nInd(Pparva)), -1]
Pparva.new.amova = poppr.amova(Pparva.new, ~population/sampling)
Pparva.new.amova

Pparva.new.amova.test = randtest(Pparva.new.amova, nrepet = 999)
Pparva.new.amova.test

plot(Pparva.new.amova.test)


#======================================
# THE END
#======================================