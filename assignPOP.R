############################################################################
#   INRA - ANR project Pseudorasbora
#
#       Population assignment of invasive to source with assignPOP
#
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com
# MSc internship

# Supervisor: Scott McCairns
#             INRA UMR ESE


# DESCRIPTION
# This script performs population assignment with the R package 'assignPOP'
# https://alexkychen.github.io/assignPOP/

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
library(assignPOP) # population assignment
# GRAPHIC LIBRARIES
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
# STATS

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
# Import data and objects from Scott McCairns, INRA ESE Epix
#---------------------------------


# Pparva is a genind object that contains genotypes, with individuals and populations names
# Individuals are diploids and markers are codominant
# Genetic markers are SNPs
load(paste(datadir,"/Pparva.3000.Rda",sep=""))

#==========================================================
#  DATA
#==========================================================
# Known populations names
native.names=c("Jap", "S1", "S2", "S3", "S4", "S6","S9",
               "S10", "S11", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "Tib")
invasive.names=c("Aus", "Bel", "Bul1", "Bul2", "Hun", "Ira", "Ita", "Pol", "Spa", "Tur", "UK")

#---------------------------------
# Format in genepop file
        # Genepop format is required for multiple R packages for genetic diversity (e.g. diveRsity)
        # genind.df=genind2df(Pparva)
        # !!!!! Missing data should be indicated with 00
        # Missing data (NA) was transformed in '0000'
        # genind.df[is.na(genind.df)]="0000"
        # # ............
        # # Open connection to new genepop files: training for AssignPOP
        # file.create("AssignPOP/Data/genepop_training.txt")
        # # print first a title line
        # write.table("Population assignment training set","AssignPOP/Data/genepop_training.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
        # # Locus names: They may be given one per line, or on the same line but separated by commas.
        # write(paste(colnames(genind.df)[-1],collapse =" , "),file="AssignPOP/Data/genepop_training.txt",append=TRUE)
        # # Then, one individual per line, with the following nomenclature 'ind_name , locus1 locus2 ...'
        # # Each sample from a different geographical original is declared by a line with a POP statement.
# Some populations had very low assignment accuracy
# Policy for populations that failed at assignment tests:
# - remove S3 (low sample size, high admixture)
# Besides, we have strong priors that none invasive individuals was coming from Tibet
# Remove Tibet from putative source population
# (otherwise it will add a noisy admixed population to the assignment test)
training.pop = native.names[-c(4,18)]

        # # Print each population at a time
        # # Add native populations
        # for (i in 1:length(training.pop)) {
        #   write("POP",file="AssignPOP/Data/genepop_training.txt",append=TRUE)
        #   print(training.pop[i])
        #   subset=genind.df[which(genind.df$pop %in% training.pop[i]),]
        #   for (s in 1:nrow(subset)) {
        #     write(paste(subset[s,1],paste(subset[s,-1],collapse=" "),sep=", "),file="AssignPOP/Data/genepop_training.txt",append=TRUE)
        #   }
        # }

# # OR transform pop in putative demes
        # file.create("AssignPOP/Data/genepop_training_demes.txt")
        # write.table("Population assignment training set","AssignPOP/Data/genepop_training_demes.txt",sep=" ",quote=FALSE,row.names=F,col.names=F)
        # write(paste(colnames(genind.df)[-1],collapse =" , "),file="AssignPOP/Data/genepop_training_demes.txt",append=TRUE)
        # genind.df.demes = genind.df
        # genind.df.demes$pop = as.character(genind.df.demes$pop)
        training.demes = c("Japan", "Coastal_Admixed_China", "Continental_Admixed_China", 'North_Central_China', "North_China", "South_China")
        # for (i in 1:nrow(genind.df.demes)) {
        #   if(genind.df.demes$pop[i] %in% c("Jap")) {genind.df.demes$pop[i] = "Japan"}
        #   if(genind.df.demes$pop[i] %in% c("S4", "S6")) {genind.df.demes$pop[i] = "Coastal_Admixed_China"}
        #   if(genind.df.demes$pop[i] %in% c("S10", "S11")) {genind.df.demes$pop[i] = "Continental_Admixed_China"}
        #   if(genind.df.demes$pop[i] %in% c("S1", "S2", "S16")) {genind.df.demes$pop[i] = "North_Central_China"}
        #   if(genind.df.demes$pop[i] %in% c("S13", "S14", "S15", "S17")) {genind.df.demes$pop[i] = "North_China"}
        #   if(genind.df.demes$pop[i] %in% c("S9", "S18", "S19", "S20")) {genind.df.demes$pop[i] = "South_China"}
        # }
        # # Print each population at a time
        # # Add native demes
        # for (i in 1:length(training.demes)) {
        #   write("POP",file="AssignPOP/Data/genepop_training_demes.txt",append=TRUE)
        #   print(training.demes[i])
        #   subset=genind.df.demes[which(genind.df.demes$pop %in% training.demes[i]),]
        #   for (s in 1:nrow(subset)) {
        #     write(paste(subset[s,1],paste(subset[s,-1],collapse=" "),sep=", "),file="AssignPOP/Data/genepop_training_demes.txt",append=TRUE)
        #   }
        # }
# Import genetic data in GENEPOP format
        # genetic.training = read.Genepop("AssignPOP/Data/genepop_training.txt", pop.names=training.pop) # Known populations to train the model
# Or use demes for training (samples considered as replicates of the same genetic population)
genetic.training = read.Genepop("AssignPOP/Data/genepop_training_demes.txt", pop.names = training.demes) # Inferred putative demes in Asia to train the model



#==========================================================
#  TRAINING
#==========================================================
# https://alexkychen.github.io/assignPOP/analyze.html

# Remove loci with low variance
# genetic.training = reduce.allele(genetic.training, p = 0.95)

#==========================================================
# Perform resampling cross-validation (training of the model)
#----------------------------------------------------------
# 2 resampling methods to evaluate baseline data

# Monte-Carlo cross-validation helps estimate the mean and variance of assignment accuracy
# through resampling random training individuals

# K-fold cross-validation helps determine membership probability across all individuals
# through using one group as test individuals and the remaining K-1 groups as training individuals
# (hence every individual is tested once).

#----------------------------------------------------------
# Monte-Carlo cross-validation by Random Forest
#----------------------------------------------------------
# A total of 360 assignment tests was performed
# (3 levels of training individuals by 4 levels of training loci by 30 iterations)
# assign.MC(genetic.training, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.1, 0.25, 0.5, 1),
#           loci.sample="fst", iterations=30, model="randomForest", dir="AssignPOP/assignPOP.MC/",
#           multiprocess = TRUE) # Beware of overheating when multiprocessing !

# To improve assignment quality:
# - add more iterations
# - do not use training loci level of 0.1 (less efficient than others, high error rate)
# (3 levels of training individuals by 3 levels of training loci by 100 iterations) = 900 assignment tests
assign.MC(genetic.training, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.25, 0.5, 1),
          loci.sample="fst", iterations=100, model="randomForest", dir="AssignPOP/assignPOP.MC.randomForest/",
          multiprocess = TRUE) # Beware of overheating when multiprocessing !
# Try the Support Vector Machine (SVM) algorithm
assign.MC(genetic.training, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.25, 0.5, 1),
          loci.sample="fst", iterations=100, model="svm", dir="AssignPOP/assignPOP.MC.svm/",
          multiprocess = TRUE) # Beware of overheating when multiprocessing !
# Try the naive Bayes algorithm
assign.MC(genetic.training, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.25, 0.5, 1),
          loci.sample="fst", iterations=100, model="naiveBayes", dir="AssignPOP/assignPOP.MC.naiveBayes/",
          multiprocess = TRUE) # Beware of overheating when multiprocessing !


##########################
# Assignment on demes

# randomForest
assign.MC(genetic.training, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.25, 0.5, 1),
          loci.sample="fst", iterations=100, model="randomForest", dir="AssignPOP/assignPOP.MC.randomForest.demes/",
          multiprocess = TRUE) # Beware of overheating when multiprocessing !
# Support Vector Machine (SVM) algorithm
assign.MC(genetic.training, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.25, 0.5, 1),
          loci.sample="fst", iterations=100, model="svm", dir="AssignPOP/assignPOP.MC.svm.demes/",
          multiprocess = TRUE) # Beware of overheating when multiprocessing !
# naive Bayes algorithm
assign.MC(genetic.training, train.inds=c(0.5, 0.7, 0.9), train.loci=c(0.25, 0.5, 1),
          loci.sample="fst", iterations=100, model="naiveBayes", dir="AssignPOP/assignPOP.MC.naiveBayes.demes/",
          multiprocess = TRUE) # Beware of overheating when multiprocessing !


#----------------------------------------------------------
# K-fold cross validation
#----------------------------------------------------------
# K-fold cross-validation in which individuals from each population are divided into K groups.
# One of the K groups is tested by the predictive model built based on the remaining K-1 groups.
# Such an assignment test repeats until every group was tested.

#  K-fold cross-validation performed a total of 99 assignment tests (99 = (3+4+5+6+7+8) folds * 3 levels of training loci)
assign.kfold(genetic.training, k.fold=c(3,4,5,6,7,8), train.loci=c(0.25, 0.5, 1), 
              loci.sample="random", model="lda", dir="AssignPOP/assignPOP.Kfold/",
             multiprocess = FALSE) # Multiprocessed failed with Kfold, unable to create multiple files in the same time
# Prop of trained loci of 0.1 failed more than others (assignment accuracy<0.4)
# 0.5 and 1 performed the best (~0.75 of mean assignment accuracy, range 0.7-0.8)


# Increase the number of K fold, from 2 to high number of pop, 14 ( 2:14 * 4 = 416 assignment tests)
# assign.kfold(genetic.training, k.fold=c(2:14), train.loci=c(0.1, 0.25, 0.5, 1), 
#              loci.sample="random", model="lda", dir="AssignPOP/assignPOP.Kfold/",
#              multiprocess = FALSE) # Multiprocessed failed with Kfold, unable to create multiple files in the same time
# FAILED

# Trainings are stored in files, in their directory
# They don't have to be run again, once correctly parameterized

#==========================================================
# Visualize results
#----------------------------------------------------------

#----------------------------------------------------------
# Calculate assignment accuracy
#----------------------------------------------------------


# Assignment matrix heatmap
# Self assignment to the known population


# At this step, assignment accuracies can help to choose a training algorithm and training parameters
# to achieve the best predictive model

# When resampling cross-validations are done, calculate assignment accuracies
accu.MC.randomForest = accuracy.MC(dir = "AssignPOP/assignPOP.MC.randomForest/")
accu.MC.svm = accuracy.MC(dir = "AssignPOP/assignPOP.MC.svm/")
accu.MC.naiveBayes = accuracy.MC(dir = "AssignPOP/assignPOP.MC.naiveBayes/")



# For demes:
accu.MC.randomForest.demes = accuracy.MC(dir = "AssignPOP/assignPOP.MC.randomForest.demes/")
accu.MC.svm.demes = accuracy.MC(dir = "AssignPOP/assignPOP.MC.svm.demes/")
accu.MC.naiveBayes.demes = accuracy.MC(dir = "AssignPOP/assignPOP.MC.naiveBayes.demes/")


# K-fold cross validation procedure
accu.Kfold = accuracy.kfold(dir = "AssignPOP/assignPOP.Kfold/")
# or read directly the results if already computed
# accu.MC = read.table("AssignPOP/assignPOP.MC/Rate_of....txt", header=T)
# accu.Kfold = read.table("AssignPOP/assignPOP.Kfold/Rate_of....txt", header=T)

#-------------------------------------------
# Assignment accuracy plots

# On populations:

# png("AssignPOP/Figures/Accuracy plot.MC.randomForest.png", res = 150, width = 800, height = 600)
# accuracy.plot(accu.MC, pop = "all") + ylim(0,1)
# dev.off()
png("AssignPOP/Figures/Accuracy plot.MC.svm.png", res = 150, width = 800, height = 600)
accuracy.plot(accu.MC.svm, pop = "all") + ylim(0,1)
dev.off()
# png("AssignPOP/Figures/Accuracy plot.MC.naiveBayes.png", res = 150, width = 800, height = 600)
# accuracy.plot(accu.MC.randomForest, pop = "all") + ylim(0,1)
# dev.off()
# png("AssignPOP/Figures/Accuracy.plot.Kfold.png", res = 150, width = 800, height = 600)
# accuracy.plot(accu.Kfold, pop = "all") + ylim(0,1)
# dev.off()

# On demes:

# png("AssignPOP/Figures/Accuracy plot.MC.randomForest.demes.png", res = 150, width = 800, height = 600)
# accuracy.plot(accu.MC, pop = "all") + ylim(0,1)
# dev.off()
png("AssignPOP/Figures/Accuracy plot.MC.svm.demes.png", res = 150, width = 800, height = 600)
accuracy.plot(accu.MC.svm.demes, pop = "all") + ylim(0,1)
dev.off()
# png("AssignPOP/Figures/Accuracy plot.MC.naiveBayes.demes.png", res = 150, width = 800, height = 600)
# accuracy.plot(accu.MC.randomForest, pop = "all") + ylim(0,1)
# dev.off()
# png("AssignPOP/Figures/Accuracy.plot.Kfold.demes.png", res = 150, width = 800, height = 600)
# accuracy.plot(accu.Kfold, pop = "all") + ylim(0,1)
# dev.off()

# Demes had better performance for assignment accuracy in training set, ranging from 0.85 to 1

#------------------------------------------
# Beautiful plots for report
accuracy.plot(accu.MC.svm, pop=c("all", training.pop)) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) + #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci")+ #Add a plot title
  theme(plot.title = element_text(size=20, face="bold")) #Edit plot title text size
ggsave("AssignPOP/Figures/Accuracy.ggplot.MC.svm.png",device = "png", width = 60, height = 10, limitsize = FALSE)

accuracy.plot(accu.MC.svm.demes, pop=c("all", training.demes)) +
  ylim(0, 1) + #Set y limit between 0 and 1
  annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) + #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  ggtitle("Monte-Carlo cross-validation using genetic loci on a training set of Asian demes")+ #Add a plot title
  theme(plot.title = element_text(size=20, face="bold")) #Edit plot title text size

acc.plot = accuracy.plot(accu.MC.svm.demes, pop=c("all")) +
  ylim(0, 1) + #Set y limit between 0 and 1
  # annotate("segment",x=0.4,xend=3.6,y=0.33,yend=0.33,colour="red",size=1) + #Add a red horizontal line at y = 0.33 (null assignment rate for 3 populations)
  # ggtitle("Monte-Carlo cross-validation using genetic loci on a training set of Asian demes")+ #Add a plot title
  theme(plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=24,hjust = 0.5),
        axis.title.x = element_text(color="black", size=24),
        axis.title.y = element_text(color="black", size=24),
        axis.text=element_text(size=24, colour="black"),
        strip.text=element_text(size=22, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=20),
        legend.title=element_text(size=20))
acc.plot
ggsave("AssignPOP/Figures/Accuracy.ggplot.MC.svm.demes.png",device = "png", width = 20, height = 20, limitsize = FALSE)


#----------------------------------------------------------
# Identify informative loci
#----------------------------------------------------------
# In some cases, using a subset of high FST loci may produce similar assignment accuracy with using all available loci.
# Identification of these loci not only could help reduce time and cost in preparing samples in the future
# but also could help identify loci that might be associated with functional genes.

# This function should only be used for the cross-validation results based on
# resampling high FST training loci (loci.sample = "fst") in the assign.MC() analysis.

# The 200 best loci on 3,000
source("Sources/check.loci.R")
check.loci(dir = "AssignPOP/assignPOP.MC.svm/", top.loci = 20)
# Results are saved in a 'High_Fst_Locus_Freq.txt' file in the directory.
# top N informative loci in N rows, and each row has a list of loci sorted by its occurrence
best.loci = read.table("AssignPOP/assignPOP.MC.svm/High_Fst_Locus_Freq.txt", sep=" ",
                       header = FALSE, skip = 1, fill = TRUE)


# A large data set can be reduced to a smaller data set of the most informative loci with lots of advantages
# including reduced calculation time and memory usage



#==========================================================
#  ASSIGNMENT
#==========================================================

# Predict source of unknown individuals (invasive populations)

# test first on a subset: population assignment of the Iranese pop only
# Import data of a single pop
# genetic.assign = read.Genepop("Data/genepop_Ira.txt", pop.names="Iran") # unknown populations to assign to native populations
# genetic.assign = read.Genepop("Data/genepop_EasternEurope.txt", pop.names="EasternEurope") # unknown populations to assign to native populations
# genetic.assign = read.Genepop("Data/genepop_WesternEurope.txt") # unknown populations to assign to native populations

# All European populations at once
genetic.assign = read.Genepop("Data/genepop_invasive.txt") # unknown populations to assign to native populations


#----------------------------------------------------------
# Use training with MC
# 1.Perform assignment test using genetic data
# Results stored in a file 'AssignmentResults.txt'
assign.X(x1=genetic.training, x2=genetic.assign, dir="AssignPOP/assignPOP.MC.svm/", model="svm", mplot = TRUE)
# membership.plot(dir = "AssignPOP/assignPOP.MC.svm/")
# Results stored in a file:
res=read.table(("AssignPOP/assignPOP.MC.svm/AssignmentResult.txt"), header=TRUE)
table(res$pred.pop, res$Ind.ID) # predicted populations as a function of invasive population

# Validate the results: print mean and standard deviation across assignment tests
assign.matrix( dir="AssignPOP/assignPOP.MC.svm/", train.inds=c(0.9), train.loci=c(0.5))

# Distribution of predicted source populations
table(res$pred.pop)
plot(res$pred.pop)

# Distribution of predicted source populations for each deme
png("AssignPOP/Figures/Predicted sources per deme.png", res = 150, width = 2200, height = 800)
par(mfrow=c(1,3))
plot(res$pred.pop[which(res$Ind.ID %in% c("Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK"))], main = "Assignment test of Western Europe")
plot(res$pred.pop[which(res$Ind.ID %in% c("Bul1", "Bul2", "Tur"))], main = "Assignment test of Eastern Europe")
plot(res$pred.pop[which(res$Ind.ID %in% c("Ira"))], main = "Assignment test of Iran")
par(mfrow=c(1,1))
dev.off()

# Distribution of predictions probabilities as error assessment 
# Same plot for posterior probability as for assignment accuracy of training 
apply(res[,-c(1,2)], 1, max) # probabilities of predicted populations for each individual
png("AssignPOP/Figures/Posterior probabilities.png", res = 150, width = 800, height = 600)
hist(apply(res[,-c(1,2)], 1, max), breaks = 20, main = "Distribution of posterior probabilities", xlab = "Posterior probability", cex = 2)
dev.off()

        # # Many times, the prediction probabilities were high and in a close range for 2-4 populations that were neighbours (i.e. in the same deme/river basin)
        # # Hence, we can measure the posterior probability of population assignment to a deme rather than a pop,
        # # with the cumulative probability
        # # It means grouping posterior probabilites by demes, to see if a target population is mostly assigned to a single deme
        # res.demes=as.data.frame(cbind(as.character(res$Ind.ID),as.character(res$pred.pop)))
        # colnames(res.demes)=c("Ind.ID","pred.pop")
        # res.demes$Jap=res$Jap
        # res.demes$North_China=res$S13+res$S14+res$S15+res$S17
        # res.demes$North_Central_China=res$S1+res$S2+res$S16 # S3 discarded from North Central China
        # res.demes$Continental_Admixed_China=res$S10+res$S11
        # res.demes$Coastal_Admixed_China=res$S4+res$S6
        # res.demes$South_China=res$S9+res$S18+res$S19+res$S20
        # res.demes$pred.deme=rep(NA,nrow(res.demes))
        # # Estimating predicted deme
        #       # Predicted deme could be the deme with the highest cumulative posterior probabilities (each pop being a replicate)
        #       # for (i in 1:nrow(res.demes)) {
        #       #   res.demes$pred.deme[i]=colnames(res.demes)[which(res.demes[i,3:8]==max(res.demes[i,3:8]))+2]
        #       # }
        # # Predicted deme is the deme of the population with the highest posterior probability
        # for (i in 1:nrow(res.demes)) {
        #   if(res.demes$pred.pop[i] %in% c("Jap")) {res.demes$pred.deme[i] = "Japan"}
        #   if(res.demes$pred.pop[i] %in% c("S4", "S6")) {res.demes$pred.deme[i] = "Coastal Admixed China"}
        #   if(res.demes$pred.pop[i] %in% c("S10", "S11")) {res.demes$pred.deme[i] = "Continental\nAdmixed\nChina"}
        #   if(res.demes$pred.pop[i] %in% c("S1", "S2", "S16")) {res.demes$pred.deme[i] = "North Central China"}
        #   if(res.demes$pred.pop[i] %in% c("S13", "S14", "S15", "S17")) {res.demes$pred.deme[i] = "North China"}
        #   if(res.demes$pred.pop[i] %in% c("S9", "S18", "S19", "S20")) {res.demes$pred.deme[i] = "South China"}
        #   }
        #       # # Renaming for legend
        #       # res.demes$pred.deme[res.demes$pred.deme %in% c("Coastal_Admixed_China")] = "Coastal Admixed China"
        #       # res.demes$pred.deme[res.demes$pred.deme %in% c("Continental_Admixed_China")] = "Continental\nAdmixed\nChina"
        #       # res.demes$pred.deme[res.demes$pred.deme %in% c("Jap")] = "Japan"
        #       # res.demes$pred.deme[res.demes$pred.deme %in% c("North_Central_China")] = "North Central China"
        #       # res.demes$pred.deme[res.demes$pred.deme %in% c("South_China")] = "South China"
        # 
        # # Replace pop name by deme name, for legend
        # res.demes$Invasive.deme = as.character(res.demes$Ind.ID)
        # res.demes$Invasive.deme[res.demes$Invasive.deme %in% c("Ira")] = "Iran"
        # res.demes$Invasive.deme[res.demes$Invasive.deme %in% c("Bul1", "Bul2", "Tur")] = "Eastern Europe"
        # res.demes$Invasive.deme[res.demes$Invasive.deme %in% c("Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK")] = "Western Europe"
        # 
        # # And rerun analyses as before, on populations grouped on demes
        # table(res.demes$pred.deme)
        # plot(as.factor(res.demes$pred.deme))
        # table(res.demes$pred.deme, res.demes$Ind.ID) # predicted source deme as a function of invasive population
        # 
        # apply(res.demes[,-c(1,2,9:10)], 1, max) # probabilities of predicted populations for each individual
        # hist(apply(res.demes[,-c(1,2,9:10)], 1, max), breaks = 20, main = "Distribution of predictions probabilities to demes")
        # 
        # #---------------------------------
        # # Beautiful plot
        # cols.demes = c("#4DAF4A", "#377EB8", "#FF7F00") # Blue, Green, Orange
        # 
        # # Fill color show the target deme in the invasive range (Europe)
        # # Number of individuals is the number of individuals predicted in a given source population (x axis)
        # # Source populations are sorted by demes
        # (PredSource.plot = ggplot(data=res.demes, aes(x=pred.pop, fill=Invasive.deme))+
        #   geom_bar(width=1) +
        #   scale_fill_manual(values = cols.demes) +
        #   facet_grid(~pred.deme, scales = "free",space="free_x") +
        #   labs(x="Predicted source population", y="Number of individuals", fill="Invasive deme") +
        #   theme(axis.line = element_blank(),
        #         # axis.line.x = element_blank(), # No x axis
        #         panel.grid.major = element_blank(),
        #         panel.grid.minor = element_blank(),
        #         panel.border = element_blank(),
        #         panel.background = element_blank(),
        #         plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        #         plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
        #         axis.title.x = element_text(color="black", size=28),
        #         axis.title.y = element_text(color="black", size=28),
        #         axis.text=element_text(size=28, colour="black"),
        #         axis.text.x=element_text(size=28, colour="black"),
        #         strip.text=element_text(size=22, colour="black"),
        #         legend.key = element_rect(fill = "white", size = 1),
        #         legend.text=element_text(size=28),
        #         legend.title=element_text(size=28)))
        # ggsave("AssignPOP/Figures/Predicted sources per deme.ggplot.png",
        #        device="png",dpi=320,units="cm",width=60,height=40)
        # 

# Individuals poorly assigned, with a low probability for the predicted populations could be due to:
# - low accuracy of the assignment 
# - similar populations in the native range, and not enough polymorphism/data to discriminate among two close putative source populations
# - admixture and gene flows in native populations, reducing the assignment accuracy
# - admixture BEFORE the founder event, as multiple introductions in the same location from multiple sources is a less parcimonious hypothesis

        # # In order to get a clearer signal of population origin
        # # remove those noisy individuals with a low probability
        # # Keep only individuals with the best predictions for target populations:
        # # i.e. indviduals with more than 0.5 of prediction to the same deme
        # res.demes.trim = data.frame()
        # for (i in 1:nrow(res.demes)) {
        #   if (max(res.demes[i,-c(1,2,9:10)]) > 0.5) {
        #     res.demes.trim = rbind(res.demes.trim, res.demes[i,])
        #   }
        # }
        # # Verification
        # png("AssignPOP/Figures/Trimmed posterior probabilities.png", res = 150, width = 800, height = 600)
        # hist(apply(res.demes.trim[,-c(1,2,9:10)], 1, max), breaks = 20, main = "Distribution of posterior probabilities after trimming", xlab = "Posterior probability")
        # dev.off()
        # 
        # # Relative probabilities, as implemented in Schmidt et al. 2019
        # # probability of membership to the most likely population divided by the probability of membership to the second most likely population
        # # Individuals are considered as correctly assigned if the relative probability > 2
        # # i.e. the first probability is at least twice the next one.
        # rel.proba = c()
        # for (i in 1:nrow(res)) {
        #   rel.proba[i] = sort(res[i,-c(1,2)], decreasing = TRUE)[1]/sort(res[i,-c(1,2)], decreasing = TRUE)[2]
        # }
        # hist(as.numeric(rel.proba), breaks = 20)
        # sum(rel.proba > 2) # Only 45/168 individuals can be considered as correctly assigned (27%)
        # res.demes$rel.proba = as.numeric(rel.proba)
        # table(res.demes$pred.pop[res.demes$rel.proba>2]) # Source pop: Jap S1 S10 S11 S16 S19 S2 S4 S6 S9
        # # Mains source was S1 (14 ind.), S11 (13 ind.) and S10, S2 (5 ind.)
        # res.demes.trim = res.demes[res.demes$rel.proba>2,] # 3 invasive demes represented (e.g. 4 successes for Iran)
        # 
        # # Plotting
        # (PredSource.plot = ggplot(data=res.demes.trim, aes(x=pred.pop, fill=Invasive.deme))+
        #     geom_bar(width=1) +
        #     scale_fill_manual(values = cols.demes) +
        #     facet_grid(~pred.deme, scales = "free",space="free_x") +
        #     labs(x="Predicted source population", y="Number of individuals", fill="Invasive deme") +
        #     theme(axis.line = element_blank(),
        #           # axis.line.x = element_blank(), # No x axis
        #           panel.grid.major = element_blank(),
        #           panel.grid.minor = element_blank(),
        #           panel.border = element_blank(),
        #           panel.background = element_blank(),
        #           plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
        #           plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
        #           axis.title.x = element_text(color="black", size=28),
        #           axis.title.y = element_text(color="black", size=28),
        #           axis.text=element_text(size=28, colour="black"),
        #           axis.text.x=element_text(size=28, colour="black"),
        #           strip.text=element_text(size=22, colour="black"),
        #           legend.key = element_rect(fill = "white", size = 1),
        #           legend.text=element_text(size=28),
        #           legend.title=element_text(size=28)))
        # ggsave("AssignPOP/Figures/Predicted sources per deme.trim.ggplot.png",
        #        device="png",dpi=320,units="cm",width=60,height=40)
        # # Not very conclusive, no much information
        # 

# May be another representation could be the proportion of individuals of a given deme coming from a source population
# Standardised for uneven sampling size in invasive demes


#----------------------------
# Visualizing errors
# Boxplots of individuals prediction probabilities
# to compare with the same plot of assignment accuracy on training dataset

# Overall plot + plot per population, to assess excessive failure of some




# 2.Perform assignment test using decision tree
# assign.X(x1=genetic.training, x2=genetic.assign, dir="AssignPOP/assignPOP.MC/", model="tree", mplot = TRUE)
# 3.Perform assignment test using random forest
# assign.X(x1=genetic.training, x2=genetic.assign, dir="AssignPOP/assignPOP.MC/", model="randomForest", ntree=100, mplot = TRUE)





#----------------------------------------------------------
# Use training with K fold
# assign.X(x1=genetic.training, x2=genetic.assign, dir="AssignPOP/assignPOP.Kfold/", model="naiveBayes", mplot = TRUE)
# membership.plot(dir = "AssignPOP/assignPOP.Kfold/")
# assign.X(x1=genetic.training, x2=genetic.assign, dir="AssignPOP/assignPOP.Kfold/", model="tree", mplot = TRUE)
# assign.X(x1=genetic.training, x2=genetic.assign, dir="AssignPOP/assignPOP.Kfold/", model="randomForest", ntree=100, mplot = TRUE)
# 
# 

#==========================================================
# ASSIGNMENT TO SOURCE DEMES
#==========================================================

# Use training with MC
# 1.Perform assignment test using genetic data
# Results stored in a file 'AssignmentResults.txt'
assign.X(x1=genetic.training, x2=genetic.assign, dir="AssignPOP/assignPOP.MC.svm.demes/", model="svm", mplot = TRUE)
# membership.plot(dir = "AssignPOP/assignPOP.MC/")
# Results stored in a file:
res=read.table(("AssignPOP/assignPOP.MC.svm.demes/AssignmentResult.txt"), header=TRUE)
table(res$pred.pop, res$Ind.ID) # predicted populations as a function of invasive population

# Validate the results: print mean and standard deviation across assignment tests
      # assign.matrix( dir="AssignPOP/assignPOP.MC.svm.demes/", train.inds=c(0.9), train.loci=c(0.5))


# Distribution of predicted source populations
      # table(res$pred.pop)
      # plot(res$pred.pop)

# Distribution of predicted source populations for each deme
png("AssignPOP/Figures/Predicted sources per deme.png", res = 150, width = 2200, height = 800)
par(mfrow=c(1,3))
plot(res$pred.pop[which(res$Ind.ID %in% c("Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK"))], main = "Assignment test of Western Europe")
plot(res$pred.pop[which(res$Ind.ID %in% c("Bul1", "Bul2", "Tur"))], main = "Assignment test of Eastern Europe")
plot(res$pred.pop[which(res$Ind.ID %in% c("Ira"))], main = "Assignment test of Iran")
par(mfrow=c(1,1))
dev.off()

# Distribution of predictions probabilities as error assessment 
# Same plot for posterior probability as for assignment accuracy of training 
apply(res[,-c(1,2)], 1, max) # probabilities of predicted populations for each individual
png("AssignPOP/Figures/Posterior probabilities for demes.png", res = 150, width = 800, height = 600)
hist(apply(res[,-c(1,2)], 1, max), breaks = 20, main = "Distribution of posterior probabilities", xlab = "Posterior probability", cex = 2)
dev.off()

# Relative probabilities, as implemented in Schmidt et al. 2019
# probability of membership to the most likely population divided by the probability of membership to the second most likely population
# Individuals are considered as correctly assigned if the relative probability > 2
# i.e. the first probability is at least twice the next one.
rel.proba = c()
for (i in 1:nrow(res)) {
  rel.proba[i] = as.numeric(sort(res[i,-c(1,2)], decreasing = TRUE)[1]/sort(res[i,-c(1,2)], decreasing = TRUE)[2])
}
hist(as.numeric(rel.proba), breaks = 20)
sum(rel.proba > 2) # 100/168 individuals can be considered as correctly assigned (59%), this is much better than assignment to populations
res$rel.proba = as.numeric(rel.proba)
table(res$pred.pop[res$rel.proba>2])
# Main sources were Continental admixed China and North Central China
res.trim = res[res$rel.proba>2,] # 3 invasive demes represented (e.g. 4 successes for Iran)
write.table(res.trim, file = "AssignPOP/assignPOP.MC.svm.demes/TrimAssignmentResult.txt", quote = F, col.names = T, row.names =F)

# Distribution of predictions probabilities of trimmed individuals
# Same plot for posterior probability as for assignment accuracy of training 
apply(res.trim[,3:8], 1, max) # probabilities of predicted populations for each individual
png("AssignPOP/Figures/Posterior probabilities for demes after trimming.png", res = 150, width = 800, height = 600)
hist(apply(res.trim[,3:8], 1, max), breaks = 20, main = "Distribution of posterior probabilities", xlab = "Posterior probability", cex = 2)
dev.off()

#---------------------------------
# Beautiful plot
# cols.demes = c("#4DAF4A", "#377EB8", "#FF7F00") # Blue, Green, Orange
# cols.invasive=brewer.pal(n = 3, name = "Set2")
cols.invasive = c("#8DA0CB", "#FC8D62", "#66C2A5")

# Replace pop name by deme name, for legend
res.trim$Invasive.deme = as.character(res.trim$Ind.ID)
res.trim$Invasive.deme[res.trim$Invasive.deme %in% c("Ira")] = "Iran"
res.trim$Invasive.deme[res.trim$Invasive.deme %in% c("Bul1", "Bul2", "Tur")] = "Eastern Europe"
res.trim$Invasive.deme[res.trim$Invasive.deme %in% c("Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK")] = "Western Europe"

# Demes are displayed in this order...
res.trim$pred.pop=factor(res.trim$pred.pop, levels = c("Japan","South_China","Coastal_Admixed_China",
                                                       "Continental_Admixed_China","North_Central_China",
                                                       "North_China"))

(PredSource.plot = ggplot(data=res.trim, aes(x=pred.pop, fill=Invasive.deme))+
    geom_bar(width=0.9) +
    scale_x_discrete(labels = c("Japan" = "Japan",
                                "South_China" = "South\n China",
                                "Coastal_Admixed_China" = "Coastal\nAdmixed\nChina",
                                "Continental_Admixed_China" = "Continental\nAdmixed\nChina",
                                "North_Central_China" = "North\nCentral\nChina",
                                "North_China" = "North\nChina")) +
    scale_fill_manual(values = cols.invasive) +
    labs(x="\nPredicted source population", y="Number of individuals", fill="Invasive deme") +
    theme(axis.line = element_blank(),
          # axis.line.x = element_blank(), # No x axis
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(color="black", size=28, face="bold.italic",hjust = 0.5),
          plot.subtitle = element_text(color="black",size=28,hjust = 0.5),
          axis.title.x = element_text(color="black", size=28),
          axis.title.y = element_text(color="black", size=28),
          axis.text=element_text(size=28, colour="black"),
          axis.text.x=element_text(size=28, colour="black"),
          strip.text=element_text(size=22, colour="black"),
          axis.ticks.x = element_blank(), 
          axis.line.y = element_line(color="black",size=1,linetype="solid"),
          # axis.line.x = element_line(color="black",size=1,linetype="solid"),
          legend.key = element_rect(fill = "white", size = 1),
          legend.text=element_text(size=28),
          legend.title=element_text(size=28)))
ggsave("AssignPOP/Figures/Predicted sources per deme.trim.ggplot.png",
       device="png",dpi=320,units="cm",width=60,height=40)

##############
# Alternative representation (inset):
# Barplot of proportions of predicted source population for each invasive deme
# Sum of each line must be 1
      # df = data.frame(inv.demes = c("Western_Europe", "Eastern_Europe", "Iran"))
      # df$Japan = c(sum(res.trim$pred.pop == "Japan" & res.trim$Invasive.deme == "Western Europe")/sum(res.trim$Invasive.deme == "Western Europe"),
      #              sum(res.trim$pred.pop == "Japan" & res.trim$Invasive.deme == "Eastern Europe")/sum(res.trim$Invasive.deme == "Eastern Europe"),
      #              sum(res.trim$pred.pop == "Japan" & res.trim$Invasive.deme == "Iran")/sum(res.trim$Invasive.deme == "Iran"))
      # df$South_China = c(sum(res.trim$pred.pop == "South_China" & res.trim$Invasive.deme == "Western Europe")/sum(res.trim$Invasive.deme == "Western Europe"),
      #              sum(res.trim$pred.pop == "South_China" & res.trim$Invasive.deme == "Eastern Europe")/sum(res.trim$Invasive.deme == "Eastern Europe"),
      #              sum(res.trim$pred.pop == "South_China" & res.trim$Invasive.deme == "Iran")/sum(res.trim$Invasive.deme == "Iran"))
      # df$Coastal_Admixed_China = c(sum(res.trim$pred.pop == "Coastal_Admixed_China" & res.trim$Invasive.deme == "Western Europe")/sum(res.trim$Invasive.deme == "Western Europe"),
      #              sum(res.trim$pred.pop == "Coastal_Admixed_China" & res.trim$Invasive.deme == "Eastern Europe")/sum(res.trim$Invasive.deme == "Eastern Europe"),
      #              sum(res.trim$pred.pop == "Coastal_Admixed_China" & res.trim$Invasive.deme == "Iran")/sum(res.trim$Invasive.deme == "Iran"))
      # df$Continental_Admixed_China = c(sum(res.trim$pred.pop == "Continental_Admixed_China" & res.trim$Invasive.deme == "Western Europe")/sum(res.trim$Invasive.deme == "Western Europe"),
      #                              sum(res.trim$pred.pop == "Continental_Admixed_China" & res.trim$Invasive.deme == "Eastern Europe")/sum(res.trim$Invasive.deme == "Eastern Europe"),
      #                              sum(res.trim$pred.pop == "Continental_Admixed_China" & res.trim$Invasive.deme == "Iran")/sum(res.trim$Invasive.deme == "Iran"))
      # df$North_Central_China = c(sum(res.trim$pred.pop == "North_Central_China" & res.trim$Invasive.deme == "Western Europe")/sum(res.trim$Invasive.deme == "Western Europe"),
      #                              sum(res.trim$pred.pop == "North_Central_China" & res.trim$Invasive.deme == "Eastern Europe")/sum(res.trim$Invasive.deme == "Eastern Europe"),
      #                              sum(res.trim$pred.pop == "North_Central_China" & res.trim$Invasive.deme == "Iran")/sum(res.trim$Invasive.deme == "Iran"))
      # df$North_China = c(sum(res.trim$pred.pop == "North_China" & res.trim$Invasive.deme == "Western Europe")/sum(res.trim$Invasive.deme == "Western Europe"),
      #                              sum(res.trim$pred.pop == "North_China" & res.trim$Invasive.deme == "Eastern Europe")/sum(res.trim$Invasive.deme == "Eastern Europe"),
      #                              sum(res.trim$pred.pop == "North_China" & res.trim$Invasive.deme == "Iran")/sum(res.trim$Invasive.deme == "Iran"))
df = res.trim[,-c(2,9)]
df$Ind = as.factor(seq(1,nrow(df)))
df.melt=melt(df)
colnames(df.melt)=c("Pop", "inv.demes","Ind","source.demes","Posterior_assignment_probability")

cols.demes=c("#377EB8", "#E41A1C", "#4DAF4A", "#FF7F00", "#984EA3", "#FFFF33", "#A65628", "#F781BF",
              "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9")

# Levels in order: Iran, Eastern Europe and Western Europe
# df.melt$inv.demes = as.factor(df.melt$inv.demes)
# df.melt$inv.demes = factor(df.melt$inv.demes, levels = c("Iran", "Eastern Europe", "Western Europe"))
df.melt$Pop = as.character(df.melt$Pop)
df.melt$Pop[df.melt$Pop %in% c("Bul1", "Bul2")] = "Bul"
df.melt$Pop = as.factor(df.melt$Pop)
df.melt$Pop = factor(df.melt$Pop, levels = c("Ira", "Bul", "Tur","Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK"),
                      labels = c("Ira", "Bul", "Tur","Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK"))

(PredSource.prop.plot = ggplot(data=df.melt, aes(x=Ind, y = Posterior_assignment_probability, fill=source.demes))+
    geom_col(width=1) +
    # facet_grid(~inv.demes, scales = "free",space="free_x") +
    facet_grid(~Pop, scales = "free",space="free_x") +
    scale_fill_manual(values = cols.demes, labels = c("Japan" = "Japan",
                                                      "South_China" = "South China",
                                                      "Coastal_Admixed_China" = "Coastal Admixed China",
                                                      "Continental_Admixed_China" = "Continental Admixed China",
                                                      "North_Central_China" = "North Central China",
                                                      "North_China" = "North China")) +
    labs(x="Individuals", y="Posterior assignment probability", fill="Native deme") +
    theme(axis.line = element_blank(),
          # axis.line.x = element_blank(), # No x axis
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(color="black", size=24, face="bold.italic",hjust = 0.5),
          plot.subtitle = element_text(color="black",size=24,hjust = 0.5),
          axis.title.x = element_text(color="black", size=24),
          axis.title.y = element_text(color="black", size=24),
          axis.text=element_text(size=24, colour="black"),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(), # No x axis
          strip.text=element_text(size=22, colour="black"),
          legend.key = element_rect(fill = "white", size = 1),
          legend.text=element_text(size=20),
          legend.title=element_text(size=20)))
ggsave("AssignPOP/Figures/Predicted prop of assignment per deme.trim.ggplot.png",
       device="png",dpi=320,units="cm",width=60,height=20)

ggarrange(acc.plot, PredSource.prop.plot, widths = c(1,2), labels = "auto",
          font.label = list(size = 26))
ggsave("AssignPOP/Figures/Final.Assignment+Accuracy.png",
       device="png",dpi=320,units="cm",width=60,height=15)




#==========================================================
# END
#==========================================================