############################################################################
#   INRA - ANR project Pseudorasbora
#
#       ABC
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
# GRAPHIC LIBRARIES
library(ggplot2)
library(ggpubr)
library(reshape2)
library(RColorBrewer)
library(cowplot)
# STATS
library(MCMCglmm)
library(coda)
library(abc) # Analysing ABC simulations
library(abcrf) # Random forest on ABC simulations
library(car) # For qq-plot with confidence interval qqPlot()

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

source("Sources/cv4postpr_nnet.R")

#---------------------------------
# Import data and objects
#---------------------------------
load(paste(datadir,"/Pparva.3000.Rda",sep=""))

#==========================================================
# SIMULATIONS WITH DIY ABC
#==========================================================

# DIY ABC performs efficient and fast simulations for SNPs and complex genealogical scenarios for ABC analysis

# Confidence in priors and SuSt must be assessed before any model choice
# This step was performed in DIY ABC GUI
# Iterative process to make simulated data close to the observed data
# (i.e. observed summary statistics must be in the distribution of summary satitstics produced by simulations)
# Iterative process of adjusting priors and topologies with a low number of simulations (e.g. 10,000 per scenario)

# PATH to the project (choose your project to analyze)
path = "SourceWesternEurope"
path = "SourceWesternEurope_w/oIta"
# path = "SourceEasternEurope"
# path = "SourceIran"
# path = "Global"

if (path == "SourceWesternEurope") {PrjDir = "/Pparva_SourceWesternEurope2_2019_5_21-1"} # 400,000 simulations
if (path == "SourceWesternEurope_w/oIta") {PrjDir = "/Pparva_SourceWesternEurope_woIta_2020_1_27-1"} # 400,000 simulations
if (path == "SourceEasternEurope") {PrjDir = "/Pparva_SourceEasternEurope_2019_4_25-1"} # 400,000 simulations 
if (path == "SourceIran") {PrjDir = "/Pparva_SourceIran_2019_4_17-1"} # 400,000 simulations
if (path == "Global") {PrjDir = "/Pparva_Global_2019_4_30-2"} # 30,000 simulations

# Observed summary statistics
obs.stats = read.table(paste(DIYABCdir, PrjDir, "/statobs.txt",sep=""), header = TRUE)

# Read reftable.bin in a DIY ABC project directory
reftable = readRefTable(filename = paste(DIYABCdir, PrjDir, "/reftable.bin",sep=""),
                        header = paste(DIYABCdir, PrjDir, "/header.txt",sep=""), N = 40000) # Keep only 40,000 simulations for model selection

scenarios = reftable$scenarios # Vector of scenarios simulated
stats = reftable$stats # A data frame of summary statistics
params = reftable$params # A data frame of parameters drawn from prior distributions
colnames(stats) = colnames(obs.stats)

# 4 scenarios to test
# Split statistics and parameters' distributions between different scenarios
stats.1 = stats[which(scenarios=="1"),]
stats.2 = stats[which(scenarios=="2"),]
stats.3 = stats[which(scenarios=="3"),]
stats.4 = stats[which(scenarios=="4"),]

params.1 = params[which(scenarios=="1"),]
params.2 = params[which(scenarios=="2"),]
params.3 = params[which(scenarios=="3"),]
params.4 = params[which(scenarios=="4"),]


#---------------------------------------------------------
# Choice of priors
#---------------------------------------------------------
# Perform a priori goodness of fit using the two first components obtained with PCA
# To see if we can discriminate among the different scenarios and if simulated SuSt encompasses the observed SuSt
# If priors and topology are good, the observed SuSt must be inside the clouds of Sut
mod.gfit.pca = gfitpca(target = obs.stats, sumstat = stats,
                             index = scenarios, cprob=0.1)
png("Figures/ABC/", path,"/Fit_PCA.png",res=150,width=800,height=600)
plot(mod.gfit.pca)
dev.off()
summary(mod.gfit.pca)

##################
# Perform a test for goodness-of-fit
# The null distribution is estimated using already performed simulations contained in sumstat as pseudo-observed datasets.
# For each pseudo-observed dataset, the rejection algorithm is performed to obtain a value of the goodness-of-fit statistic

# !!!!! Warning ! This test is usually performed on lots of simulations (>> 10,000 required for random forest)
# Hence, this is just a 'rough' preliminary assessment than models are getting close to the data
mod1.gfit = gfit(target = obs.stats, sumstat = stats.1, nb.replicate = 100, tol = 0.01)
mod2.gfit = gfit(target = obs.stats, sumstat = stats.2, nb.replicate = 100, tol = 0.01)
mod3.gfit = gfit(target = obs.stats, sumstat = stats.3, nb.replicate = 100, tol = 0.01)
mod4.gfit = gfit(target = obs.stats, sumstat = stats.4, nb.replicate = 100, tol = 0.01)

# Plot
png(paste("Figures/ABC/",path,"/Goodness of fit distribution.png", sep=""),res=150,width=2400,height=2*600)
par(mfrow=c(2,2))
plot(mod1.gfit, sub = paste("Scenario 1 with p = ",summary(mod1.gfit)$pvalue,sep=""))
plot(mod2.gfit, sub = paste("Scenario 2 with p = ",summary(mod2.gfit)$pvalue,sep=""))
plot(mod3.gfit, sub = paste("Scenario 3 with p = ",summary(mod3.gfit)$pvalue,sep=""))
plot(mod4.gfit, sub = paste("Scenario 4 with p = ",summary(mod4.gfit)$pvalue,sep=""))
par(mfrow=c(1,1))
dev.off()

# Summary for p-value
summary(mod1.gfit)
summary(mod2.gfit)
summary(mod3.gfit)
summary(mod4.gfit)

# Prior calibration was performed both with DIY ABC implemented methods and 'abc' gfit
# Calibration was run iteratively until a good prior distribution was found, fitting well the data
# It was never possible to reach the median of the simulated statistics distribution with the observed data
# Single population statistics were failing more often than others, but two-pop and 3-pop statistics were pretty robust
# but, we achieved to fit significantly the data
# Once the models fit well the data, simulations can be run to test scenarios probabilities...

#==========================================================
# MODEL SELECTION WITH RANDOM FOREST
#==========================================================

# data(snp) # A simulated example in population genetics

# The 'abcrf' package allows to select a model and estimates parameters with Random Forest
# And with a low number of simulations (around 10,000 per scenario)

# Random forest don't need a lot of simulations: 10,000 simulations per scenario is enough 
data = data.frame(scenarios, stats)
# Construct a random forest from a reference table towards performing an ABC model choice
abc.mod = abcrf(scenarios~., data = data, lda = FALSE, ntree = 1000,
          sampsize=min(100000, nrow(data)), paral = TRUE, ncores = 4)
abc.mod
# The trained model was saved in a R object
save(abc.mod, file = paste(DIYABCdir, PrjDir, "/RandomForest.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/RandomForest.Rda",sep=""))

# Calculate and plot for different numbers of tree, the out-of-bag errors associated with an ABC-RF object
png(paste("Figures/ABC/",path,"/Error~ntrees.png", sep = ""),res=150,width=1000,height=1000)
err.abcrf(abc.mod, training = data, paral = TRUE, ncores = 4)
dev.off()
# Variable importance plot from a random forest: contribution of the most important variables
png(paste("Figures/ABC/",path,"/VariableImportance.png", sep = ""),res=150,width=1000,height=1000)
variableImpPlot(abc.mod)
dev.off()
# Visualize variable importance plot of a model choice ABC-RF object and the projection of the reference table on the LDA axes
# plot(abc.mod, training = data)

# Predict and evaluate the posterior probability of the MAP for new data using an ABC-RF object
abc.pred = predict(abc.mod, obs = obs.stats, training = data, paral = TRUE, ncores = 4)
save(abc.pred, file = paste(DIYABCdir, PrjDir, "/RandomForest.pred.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/RandomForest.pred.Rda",sep=""))
abc.pred

#==========================================================
# MODEL SELECTION WITH NEURAL NETWORK
#==========================================================

#----------------------------------------------------------
# POSTERIOR PROBABILITIES
# Selecting which model is the best with a Neural Network algorithm
# Tolerance must be calibrated
# If tolerance too stringent, too much rejection and
# Neural Net's performances declines quickly
# But i
stats.clean = as.matrix(stats)
colnames(stats.clean) = NULL
mod.sel = postpr(target = obs.stats, index = as.vector(scenarios),
                   sumstat = stats.clean, tol=0.05, method="neuralnet", corr = TRUE, kernel = "epanechnikov",
                   numnet = 1000, sizenet = 12, lambda = c(0.0001, 0.001, 0.01), trace = TRUE, maxit = 1000, MaxNWts = 100000)

# Model selection with Neural network was computationnally expensive
# The trained model was saved in a R object
save(mod.sel, file = paste(DIYABCdir, PrjDir, "/NeuralNet.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/NeuralNet.Rda",sep=""))

summary(mod.sel)

#----------------------------------------------------------
# Cross-validation: testing for model selection performance
# Pseudo-observed data (i.e. simulations) were reassigned to a scenario
# Accuracy is given by the proportion of re-assignment to the good scenario (under which the pseudo-observed dataset was simulated)
# Long time for computation

# Bug in the original function cv4postpr(), not passing supplementary arguments to neural network nnet()
# Modified function to pass efficiently arguments to nnet (e.g. MaxNWts)
cv.modsel = cv4postpr_nnet(index = as.vector(scenarios), sumstat = stats.clean, nval = 100, tol = 0.2, method = "neuralnet",
                        corr = TRUE, kernel = "epanechnikov", numnet = 1000, sizenet = 12,
                        lambda = c(0.0001, 0.001, 0.01), trace = TRUE, maxit = 1000, MaxNWts = 100000)

# Cross-validation with Neural network was computationnally expensive
# Results were saved in a R object
save(cv.modsel, file = paste(DIYABCdir, PrjDir, "/cv.NeuralNet.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/cv.NeuralNet.Rda",sep=""))
# ?plot.cv4abc
plot(cv.modsel)
png(paste("Figures/ABC/", path,"/Neuralnet_Confusion_matrix.png", sep=""),res=150,width=800,height=600)
plot(cv.modsel)
dev.off()
summary(cv.modsel)

# Prior error rate is simply the proportion of misclassification among all cross-validation iterations
cv.modsel$true
unlist(cv.modsel$estim)
1-sum(cv.modsel$true == unlist(cv.modsel$estim))/length(cv.modsel$true)


#----------------------------------------------------------
# RESULTS
#----------------------------------------------------------

###########################
# SOURCE OF WESTERN INVASION
###########################
# Out-of-bag prior error rate: 5.95% (i.e. re-assignment to the good scenario error rate)
# We had a good but not excellent power of random forest classification

# Confusion matrix:
#   1    2    3    4 class.error
# 1 9711  164  146   97  0.04022534
# 2  136 9050  812   34  0.09788676
# 3  102  580 9161   28  0.07192787
# 4  122   87   72 9698  0.02815913
# Most of errors came from confusion between model 2 & 3, i.e. single introduction from admixed coastal China or admixed continental China

# The Out-of-bag prior error rate as a function of number of growing trees stabilized around the mean value
# Hence, we achieved a sufficient number of growing trees for statistical precision (ntree = 1000)

# Model prediction with Random Forest: which scenario fit best the data?
# selected model votes model1 votes model2 votes model3 votes model4 post.proba
#        4           96           68          222          614        0.70255

# Model 4 was selected: single introduction from admixture in the native range between 3 most probable sources
# 614 votes on 1,000 possible votes (posterior probability 0.70)

# With the Neural Network:
# Prior error rate = 0.03

# Posterior model probabilities (neuralnet):
#   1      2      3      4 
# 0.3939 0.1566 0.1734 0.2761 
# Neural Network was contradictory with Random Forest,
# giving a posterior probability of 39% for scenario 1 (North Central China origin)
# rather than Admixed origin (27%)

# Bayes factors:
#   1      2      3      4
# 1 1.0000 2.5159 2.2713 1.4268
# 2 0.3975 1.0000 0.9028 0.5671
# 3 0.4403 1.1077 1.0000 0.6282
# 4 0.7009 1.7633 1.5919 1.0000

# Cross-validation of the Neural Network algorithm accuracy:




###########################
# SOURCE OF EASTERN INVASION
###########################
# Out-of-bag prior error rate: 5.6075%
# Low, but not perfect prior error rate. Yet, we had a good power to discriminate between scenarios

# Confusion matrix:
#   1    2    3    4 class.error
# 1 9652  157  108   90  0.03547517
# 2  147 9426  622   40  0.07904250
# 3   96  684 9073   33  0.08223751
# 4  104   82   80 9606  0.02694489
# Once again (see Western Invasion), the algorithm didn't performed well for discriminating scenarios 2 & 3
# But classification errors were not dramatic

# selected model votes model1 votes model2 votes model3 post.proba
#        1          859          124           17         0.4755167

# Scenario 4 was selected without uncertainty, with 859 votes on 1,000 and a posterior probability
# depsite a weaker posterior probability than other scenario sets (0.47)
# single origin from an admixture of the 3 probable putative sources

# With the Neural Network:
# Prior error rate = 0.04

# Posterior model probabilities (neuralnet):
#   1      2      3      4 
# 0.0220 0.0241 0.0316 0.9223 
# 
# Bayes factors:
#   1       2       3       4
# 1  1.0000  0.9142  0.6963  0.0239
# 2  1.0938  1.0000  0.7616  0.0261
# 3  1.4363  1.3130  1.0000  0.0343
# 4 41.8989 38.3042 29.1721  1.0000

# Cross-validation of the Neural Network algorithm accuracy:


###########################
# SOURCE OF IRAN INVASION
###########################
# Out-of-bag prior error rate: 6.9425%
# Prior error rate similar to other invasive demes.

# Confusion matrix:
#   1    2    3    4 class.error
# 1 9511  219   84  105  0.04113318
# 2  170 8981  907   50  0.11149584
# 3   92  842 9035   22  0.09568612
# 4  111   95   80 9696  0.02865157
# Noticeable confusion between scenarios 2 & 3 but nothing dramatic

# selected model votes model1 votes model2 votes model3 votes model4 post.proba
#       4           60          210           36             694       0.8047

# The fourth model was chosen, with 69% of votes:
# single origin from an admixture of the 3 probable putative sources

# With the Neural Network:
# Prior error rate = 0.045

#   Posterior model probabilities (neuralnet):
#   1      2      3      4 
# 0.2909 0.1671 0.1144 0.4276 
# 42% of posterior probability in favor to the model 4

# Bayes factors:
#   1      2      3      4
# 1 1.0000 1.7410 2.5414 0.6802
# 2 0.5744 1.0000 1.4598 0.3907
# 3 0.3935 0.6850 1.0000 0.2676
# 4 1.4702 2.5596 3.7364 1.0000




#----------------------------------------------------------
# SECOND STEP OF HIERARCHICAL ANALYSES
# Since there was admixture in invasive demes, what can we infer at a global scale?
# Was it a single introduction that splitted in 3 invasive demes
# or 3 independent introductions that colonised different regions?

###########################
# WORLDWIDE ROUTES OF INVASION
###########################
# Out-of-bag prior error rate: 3.4967%
# 
# Confusion matrix:
#         1    2    3  class.error
# 1 9950    2    0 0.0002009646
# 2   12 9392  618 0.0628617043
# 3    2  415 9609 0.0415918612

# selected model votes model1 votes model2 votes model3 votes model4 post.proba
#       4           60          210           36             694       0.8047

# The fourth model was chosen, with 69% of votes:
# single origin from an admixture of the 3 probable putative sources

# With the Neural Network:
# Prior error rate = 0.09333333

# Posterior model probabilities (neuralnet):
#         1      2      3 
# 0.9897 0.0059 0.0044 
# Model 1 accepted with high confidence
# 3 independent introductions from 3 independent admixture events

# Bayes factors:
#         1        2        3
# 1   1.0000 166.9438 225.4166
# 2   0.0060   1.0000   1.3503
# 3   0.0044   0.7406   1.0000



#==========================================================
# PARAMETER ESTIMATES WITH RANDOM FOREST
#==========================================================

# A reference table including 100,000 datasets is a good starting point for this step of parameter estimates (Raynal et al. 2017)
# 500 growing trees is a default choice
# Perform a random forest regression on the response variable (1 parameter to estimate)
# For the selected scenario only
# 100,000 simulations wasn't performing well for parameters' estimates
# So, re-import 1,000,000 simulations instead of 100,000 (N=0 means all simulations)

# Trial on selected scenarios at step 1: Western Europe
PrjDir = "/Pparva_SourceWesternEurope2_Params_2019_5_21-3" # Works for Western Europe scenario
# PrjDir = "/Pparva_SourceEasternEurope_Params_2019_4_25-1"
# PrjDir = "/Pparva_SourceIran_Params_2019_4_17-1" # 1,000,000 simulations

# Parameters were also estimated on the selected scenario at step 2
# PrjDir = "/Pparva_Global_Params_2019_5_6-2" # Pbm with NA values in parameters
# PrjDir = "/Pparva_Global_Params_2019_5_9-1" # Test for a solution to NA values in parameters

# Observed summary statistics
obs.stats = read.table(paste(DIYABCdir, PrjDir, "/statobs.txt",sep=""), header = TRUE)
# Reference table
reftable = readRefTable(filename = paste(DIYABCdir, PrjDir, "/reftable.bin",sep=""),
                        header = paste(DIYABCdir, PrjDir, "/header.txt",sep=""), N = 0)

stats = reftable$stats # A data frame of summary statistics
params = reftable$params # A data frame of parameters drawn from prior distributions
colnames(stats) = colnames(obs.stats)

#-----------------------------------------------------------
#   N1b, EFFECTIVE POPULATION SIZE DURING BOTTLENECK
# Here, estimate the effective population size during the bottleneck, in the scenario 4
rm(data.reg)
params = as.data.frame(params)
data.reg = data.frame(N1b = params$N1b,
                      stats)
# min.node.size = 5 for regression (?ranger)
reg.abc.N1b = regAbcrf(N1b~., data = data.reg, ntree = 1000,
                       sampsize=min(100000, nrow(data.reg)), paral = TRUE, ncores = 4,
                       min.node.size = 5)
reg.abc.N1b
save(reg.abc.N1b, file = paste(DIYABCdir, PrjDir, "/reg.abc.N1b.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/reg.abc.N1b.Rda",sep=""))

# plot.regAbcrf provides a variable importance plot used to construct the reg-ABC-RF object
plot(reg.abc.N1b)

# Based on a reg-ABC-RF object this function predicts the posterior expectation,
# median, variance, quantiles for the corresponding parameter given new dataset
pred.reg.N1b = predict(reg.abc.N1b, obs = obs.stats, training = data.reg, quantiles=c(0.025,0.975),
                       paral = TRUE, ncores = 4)
pred.reg.N1b
err.regAbcrf(reg.abc.N1b, training = data.reg,
             paral = TRUE, ncores = 4)

#-----------------------
# Plot The Posterior Density Given A New Summary Statistic
# Given a reg-ABC-RF object and a new value of the summary statistics,
# densityPlot gives the corresponding posterior density plot of the parameter,
# as well as the prior
densityPlot(reg.abc.N1b, obs = obs.stats, training = data.reg,
            paral = TRUE, ncores = 4)

#-----------------------
# EVALUATE CONFIDENCE IN PARAMETER ESTIMATES
# Predict Out-Of-Bag Posterior Expectation, Median, Variance, Quantiles And Error Measures Using A Reg-ABC-RF Object
# Based on a reg-ABC-RF object this function predicts the out-of-bag posterior expectation, median, variance, quantiles, mean squared errors, normalized mean absolute errors, credible interval coverage and relative mean range,
# for the corresponding parameter using the out-of-bag observations of the training data set.
# Mean squared errors and normalized mean absolute errors are computed both with mean and median of the response variable.
# Memory allocation issues might be encountered when the size of the training data set is large.
confidence.reg.N1b = predictOOB(reg.abc.N1b, training = data.reg,
                                paral = TRUE, ncores = 4)
confidence.reg.N1b
save(confidence.reg.N1b, file = paste(DIYABCdir, PrjDir, "/confidence.reg.N1b.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/confidence.reg.N1b.Rda",sep=""))


#==========================================================
# PARAMETERS' ESTIMATES WITH LOCAL LINEAR REGRESSION
#==========================================================
?abc
# Tolerance: 0.05 and 0.01

# This function performs multivariate parameter estimation based on summary statistics using an ABC algorithm
# And return an object of class 'abc'
stats.C = as.data.frame(stats)
colnames(stats.C) = colnames(obs.stats)

# Local linear regression is the faster method to compute parameters' estimates,
# hence it was used before computationnaly intensive Neural Network
# Tolerance rate allows to adjust the number of simulations retained for regression
# Only simulations closest to the dataset were retained,
# and design was to retain not too much (inaccurate estimates) or not too few (biased estimates)
# Starting with local linear regression with 1% tolerance (a common value in ABC parameters' esitmates)
mod.abc=abc(target = obs.stats, param = params[,c(7,9,10,12)], sumstat = stats.C, tol = 0.005,
             method = "loclinear", transf = "none", corr = TRUE, kernel = "epanechnikov")
# Local linear regression with 0.1% tolerance (1000 simulations closest to data)
# mod.abc=abc(target = obs.stats, param = params, sumstat = stats.C, tol = 0.001,
#             method = "loclinear", transf = "none", corr = TRUE, kernel = "epanechnikov")
# summary(mod.abc)
# 0.1% gave better estimates with lower CI within the wide CI at 1% but residuals wasn't gaussian
# Keeping only 1% tolerance rate for further analyses (10,000 sampled simulations)

save(mod.abc, file = paste(DIYABCdir, PrjDir, "/abcparam.LocLinear.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/abcparam.LocLinear.Rda",sep=""))
# Assessing parameters' estimates
summary(mod.abc)
# Checking quality of predictions
plot(mod.abc, param = params[,c(7,9,10,12)], ask = FALSE)

#--------------------------------------------
# POSTERIOR DISTRIBUTIONS
# Histogram of posterior samples
# ?hist.abc
hist(x=mod.abc, ask = FALSE)
# hist(x=mod4.abc, file = paste(DIYABCdir, PrjDir, "/Figures/abcparam.LocLinear.pdf",sep=""))

# Diagnostic plots and parameters density
# ?plot.abc
png(file = paste(DIYABCdir, PrjDir, "/Figures/abcparam.LocLinear.N1b.png",sep=""),width=800,height=800)
plot(density(mod.abc$adj.values[,2]), main = "Parameter's estimate with local linear regression",
     xlab = "Population effective size (N1b)", cex = 2, xlim = c(0, 500), col = "red") # N1b density plot of posterior with LocLinear
lines(density(mod.abc$unadj.values[,2]), col = "black") # add N1b posterior density plot
lines(density(params[,9]), col = "black", lty = 2) # add N1b priors density plot
dev.off()

# Residuals of the model lsfit()
png(file = paste(DIYABCdir, PrjDir, "/Figures/abcparam.res.LocLinear.N1b.png",sep=""),width=800,height=800)
qqPlot(mod.abc$residuals[,2], main = "Residuals of local linear regression",
       xlab = "Theoretical", ylab = "Residuals", cex = 2) # N1b density plot of priors
dev.off()

# Density of T1, time of introduction
png(file = paste(DIYABCdir, PrjDir, "/Figures/abcparam.LocLinear.t1.png",sep=""),width=800,height=800)
plot(density(mod.abc$adj.values[,1]), main = "Parameter's estimate with local linear regression",
     xlab = "Time of introduction (T1)", cex = 2, xlim = c(0, 100), col = "red") # T1 density plot of posterior with LocLinear
lines(density(mod.abc$unadj.values[,1]), col = "black") # add T1 posterior density plot
lines(density(params[,7]), col = "black", lty = 2) # add T1 priors density plot
dev.off()

# Residuals of the model lsfit()
png(file = paste(DIYABCdir, PrjDir, "/Figures/abcparam.res.LocLinear.t1.png",sep=""),width=800,height=800)
qqPlot(mod.abc$residuals[,1], main = "Residuals of local linear regression",
       xlab = "Theoretical", ylab = "Residuals", cex = 2) # N1b density plot of priors
dev.off()

# Residuals for 0.001 tolerance rate wasn't following normal law, strong departure from normality on edges
# Hence, going back to 0.01
# For 0.01, little departure on edges, but globally following normal law
# estimates within a restricted range
# Testing 0.05: failing to compute estimates and non gaussian residuals 
# Testing 0.005: small confidence interval and good but not perfect residuals
# After all, 0.005 gave the best results although it didn't perfectly fitted the data (departure from gaussian law in residuals)


# Summaries of posterior samples by ABC
# Calculates simple summaries of posterior samples: the minimum and maximum, the weighted mean, median, mode, and credible intervals
# ?summary.abc
summary(mod.abc)

# Cross-validation of parameters' estimates precision
cv.mod.abc = cv4abc(param = params[,c(7,9,10,12)], sumstat = stats.C, abc.out = NULL, nval = 1000, tols = c(0.01, 0.005),
                    statistic = "median", method = "loclinear", hcorr = TRUE, transf = "none",
                    kernel = "epanechnikov")
save(cv.mod.abc, file = paste(DIYABCdir, PrjDir, "/cv.abcparam.LocLinear.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/cv.abcparam.LocLinear.Rda",sep=""))
plot(cv.mod.abc)
png(paste(DIYABCdir, PrjDir, "/Figures/cv.abcparam.LocLinear.png", sep=""),res=150,width=800,height=600)
plot(cv.mod.abc)
dev.off()
summary(cv.mod.abc)

# Prior error rate is simply the proportion of misclassification among all cross-validation iterations
cv.mod.abc$true
unlist(cv.mod.abc$estim)
1-sum(cv.mod.abc$true == unlist(cv.mod.abc$estim))/length(cv.mod.abc$true)



#==========================================================
# PARAMETERS' ESTIMATES WITH NEURAL NETWORK
#==========================================================

# Given results on local linear regression, keeping only 0.1% tolerance rate
mod.nnet.abc=abc(target = obs.stats, param = params, sumstat = stats.C, tol = 0.01,
             method = "neuralnet", transf = "none", corr = TRUE, kernel = "epanechnikov",
             numnet = 500, sizenet = 4, lambda = c(0.0001, 0.001, 0.01), trace = TRUE,
             maxit = 1000, MaxNWts = 100000)
save(mod.nnet.abc, file = paste(DIYABCdir, PrjDir, "/abcparam.NeuralNet.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/abcparam.NeuralNet.Rda",sep=""))

#--------------------------------------------
# POSTERIOR DISTRIBUTIONS
# Histogram of posterior samples
# ?hist.abc
hist(x=mod.nnet.abc, ask = FALSE)
# hist(x=mod4.abc, file = paste(DIYABCdir, PrjDir, "/Figures/abcparam.LocLinear.pdf",sep=""))

# Diagnostic plots
# ?plot.abc
# plot(mod4.abc, param = params.4[,-c(1,2,3,8,9)], ask = FALSE)
png(file = paste(DIYABCdir, PrjDir, "/Figures/abcparam.NeuralNet.png",sep=""),width=800,height=800)
plot(density(params[,9]), main = "Parameter's estimate with neural network regression",
     xlab = "Population effective size (N1b)", cex = 2, ylim = c(0, 0.1)) # N1b density plot of priors
lines(density(mod.nnet.abc$adj.values[,9]), col = "red") # add N1b posterior density plot
lines(density(mod.nnet.abc$unadj.values[,9]), col = "black") # add N1b unadjusted posterior density plot
dev.off()

# Residuals of the model lsfit()
png(file = paste(DIYABCdir, PrjDir, "/Figures/abcparam.res.NeuralNet.png",sep=""),width=800,height=800)
qqPlot(mod.abc$residuals[,9], main = "Residuals of neural network regression",
       xlab = "Theoretical", ylab = "Residuals", cex = 2) # N1b density plot of priors
dev.off()

# Summaries of posterior samples by ABC
# Calculates simple summaries of posterior samples: the minimum and maximum, the weighted mean, median, mode, and credible intervals
# ?summary.abc
summary(mod.nnet.abc)

# Cross-validation of parameters' estimates precision
cv.mod.nnet.abc = cv4abc(param = params[,c(7,9)], sumstat = stats.C, abc.out = NULL, nval = 1000, tols = c(0.01, 0.005),
                         method = "neuralnet", transf = "none", corr = TRUE, kernel = "epanechnikov",
                         numnet = 500, sizenet = 4, lambda = c(0.0001, 0.001, 0.01), trace = TRUE,
                         maxit = 1000, MaxNWts = 100000)
save(cv.mod.nnet.abc, file = paste(DIYABCdir, PrjDir, "/cv.abcparam.NeuralNet.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/cv.abcparam.NeuralNet.Rda",sep=""))
plot(cv.mod.nnet.abc)
png(paste(DIYABCdir, PrjDir, "/Figures/cv.abcparam.NeuralNet.png", sep=""),res=150,width=800,height=600)
plot(cv.mod.nnet.abc)
dev.off()
summary(cv.mod.nnet.abc)

# Prior error rate is simply the proportion of misclassification among all cross-validation iterations
cv.mod.nnet.abc$true
unlist(cv.mod.nnet.abc$estim)
1-sum(cv.mod.nnet.abc$true == unlist(cv.mod.nnet.abc$estim))/length(cv.mod.nnet.abc$true)




#----------------------------------------------------------
# RESULTS
#----------------------------------------------------------

###########################
# SOURCE OF WESTERN INVASION
###########################

#----------------------------------------
# LOCAL LINEAR REGRESSION
# N1b (population effective size during bottleneck)
# Mean = 383.6812 (CI = 98.2315; 630.3740)


# t1, time of introduction
# Mean = 40.5615 (CI = 19.9269; 62.8547)


#----------------------------------------
# NEURAL NETWORK




###########################
# SOURCE OF EASTERN INVASION
###########################

#----------------------------------------
# LOCAL LINEAR REGRESSION
# N1b (population effective size during bottleneck)
# Mean = 89.7669 (CI = 32.4534; 148.2742)

# t1, time of introduction
# Mean = 38.8207 (CI = 18.2563; 61.5819)


#----------------------------------------
# NEURAL NETWORK

###########################
# SOURCE OF IRAN INVASION
###########################

#----------------------------------------
# LOCAL LINEAR REGRESSION
# N1b (population effective size during bottleneck)
# Mean = 217.2987 (CI = 185.2318; 253.4952)

# t1, time of introduction
# Mean = 65.7185 (CI = 44.3454; 89.9336)

#----------------------------------------
# NEURAL NETWORK










#----------------------------------------------------------
# FIGURE FOR REPORT
# 3 successive plots of density functions of prior and posterior
# (unadj. regression/LocLinear adj. regression) parameter N1b (effective population size during bottleneck)
load(file = paste(DIYABCdir, "/Pparva_SourceWesternEurope_2019_4_25-1", "/abcparam.LocLinear.Rda",sep=""))
reftable = readRefTable(filename = paste(DIYABCdir, PrjDir, "/reftable.bin",sep=""),
                        header = paste(DIYABCdir, PrjDir, "/header.txt",sep=""), N = 0) # Keep only 40,000 simulations for model selection
scenarios = reftable$scenarios # Vector of scenarios simulated
params = reftable$params # A data frame of parameters drawn from prior distributions
params.4 = params[which(scenarios=="4"),]
df.W = data.frame(adj = mod4.abc$adj.values[,4], unadj = mod4.abc$unadj.values[,4],
                prior = params.4[1:length(mod4.abc$adj.values[,4]),7],
                res = mod4.abc$residuals[,4],
                deme = rep("Western", length(mod4.abc$adj.values[,4])))
############
load(file = paste(DIYABCdir, "/Pparva_SourceEasternEurope_2019_4_25-1", "/abcparam.LocLinear.Rda",sep=""))
reftable = readRefTable(filename = paste(DIYABCdir, PrjDir, "/reftable.bin",sep=""),
                        header = paste(DIYABCdir, PrjDir, "/header.txt",sep=""), N = 0) # Keep only 40,000 simulations for model selection
scenarios = reftable$scenarios # Vector of scenarios simulated
params = reftable$params # A data frame of parameters drawn from prior distributions
params.4 = params[which(scenarios=="4"),]
df.E = data.frame(adj = mod4.abc$adj.values[,4], unadj = mod4.abc$unadj.values[,4],
                  prior = params.4[1:length(mod4.abc$adj.values[,4]),7],
                  res = mod4.abc$residuals[,4],
                  deme = rep("Eastern", length(mod4.abc$adj.values[,4])))
###########
load(file = paste(DIYABCdir, "/Pparva_SourceIran_2019_4_29-1", "/abcparam.LocLinear.Rda",sep=""))
reftable = readRefTable(filename = paste(DIYABCdir, PrjDir, "/reftable.bin",sep=""),
                        header = paste(DIYABCdir, PrjDir, "/header.txt",sep=""), N = 0) # Keep only 40,000 simulations for model selection
scenarios = reftable$scenarios # Vector of scenarios simulated
params = reftable$params # A data frame of parameters drawn from prior distributions
params.4 = params[which(scenarios=="4"),]
df.I = data.frame(adj = mod4.abc$adj.values[,4], unadj = mod4.abc$unadj.values[,4],
                  prior = params.4[1:length(mod4.abc$adj.values[,4]),7],
                  res = mod4.abc$residuals[,4],
                  deme = rep("Iran", length(mod4.abc$adj.values[,4])))
df = rbind(df.W, df.E, df.I)

N1b.param.plot = ggplot(data=df, aes(x = prior)) +
                geom_density(colour = "black", linetype = "dashed") +
                geom_density(aes(x = adj), colour = "red", linetype = "solid") +
                geom_density(aes(x = unadj), colour = "black", linetype = "solid") +
                xlab("Effective population size during bottleneck") + ylab("Density") +
                xlim(0, 100) +
                facet_wrap(~ deme) +
                theme(axis.line = element_line(colour = "black"),
                      panel.grid.major = element_blank(),
                      panel.grid.minor = element_blank(),
                      panel.border = element_blank(),
                      panel.background = element_blank(),
                      plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.5),
                      plot.subtitle = element_text(color="black",size=20,hjust = 0.5),
                      axis.title.x = element_text(color="black", size=20),
                      axis.title.y = element_text(color="black", size=20),
                      axis.text=element_text(size=20, colour="black"),
                      legend.key = element_rect(fill = "white", size = 1),
                      legend.text=element_text(size=20),
                      legend.title=element_text(size=20),
                      strip.text = element_text(size=20))
N1b.param.plot
ggsave(paste(figuresdir,"/ABC/Final/N1b.parameter.LocLinear.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=20)

# Plots of residuals for local linear regression
N1b.res.plot = ggplot(data=df, aes(sample = res)) +
        stat_qq() +
        stat_qq_line() +
        # geom_density(aes(x = adj), colour = "red", linetype = "solid") +
        # geom_density(aes(x = unadj), colour = "black", linetype = "solid") +
        # xlab("Effective population size during bottleneck") + ylab("Density") +
        facet_grid(deme ~ .) +
        theme(axis.line = element_line(colour = "black"),
              panel.grid.major = element_blank(),
              panel.grid.minor = element_blank(),
              panel.border = element_blank(),
              panel.background = element_blank(),
              plot.title = element_text(color="black", size=20, face="bold.italic",hjust = 0.5),
              plot.subtitle = element_text(color="black",size=20,hjust = 0.5),
              axis.title.x = element_text(color="black", size=20),
              axis.title.y = element_text(color="black", size=20),
              axis.text=element_blank(),
              axis.ticks = element_blank(),
              legend.key = element_rect(fill = "white", size = 1),
              legend.text=element_text(size=20),
              legend.title=element_text(size=20),
              strip.text = element_text(size=20))
N1b.res.plot
ggsave(paste(figuresdir,"/ABC/Final/N1b.res.LocLinear.png",sep=""),
       device="png",dpi=320,units="cm",width=20,height=60)

ggarrange(N1b.param.plot, N1b.res.plot, widths=c(6,1), labels = "auto",
          font.label = list(size = 24))
ggsave(paste(figuresdir,"/ABC/Final/N1b.LocLinear.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=20)







#==========================================================
# CROSS-VALIDATION OF POSTERIORS
#==========================================================
# 'cv4abc' performs a leave-one out cross-validation to evaluate the accuracy of parameter estimates
# and the robustness of the estimates to the tolerance rate
cv.mod4.abc.05 = cv4abc(param = params.4[,-c(1,2,3,8,9)], sumstat = stats.4.C, abc.out = NULL, nval = 100, tols = 0.05, statistic = "median",
       prior.range = NULL, method = "neuralnet", hcorr = TRUE, transf = "none",
       logit.bounds = c(0,0), subset = NULL, kernel = "epanechnikov", numnet = 500, sizenet = 4,
       lambda = c(0.0001,0.001,0.01), trace = TRUE, maxit = 1000, MaxNWts = 100000)
save(cv.mod4.abc.05, file = paste(DIYABCdir, PrjDir, "/cv05.abcparam.NeuralNet.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/cv05.abcparam.NeuralNet.Rda",sep=""))
# 2 tolerance values
cv.mod4.abc.01 = cv4abc(param = params.4[,-c(1,2,3,8,9)], sumstat = stats.4.C, abc.out = NULL, nval = 100, tols = 0.01, statistic = "median",
       prior.range = NULL, method = "neuralnet", hcorr = TRUE, transf = "none",
       logit.bounds = c(0,0), subset = NULL, kernel = "epanechnikov", numnet = 500, sizenet = 4,
       lambda = c(0.0001,0.001,0.01), trace = TRUE, maxit = 1000, MaxNWts = 100000)
save(cv.mod4.abc.01, file = paste(DIYABCdir, PrjDir, "/cv01.abcparam.NeuralNet.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/cv01.abcparam.NeuralNet.Rda",sep=""))

#--------------------------------------------
# Local regression cross-validation
cv.mod4.abc.05 = cv4abc(param = params.4[,-c(1,2,3,8,9)], sumstat = stats.4.C, abc.out = NULL, nval = 100, tols = 0.05, statistic = "mean",
                        prior.range = NULL, method = "loclinear", hcorr = TRUE, transf = "none",
                        logit.bounds = c(0,0), subset = NULL, kernel = "epanechnikov")
save(cv.mod4.abc.05, file = paste(DIYABCdir, PrjDir, "/cv05.abcparam.LocLinear.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/cv05.abcparam.LocLinear.Rda",sep=""))
# 2 tolerance values
cv.mod4.abc.01 = cv4abc(param = params.4[,-c(1,2,3,8,9)], sumstat = stats.4.C, abc.out = NULL, nval = 100, tols = 0.01, statistic = "mean",
                       prior.range = NULL, method = "loclinear", hcorr = TRUE, transf = "none",
                       logit.bounds = c(0,0), subset = NULL, kernel = "epanechnikov")
save(cv.mod4.abc.01, file = paste(DIYABCdir, PrjDir, "/cv01.abcparam.LocLinear.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/cv01.abcparam.LocLinear.Rda",sep=""))

cv.mod4.abc.05
cv.mod4.abc.01

plot(cv.mod4.abc.05, caption = "N1b")
plot(cv.mod4.abc.01, caption = "N1b")




#==========================================================
# POSTERIOR MODEL VALIDATION AND ERRORS
#==========================================================

# The quality of the selected model was evaluated  y posterior predictive check
# Simulate 10,000 times the selected scenario with priors drawn from posterior distribution
# Compare goodness-of-fit of simulations to the observed data

PrjDir = "/Pparva_SourceWesternEurope2_Params_2019_5_21-3" # Works for Western Europe scenario
# PrjDir = "/Pparva_SourceEasternEurope_Params_2019_4_25-1"
# PrjDir = "/Pparva_SourceIran_Params_2019_4_17-1" # 1,000,000 simulations

# Parameters were also estimated on the selected scenario at step 2
# PrjDir = "/Pparva_Global_Params_2019_5_6-2" # Pbm with NA values in parameters
# PrjDir = "/Pparva_Global_Params_2019_5_9-1" # Test for a solution to NA values in parameters

# Observed summary statistics
obs.stats = read.table(paste(DIYABCdir, PrjDir, "/statobs.txt",sep=""), header = TRUE)
# Reference table
reftable = readRefTable(filename = paste(DIYABCdir, PrjDir, "/reftable.bin",sep=""),
                        header = paste(DIYABCdir, PrjDir, "/header.txt",sep=""), N = 0)

stats = reftable$stats # A data frame of summary statistics
params = reftable$params # A data frame of parameters drawn from prior distributions
colnames(stats) = colnames(obs.stats)
# This function performs multivariate parameter estimation based on summary statistics using an ABC algorithm
# And return an object of class 'abc'
stats.C = as.data.frame(stats)
colnames(stats.C) = colnames(obs.stats)

# Parameters estimates with rejection algorithm: the most conservative one
mod.abc=abc(target = obs.stats, param = params, sumstat = stats.C, tol = 0.005,
            method = "rejection", transf = "none", corr = TRUE, kernel = "epanechnikov")
# Local linear regression with 0.1% tolerance (1000 simulations closest to data)
# mod.abc=abc(target = obs.stats, param = params, sumstat = stats.C, tol = 0.001,
#             method = "loclinear", transf = "none", corr = TRUE, kernel = "epanechnikov")
# summary(mod.abc)
# 0.1% gave better estimates with lower CI within the wide CI at 1% but residuals wasn't gaussian
# Keeping only 1% tolerance rate for further analyses (10,000 sampled simulations)

save(mod.abc, file = paste(DIYABCdir, PrjDir, "/abcparam.rejection.Rda",sep=""))
load(file = paste(DIYABCdir, PrjDir, "/abcparam.rejection.Rda",sep=""))
# Assessing parameters' estimates
summary(mod.abc)


##################
# Perform a priori goodness of fit using the two first components obtained with PCA
# To see if we can discriminate among the different scenarios and if simulated SuSt encompasses the observed SuSt
# If priors and topology are good, the observed SuSt must be inside the clouds of Sut
mod.gfit.pca = gfitpca(target = obs.stats, sumstat = stats,
                       index = scenarios, cprob=0.1)
png("Figures/ABC/", path,"/Posterior Model checking PCA.png",res=150,width=800,height=600)
plot(mod.gfit.pca)
dev.off()
summary(mod.gfit.pca)

##################
# Perform a test for goodness-of-fit
# The null distribution is estimated using already performed simulations contained in sumstat as pseudo-observed datasets.
# For each pseudo-observed dataset, the rejection algorithm is performed to obtain a value of the goodness-of-fit statistic
mod.gfit = gfit(target = obs.stats, sumstat = stats, nb.replicate = 100, tol = 0.01)
# Plot
png(paste("Figures/ABC/",path,"/Posterior Model checking Goodness of fit.png", sep=""),res=150,width=2400,height=2*600)
plot(mod.gfit, sub = paste("p = ",summary(mod.gfit)$pvalue,sep=""))
dev.off()






#==========================================================
# END
#==========================================================