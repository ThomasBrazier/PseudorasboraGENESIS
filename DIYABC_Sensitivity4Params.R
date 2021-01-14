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
# load(paste(datadir,"/Pparva.3000.Rda",sep=""))


#==========================================================
# PARAMETER ESTIMATES SENSITIVITY TESTS
#==========================================================

PrjDir = "/Pparva_Sensitivity4Params_2019_10_30-1" # 1,000,000 simulations based on a simulated dataset of 1,500 SNPs
# PrjDir = "/Pparva_Sensitivity4Params_2019_10_30-2" # 2,000,000 simulations based on a simulated dataset of 1,500 SNPs

# Observed summary statistics
obs.stats = read.table(paste(DIYABCdir, PrjDir, "/statobs.txt",sep=""), header = TRUE)
# Reference table
reftable = readRefTable(filename = paste(DIYABCdir, PrjDir, "/reftable.bin",sep=""),
                        header = paste(DIYABCdir, PrjDir, "/header.txt",sep=""), N = 0)

stats = reftable$stats # A data frame of summary statistics
params = reftable$params # A data frame of parameters drawn from prior distributions
colnames(stats) = colnames(obs.stats)

#==========================================================
# PARAMETERS' ESTIMATES WITH LOCAL LINEAR REGRESSION
#==========================================================
?abc
# Sensitivity to:
#       tolerance
#       method
#       kernel
#       sizenet
#       numnet
#       Number of simulations (default = 1,000,000)

# This function performs multivariate parameter estimation based on summary statistics using an ABC algorithm
# And return an object of class 'abc'
stats.C = as.data.frame(stats)
colnames(stats.C) = colnames(obs.stats)

#-----------------------------------------------------------
# TOLERANCE
#-----------------------------------------------------------
# Local linear regression is the faster method to compute parameters' estimates,
# hence it was used before computationnaly intensive Neural Network
# Tolerance rate allows to adjust the number of simulations retained for regression
# Only simulations closest to the dataset were retained,
# and design was to retain not too much (inaccurate estimates) or not too few (biased estimates)
list_sensitivity = c("LocLinear_tol0_1",
                     "LocLinear_tol0_05",
                     "LocLinear_tol0_01",
                     "LocLinear_tol0_005",
                     "LocLinear_tol0_001")
tol_sensitivity = c(0.1, 0.05, 0.01, 0.005, 0.001)

for (i in 1:length(list_sensitivity)) {
        test = list_sensitivity[i]
        print(test)
        if (!dir.exists(paste(DIYABCdir, PrjDir, "/", test,sep=""))) {dir.create(paste(DIYABCdir, PrjDir, "/", test,sep=""))}
        
        mod.abc=abc(target = obs.stats, param = params[,c(7,9,10,12)], sumstat = stats.C, tol = tol_sensitivity[i],
                    method = "loclinear", transf = "none", corr = TRUE, kernel = "epanechnikov")
        save(mod.abc, file = paste(DIYABCdir, PrjDir, "/", test, "/abcparam.LocLinear.Rda",sep=""))
        load(file = paste(DIYABCdir, PrjDir, "/", test,"/abcparam.LocLinear.Rda",sep=""))
        # Assessing parameters' estimates
        summary(mod.abc)
        # Checking quality of predictions
        png(paste(DIYABCdir, PrjDir, "/", test, "/abcparam.LocLinear.png",sep=""),width=500,height=500)
        plot(mod.abc, param = params[,c(7,9,10,12)], ask = FALSE)
        dev.off()
        
        #--------------------------------------------
        # POSTERIOR DISTRIBUTIONS
        # Histogram of posterior samples
        # ?hist.abc
        # hist(x=mod.abc, ask = FALSE)
        # hist(x=mod4.abc, file = paste(DIYABCdir, PrjDir, "/Figures/abcparam.LocLinear.pdf",sep=""))
        # Diagnostic plots and parameters density
        # ?plot.abc
        png(file = paste(DIYABCdir, PrjDir, "/", test, "/abcparam.LocLinear.N1b.png",sep=""),width=800,height=800)
        plot(density(mod.abc$adj.values[,2]), main = "Parameter's estimate with local linear regression",
             xlab = "Population effective size (N1b)", cex = 2, xlim = c(0, 500), col = "red") # N1b density plot of posterior with LocLinear
        lines(density(mod.abc$unadj.values[,2]), col = "black") # add N1b posterior density plot
        lines(density(params[,9]), col = "black", lty = 2) # add N1b priors density plot
        dev.off()
        # Residuals of the model lsfit()
        png(file = paste(DIYABCdir, PrjDir, "/", test, "/abcparam.res.LocLinear.N1b.png",sep=""),width=800,height=800)
        qqPlot(mod.abc$residuals[,2], main = "Residuals of local linear regression",
               xlab = "Theoretical", ylab = "Residuals", cex = 2) # N1b density plot of priors
        dev.off()
        # Density of T1, time of introduction
        png(file = paste(DIYABCdir, PrjDir, "/", test, "/abcparam.LocLinear.t1.png",sep=""),width=800,height=800)
        plot(density(mod.abc$adj.values[,1]), main = "Parameter's estimate with local linear regression",
             xlab = "Time of introduction (T1)", cex = 2, xlim = c(0, 100), col = "red") # T1 density plot of posterior with LocLinear
        lines(density(mod.abc$unadj.values[,1]), col = "black") # add T1 posterior density plot
        lines(density(params[,7]), col = "black", lty = 2) # add T1 priors density plot
        dev.off()
        # Residuals of the model lsfit()
        png(file = paste(DIYABCdir, PrjDir, "/", test, "/abcparam.res.LocLinear.t1.png",sep=""),width=800,height=800)
        qqPlot(mod.abc$residuals[,1], main = "Residuals of local linear regression",
               xlab = "Theoretical", ylab = "Residuals", cex = 2) # N1b density plot of priors
        dev.off()
        # Summaries of posterior samples by ABC
        # Calculates simple summaries of posterior samples: the minimum and maximum, the weighted mean, median, mode, and credible intervals
        # ?summary.abc
        summary(mod.abc)
}

        for (i in 1:length(list_sensitivity)) {
        test = list_sensitivity[i]
        print(test)
        if (!dir.exists(paste(DIYABCdir, PrjDir, "/", test,sep=""))) {dir.create(paste(DIYABCdir, PrjDir, "/", test,sep=""))}
        
        # Cross-validation of parameters' estimates precision
        cv.mod.abc = cv4abc(param = params[,c(7,9,10,12)], sumstat = stats.C, abc.out = NULL, nval = 1000, tols = c(tol_sensitivity[i]),
                            statistic = "median", method = "loclinear", hcorr = TRUE, transf = "none",
                            kernel = "epanechnikov")
        save(cv.mod.abc, file = paste(DIYABCdir, PrjDir, "/", test, "/cv.abcparam.LocLinear.Rda",sep=""))
        load(file = paste(DIYABCdir, PrjDir, "/", test, "/cv.abcparam.LocLinear.Rda",sep=""))
        plot(cv.mod.abc)
        png(paste(DIYABCdir, PrjDir, "/", test, "/cv.abcparam.LocLinear.png", sep=""),res=150,width=800,height=600)
        plot(cv.mod.abc)
        dev.off()
        summary(cv.mod.abc)
        # Prior error rate is simply the proportion of misclassification among all cross-validation iterations
        cv.mod.abc$true
        unlist(cv.mod.abc$estim)
        1-sum(cv.mod.abc$true == unlist(cv.mod.abc$estim))/length(cv.mod.abc$true)
}




#-----------------------------------------------------------
# KERNEL
#-----------------------------------------------------------









#==========================================================
# PARAMETERS' ESTIMATES WITH NEURAL NETWORK
#==========================================================



#-----------------------------------------------------------
# TOLERANCE
#-----------------------------------------------------------
list_sensitivity = c("NNet_tol0_5",
                     "NNet_tol0_1",
                     "NNet_tol0_05",
                     "NNet_tol0_01",
                     "NNet_tol0_005",
                     "NNet_tol0_001")
tol_sensitivity = c(0.5, 0.1, 0.05, 0.01, 0.005, 0.001)


list_sensitivity = c("NNet_tol0_005")
tol_sensitivity = c(0.005)


for (i in 1:length(list_sensitivity)) {
        test = list_sensitivity[i]
        print(test)
        if (!dir.exists(paste(DIYABCdir, PrjDir, "/", test,sep=""))) {dir.create(paste(DIYABCdir, PrjDir, "/", test,sep=""))}
        
        mod.nnet.abc=abc(target = obs.stats, param = params, sumstat = stats.C, tol = tol_sensitivity[i],
                         method = "neuralnet", transf = "none", corr = TRUE, kernel = "epanechnikov",
                         numnet = 500, sizenet = 4, lambda = c(0.0001, 0.001, 0.01), trace = TRUE,
                         maxit = 1000, MaxNWts = 10000000)
        save(mod.nnet.abc, file = paste(DIYABCdir, PrjDir, "/", test, "/abcparam.NeuralNet.Rda",sep=""))
        load(file = paste(DIYABCdir, PrjDir, "/", test, "/abcparam.NeuralNet.Rda",sep=""))
        
        #--------------------------------------------
        # POSTERIOR DISTRIBUTIONS
        # Histogram of posterior samples
        # ?hist.abc
        # hist(x=mod.nnet.abc, ask = FALSE)
        # hist(x=mod4.abc, file = paste(DIYABCdir, PrjDir, "/Figures/abcparam.LocLinear.pdf",sep=""))
        
        # Diagnostic plots
        # ?plot.abc
        # plot(mod4.abc, param = params.4[,-c(1,2,3,8,9)], ask = FALSE)
        # png(file = paste(DIYABCdir, PrjDir, "/", test, "/abcparam.NeuralNet.png",sep=""),width=800,height=800)
        # plot(density(params[,9]), main = "Parameter's estimate with neural network regression",
        #      xlab = "Population effective size (N1b)", cex = 2, ylim = c(0, 0.1)) # N1b density plot of priors
        # lines(density(mod.nnet.abc$adj.values[,9]), col = "red") # add N1b posterior density plot
        # lines(density(mod.nnet.abc$unadj.values[,9]), col = "black") # add N1b unadjusted posterior density plot
        # dev.off()
        
        # Residuals of the model lsfit()
        # png(file = paste(DIYABCdir, PrjDir, "/", test, "/abcparam.res.NeuralNet.png",sep=""),width=800,height=800)
        # qqPlot(mod.abc$residuals[,9], main = "Residuals of neural network regression",
        #        xlab = "Theoretical", ylab = "Residuals", cex = 2) # N1b density plot of priors
        # dev.off()
        
        # Summaries of posterior samples by ABC
        # Calculates simple summaries of posterior samples: the minimum and maximum, the weighted mean, median, mode, and credible intervals
        # ?summary.abc                  
        summary(mod.nnet.abc)
        rm(mod.nnet.abc)
}

for (i in 1:length(list_sensitivity)) {
        test = list_sensitivity[i]
        print(test)
        if (!dir.exists(paste(DIYABCdir, PrjDir, "/", test,sep=""))) {dir.create(paste(DIYABCdir, PrjDir, "/", test,sep=""))}
        
        # Cross-validation of parameters' estimates precision
        cv.mod.nnet.abc = cv4abc(param = params[,c(7,9)], sumstat = stats.C, abc.out = NULL, nval = 1000, tols = c(tol_sensitivity[i]),
                                 method = "neuralnet", transf = "none", corr = TRUE, kernel = "epanechnikov",
                                 numnet = 500, sizenet = 4, lambda = c(0.0001, 0.001, 0.01), trace = TRUE,
                                 maxit = 1000, MaxNWts = 100000)
        save(cv.mod.nnet.abc, file = paste(DIYABCdir, PrjDir, "/", test, "/cv.abcparam.NeuralNet.Rda",sep=""))
        load(file = paste(DIYABCdir, PrjDir, "/", test, "/cv.abcparam.NeuralNet.Rda",sep=""))
        plot(cv.mod.nnet.abc)
        png(paste(DIYABCdir, PrjDir, "/", test, "/cv.abcparam.NeuralNet.png", sep=""),res=150,width=800,height=600)
        plot(cv.mod.nnet.abc)
        dev.off()
        summary(cv.mod.nnet.abc)
        
        # Prior error rate is simply the proportion of misclassification among all cross-validation iterations
        cv.mod.nnet.abc$true
        unlist(cv.mod.nnet.abc$estim)
        1-sum(cv.mod.nnet.abc$true == unlist(cv.mod.nnet.abc$estim))/length(cv.mod.nnet.abc$true)
}


# Display summaries
for (i in 1:length(list_sensitivity)) {
        test = list_sensitivity[i]
        print(test)
        load(file = paste(DIYABCdir, PrjDir, "/", test, "/abcparam.NeuralNet.Rda",sep=""))
        summary(mod.nnet.abc)
        }


#--------------------------------------------
# POSTERIOR DISTRIBUTIONS
# Histogram of posterior samples
# ?hist.abc
# hist(x=mod.nnet.abc, ask = FALSE)
# hist(x=mod4.abc, file = paste(DIYABCdir, PrjDir, "/Figures/abcparam.LocLinear.pdf",sep=""))

# Diagnostic plots
# ?plot.abc
# plot(mod4.abc, param = params.4[,-c(1,2,3,8,9)], ask = FALSE)
# png(file = paste(DIYABCdir, PrjDir, "/", test, "/abcparam.NeuralNet.png",sep=""),width=800,height=800)
# plot(density(params[,9]), main = "Parameter's estimate with neural network regression",
#      xlab = "Population effective size (N1b)", cex = 2, ylim = c(0, 0.1)) # N1b density plot of priors
# lines(density(mod.nnet.abc$adj.values[,9]), col = "red") # add N1b posterior density plot
# lines(density(mod.nnet.abc$unadj.values[,9]), col = "black") # add N1b unadjusted posterior density plot
# dev.off()

# Residuals of the model lsfit()
# png(file = paste(DIYABCdir, PrjDir, "/", test, "/abcparam.res.NeuralNet.png",sep=""),width=800,height=800)
# qqPlot(mod.abc$residuals[,9], main = "Residuals of neural network regression",
#        xlab = "Theoretical", ylab = "Residuals", cex = 2) # N1b density plot of priors
# dev.off()

# Summaries of posterior samples by ABC
# Calculates simple summaries of posterior samples: the minimum and maximum, the weighted mean, median, mode, and credible intervals
# ?summary.abc                  
summary(mod.nnet.abc)


#-----------------------------------------------------------
# KERNEL
#-----------------------------------------------------------




#-----------------------------------------------------------
# SIZENET & NUMNET
#-----------------------------------------------------------





#==========================================================
# Figures
#==========================================================

#-----------------------------------------------------------------------------------------------------
# Mean and CI of the posterior estimate as a function of the tolerance rate
#-----------------------------------------------------------------------------------------------------
# Statistic used is the mean estimate of t1
# True values are
# Prior: mean = 120, lower = 1, upper = 200


df = read.table("Tables/Sensitivity_DIYABCParameter_t1ad.csv", header = TRUE, sep = ",")
# Transform negative values of CI to 0
df$lower[df$lower < 0] = 0
colnames(df)[1] = "Method"
plot.tol.t1ad = ggplot(data = df, aes(x = as.factor(tolerance), y = mean, colour = Method)) +
        geom_point(size = 1, position = position_dodge(width = 0.5)) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
        geom_hline(yintercept = 120, linetype="dashed") +
        geom_hline(yintercept = 200, linetype="dashed", colour = "Grey") +
        geom_hline(yintercept = 0, linetype="dashed", colour = "Grey") +
        xlab("Tolerance") + ylab("Mean posterior estimate (± C.I.)") +
        ylim(0, 220) +
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
              # legend.position = "none",
              legend.key = element_rect(fill = "white", size = 1),
              legend.text=element_text(size=18),
              legend.title=element_text(size=18))
plot.tol.t1ad

df = read.table("Tables/Sensitivity_DIYABCParameter_N1b.csv", header = TRUE, sep = ",")
# Transform negative values of CI to 0
df$lower[df$lower < 0] = 0
colnames(df)[1] = "Method"
plot.tol.N1b = ggplot(data = df, aes(x = as.factor(tolerance), y = mean, colour = Method)) +
        geom_point(size = 1, position = position_dodge(width = 0.5)) +
        geom_pointrange(aes(ymin = lower, ymax = upper), position = position_dodge(width = 0.5)) +
        geom_hline(yintercept = 1200, linetype="dashed") +
        geom_hline(yintercept = 2000, linetype="dashed", colour = "Grey") +
        geom_hline(yintercept = 0, linetype="dashed", colour = "Grey") +
        xlab("Tolerance") + ylab("Mean posterior estimate (± C.I.)") +
        ylim(0, 2400) +
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
              # legend.position = "none",
              legend.key = element_rect(fill = "white", size = 1),
              legend.text=element_text(size=18),
              legend.title=element_text(size=18))
plot.tol.N1b

ggarrange(plot.tol.t1ad, plot.tol.N1b, ncol=2, nrow=1, common.legend = TRUE, legend="bottom", labels = "auto", font.label = list(size = 20))
ggsave(paste(figuresdir,"/Sensitivity_DIYABCParameter.png",sep=""),
       device="png",dpi=320,units="cm",width=40,height=15)





#==========================================================
# END
#==========================================================