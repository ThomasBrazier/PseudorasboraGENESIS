############################################################################
#   INRA - ANR project Pseudorasbora
#
#       Discriminant Analysis of Principal Component (DAPC)
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


# Load Pparva object (8.7 Mb object)
load(paste(datadir,"/Pparva.3000.Rda",sep=""))


#======================================
# From this point, the script is freely adapted from Scott McCairns (INRA UMR ESE)
#======================================


#########################################################################################################################################################
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
###===================================================================================================================================================###
###							Clustering of Native Populations
###===================================================================================================================================================###
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
#########################################################################################################################################################
x<-c("Jap", "S1", "S2", "S3", "S4", "S6", "S9", "S10", "S11", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "Tib")
i=1
y<-which(Pparva$pop==x[i])
# create an index y of all individuals of Asian populations
for (i in 2:length(x)){ y<-c(y, which(Pparva$pop==x[i])) }
rm(i)
rm(x)

native<-Pparva
# Retain only individuals from Asian populations
native$tab<-native$tab[y,]
native$pop<-native$pop[y]
native$ploidy<-native$ploidy[y]
rm(y)

save(native,file="Data/DAPC/native.Rda")
load(file="Data/DAPC/native.Rda")


###===================================================================================================================================================###
###  Preliminary Clustering of Data -- k Means
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# run interactively to select the number of clusters (based on minimizing BIC) and to estimate variance explained
# number of retained PCs will be set to maximum allowed for model "stability" (< N/3; 100< 300/3)
# 100 retained PCs appears to exlain 75-80% of cumulative variance
# BIC minimum occurs at k=7
#!# HOWEVER, true inflection point may be between 6 and 8 based on graphics

k.native<-find.clusters(native, max.n=19)

# Compare if clusters corresponded to observed groups (i.e. sampled locations)
table.value(table(pop(native), k.native$grp), col.lab=paste("inf", 1:19),
            row.lab=paste("ori", 1:19))
# Observed groups globally corresponded to inferred groups

###---------------------------------------------------------------------------------------------------------------------------------------------------###
###  Initial Discriminant Analysis of Principal Components (DAPC)
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# again set to retain maximum numbers of PCs & DA's
# again run interactively to examine scree-plot of F-statistics

dapc.native<-dapc(native, pop=k.native$grp, n.pca=100,n.da=7)

###---------------------------------------------------------------------------------------------------------------------------------------------------###
compoplot(dapc.native, posi="top", col=c(rainbow(6,s=1,v=1), rainbow(6,s=0.5,v=0.75)))
dapc.native$grp
compoplot(dapc.native, posi="top", col=c(rainbow(6,s=1,v=1), rainbow(6,s=0.5,v=0.75)),show.lab=TRUE,lab=dapc.native$grp)
# no individual show evidence of admixture
# Japanese indiv. globally clustered together
# Tibetan and some S10-S11 indiv. clustered together
# Clusters seemed to correspond to existing sampled populations


#compoplot(dapc.native, posi="top", col=c(rainbow(6,s=1,v=1), rainbow(6,s=0.5,v=0.75)), lab="", only.grp="1") # Currently not working due to unresolved error
# cluster specific plotting to quickly assess group characteristics


###===================================================================================================================================================###
###  Alternate Strategy for Clustering of Data:  No. PCs Retained to Achieve Set Cumulative Variance Explained
###---------------------------------------------------------------------------------------------------------------------------------------------------###
#!# am I over-fitting?
# perhaps retain fewer PCs initially
# use cumulative variance as a threshold
###---------------------------------------------------------------------------------------------------------------------------------------------------###
### 75% of variance
###---------------------------------------------------------------------------------------------------------------------------------------------------###
k.native<-find.clusters(native, max.n=19, pca.select="percVar", perc.pca=75)

# run interactively to select the number of clusters
# BIC minimum occurs at k=9
#!# HOWEVER, suspect this is also an over-fitting given previous attempts, plus:
# 1 "uptick" at 10
# near zero slope between 6 & 8

###---------------------------------------------------------------------------------------------------------------------------------------------------###
### 70% of variance
###---------------------------------------------------------------------------------------------------------------------------------------------------###
k.native<-find.clusters(native, max.n=19, pca.select="percVar", perc.pca=70)

# run interactively to select the number of clusters
# BIC minimum occurs at k=8
# 1 "downtick" at K=8

# checking k=8
# same pattern as before re: scree plot of F-stats
dapc.native<-dapc(native, pop=k.native$grp, pca.select="percVar", perc.pca=70, n.da=7)

compoplot(dapc.native, posi="top", col=c(rainbow(6,s=1,v=1), rainbow(6,s=0.5,v=0.75)))

#!#!# NB! there is NO consistency between clusters if multiples iterations are compared
# how then does one evaluate which model is truly the best representation of the data?!?
# perhaps I need to do so iteratively

###===================================================================================================================================================###
###  Alternate Strategy for Clustering of Data:  Using a "best fit" Selection Criterion (in addition to cum. var.)
###---------------------------------------------------------------------------------------------------------------------------------------------------###
test1.k<-list()
for (i in 1:100){ print(paste("Iteration No.", i, "of 100"))
  k.native<-find.clusters(native, max.n=19, pca.select="percVar", perc.pca=70, choose.n.clust=F, criterion="goodfit", n.iter=1000, n.start=sample(10:21,1))
  test1.k[[i]]<-k.native
  rm(k.native) }
rm(i)

i=1
x<-length(test1.k[[i]]$size)
for(i in 2:100){ x<-c(x, length(test1.k[[i]]$size))}
rm(i)
table(x)
hist(x)
# 100 iterations is enough
# Clear score of 97% for K=5
# BUT "Goodfit" criterion tended to favor a small K, that was not informative in our case
# Lack of resolution in native populations

# Same iterative test but with the criterion "min" of BIC
# As BIC~Nb of clusters decreases up to a plateau and stop decreasing,
# in an exponential negative shape like an island model (see details of ?find.clusters())
test1b.k<-list()
for (i in 1:100){ print(paste("Iteration No.", i, "of 100"))
  k.native<-find.clusters(native, max.n=19, pca.select="percVar", perc.pca=70, choose.n.clust=F, criterion="min", n.iter=1000, n.start=sample(10:21,1))
  test1b.k[[i]]<-k.native
  rm(k.native) }
rm(i)

i=1
x<-length(test1b.k[[i]]$size)
for(i in 2:100){ x<-c(x, length(test1b.k[[i]]$size))}
rm(i)
table(x)
hist(x)
# "min" criterion tended to over-estimate K

# Same iterative test but with the criterion "diffNgroup" of BIC
test1.k<-list()
for (i in 1:100){ print(paste("Iteration No.", i, "of 100"))
  k.native<-find.clusters(native, max.n=19, pca.select="percVar", perc.pca=70, choose.n.clust=F, criterion="diffNgroup", n.iter=1000, n.start=sample(10:21,1))
  test1.k[[i]]<-k.native
  rm(k.native) }
rm(i)

i=1
x<-length(test1.k[[i]]$size)
for(i in 2:100){ x<-c(x, length(test1.k[[i]]$size))}
rm(i)
table(x)
hist(x)
# This criterion gave K=4 as the most probable K (100/100 iterations)

###---------------------------------------------------------------------------------------------------------------------------------------------------###

# probably best to run at least 1,000 total iterations though if I want to be somewhat "Structure-esque"
###---------------------------------------------------------------------------------------------------------------------------------------------------###
for (i in 1:1000){ print(paste("Iteration No.", i, "of 1000"))
  k.native<-find.clusters(native, max.n=19, pca.select="percVar", perc.pca=70, choose.n.clust=F, criterion="goodfit", n.iter=1000, n.start=sample(10:21,1))
  test1.k[[i]]<-k.native
  rm(k.native) }
rm(i)

i=1
x<-length(test1.k[[i]]$size)
for(i in 1:1000){ x<-c(x, length(test1.k[[i]]$size))}
rm(i)
table(x)
hist(x)
# 5 clusters (n=958/1000) 

for (i in 1:1000){ print(paste("Iteration No.", i, "of 1000"))
  k.native<-find.clusters(native, max.n=19, pca.select="percVar", perc.pca=70, choose.n.clust=F, criterion="min", n.iter=1000, n.start=sample(10:21,1))
  test1b.k[[i]]<-k.native
  rm(k.native) }
rm(i)

i=1
x<-length(test1b.k[[i]]$size)
for(i in 1:1000){ x<-c(x, length(test1b.k[[i]]$size))}
rm(i)
table(x)
hist(x)
# K=9-10 were favored for minimum criterion (n=303 for K=9 & n=390 for K=10 ; n=251 for K=11) 




# BUT, perhaps should test for over-fitting

###---------------------------------------------------------------------------------------------------------------------------------------------------###
###  Step 2 in Strategy:  Resolving the Over-fitting Problem
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# cross-validation to estimate the optimal number of retained PCs
# given size of dataset, consider 100 as an upper limit
# maximum "allowed" for model stability (< N/3; 100< 259/3)
mat<-tab(native, NA.method="mean")

# 80 PCs to retain to discriminate between the 19 sampling sites, although it seemed that at least 20 is sufficient
xval<-xvalDapc(mat, grp=native@pop, n.pca=seq(5,85,5))

# 20 was enough to discriminate amongst k=5 clusters
xval<-xvalDapc(mat, grp=test1.k[[10]]$grp, n.pca=seq(5,85,5))

#!# NB
# NONE of the plots displayed a nice peaked pattern
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# a-score optimization
# Linear decline of a-score from 0 to 85
# suggests far fewer (<10) PCs to avoid overfitting
dapc.1<-dapc(native, pop=test1.k[[10]]$grp, pca.select="percVar", perc.pca=70, n.da=6)
scatter(dapc.1)
optim.a.score(dapc.1, n.pca=seq(5,85,5))


dapc.2<-dapc(native, pop=test1.k[[10]]$grp, pca.select="percVar", n.pca=85, n.da=6)
scatter(dapc.2)
dapc.3<-dapc(native, pop=test1.k[[10]]$grp, pca.select="percVar", n.pca=10, n.da=6)
scatter(dapc.3)
# scatterplot shows greater degree of separation with 85 PCs than with 10 PCs


###===================================================================================================================================================###
###  Starting Again....using what we've learned so far
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# k.native<-find.clusters(native, max.n=19, n.pca=100, choose.n.clust=F, criterion="goodfit", n.iter=1000, n.start=sample(10:19,1))
# k.native$size
# 
# xval<-xvalDapc(mat, grp=k.native$grp, n.pca=seq(5,100,5))
# 
# dapc.native<-dapc(native, pop=k.native$grp, n.pca=40, n.da=9)
# scatter(dapc.native)
# 
# compoplot(dapc.native, posi="top", col=c(rainbow(6,s=1,v=1), rainbow(6,s=0.5,v=0.75)))
# 
# rm(k.native)
# rm(dapc.native)
# rm(xval)

# preliminary test seemed sensible
# must now perform multiple iterations
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# test2.k<-list()
# for (i in 1:1000){ print(paste("Iteration No.", i, "of 1000"))
#   k.native<-find.clusters(native, max.n=21, n.pca=110, choose.n.clust=F, criterion="goodfit", n.iter=1000, n.start=sample(5:16,1))
#   test2.k[[i]]<-k.native
#   rm(k.native) }
# rm(i)
# 
# i=1
# k_freqs2<-length(test2.k[[i]]$size)
# for(i in 2:1000){ k_freqs2<-c(k_freqs2, length(test2.k[[i]]$size))}
# rm(i)
# table(k_freqs2)
# hist(k_freqs2, las=1)
# 
# i=1
# x<-test2.k[[i]]$stat
# for(i in 2:1000){ x<-c(x, test2.k[[i]]$stat)}
# rm(i)
# 
# length(which(names(sort(x))=="K=4"))
# length(which(names(sort(x))=="K=5"))
# length(which(names(sort(x))=="K=6"))
# # next smallest was k=6 with 1 iteration only
# 
# # K=5 has been clearly assessed
# # K=6 was rare but BIC was very close to K=5
# # Now search for the clusters assignation with K=5 and the lowest BIC (best working clusters for K=5)
# which(x==min(x)) # Identify the clusters with the lowest BIC for K=5
# test2.k[[1]]$stat # lowest BIC was here 890.47
# 
# which(names(sort(x))=="K=4")
# which(names(sort(x))=="K=5")

###===================================================================================================================================================###
###  Working Clusters
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# K=5 was a parcimonious number of K and the lowest BIC value
# BIC for K=6 was very close to the BIC value of K=5, even if rare, and clusters were more consistent with geography
# But clusters with K=5 did not gave us useful information about source populations
# We looked for more clusters, for a finer structuring resolution
# First "naive" graphical approach tended to favor K=8-9; Good, as we searched for a fine-scale genetic structure of source populations
# Selection by the "min" criterion favored K=9-10, but tended to produce overclustering with frequent admixture in a large pan-asian region
# So we chose to use K=7, producing more clusters and more consistent ones with the geography (Tib, Jap)
nclust=9
k.native<-find.clusters(native, max.n=19,n.pca=100,n.clust=nclust)
k.native$size
# save(k.native,file="Data/DAPC/k.native.rda")

xval<-xvalDapc(mat, grp=k.native$grp, n.pca=seq(5,100,5))

dapc.native<-dapc(native, pop=k.native$grp, n.pca=40,n.da=(nclust-1))
# save(dapc.native,file="Data/DAPC/dapc.native.rda")
# load("Data/DAPC/dapc.native.rda")

scatter(dapc.native)

compoplot(dapc.native, posi="top", col=c(rainbow(6,s=1,v=1), rainbow(6,s=0.5,v=0.75)))

# x=1
# y=2
# scatter(dapc.native, xax=x, yax=y, sub=paste("PC", y, " v. PC", x, sep=""), csub = 2, 
#         col=c(rainbow(6,s=1,v=1), rainbow(3,s=0.5,v=0.75)), posi.da="bottomright")
# rm(x)
# rm(y)


# plotting somewhat geographically-close clusters together
# layout used to better visualize larger clusters
# layout(rbind(c(1, 2), c(3, 3)))
# for(i in c(1,9,3)){ compoplot(dapc.native, legend=F, col=c(rainbow(6,s=1,v=1), rainbow(3,s=0.5,v=0.75)), only.grp=as.character(i), main=paste("Cluster", i)) }
# rm(i)
# par(mfrow=c(1,1))
# 
# 
# layout(rbind(c(1, 2), c(3, 3)))
# for(i in c(7,2,6)){ compoplot(dapc.native, legend=F, col=c(rainbow(6,s=1,v=1), rainbow(3,s=0.5,v=0.75)), only.grp=as.character(i), main=paste("Cluster", i)) }
# rm(i)
# par(mfrow=c(1,1))
# 
# 
# layout(rbind(c(1, 2), c(3, 3)))
# for(i in c(8,4,5)){ compoplot(dapc.native, legend=F, col=c(rainbow(6,s=1,v=1), rainbow(3,s=0.5,v=0.75)), only.grp=as.character(i), main=paste("Cluster", i)) }
# rm(i)
# par(mfrow=c(1,1))

###---------------------------------------------------------------------------------------------------------------------------------------------------###
### Study of admixture between inferred clusters
###---------------------------------------------------------------------------------------------------------------------------------------------------###

### Distribution of posterior probabilities across clusters
hist(dapc.native$posterior,breaks=50)
# Posterior probabilities distribution show two extreme values, close to 0 and 1
# Suggesting no admixture between clusters

assignplot(dapc.native, subset=1:300)
# No evidence of admixed individuals with DAPC


###---------------------------------------------------------------------------------------------------------------------------------------------------###
### Identifying Members of Each Cluster
###---------------------------------------------------------------------------------------------------------------------------------------------------###
native.cluster.df<-data.frame(Cluster=factor(), Ind.Site=factor(), Ind=factor())

for(g in 1:7){ x<-sort(names(k.native$grp)[which(k.native$grp==g)])

for(i in 1:length(x)){ temp.df<-data.frame(Cluster=paste("Cluster", g, sep="_"), Ind.Site=strsplit(x[i], "_")[[1]][1], Ind=x[i])

native.cluster.df<-rbind(native.cluster.df, temp.df)

rm(i)
rm(temp.df)}
rm(x) }
rm(g)

###---------------------------------------------------------------------------------------------------------------------------------------------------###
# Create a data frame with all populations allocated to each cluster
unique(native.cluster.df[,1:2])

write.csv(unique(native.cluster.df[,1:2]), "Tables/DAPC/native_cluster.csv", row.names=F,quote=F)

###---------------------------------------------------------------------------------------------------------------------------------------------------###
# To look for a given cluster

subset((unique(native.cluster.df[,1:2])), Cluster=="Cluster_3")
barplot(table(subset(native.cluster.df, Cluster=="Cluster_3")$Ind.Site), las=1)

barplot(sort(table(subset(native.cluster.df, Cluster=="Cluster_3")$Ind.Site)), las=1)

###---------------------------------------------------------------------------------------------------------------------------------------------------###
# Ou bien pour identifier chaque cluster ou se trouvent des echantillons d'une site en particuliere
subset((unique(native.cluster.df[,1:2])), Ind.Site=="Jap")

x<-c("Jap", "S1", "S2", "S3", "S4", "S6", "S9", "S10", "S11", "S13", "S14", "S15", "S16", "S17", "S18", "S19", "S20", "Tib")
for (i in 1:length(x)){ print(i)
  print(subset((unique(native.cluster.df[,1:2])), Ind.Site==x[i])) }
rm(i)
rm(x)

#========================================================================================================================================================

###===================================================================================================================================================###
###  Using what we've learned from STRUCTURE (K=6), we investigated furthermore K=6 with DAPC...
###---------------------------------------------------------------------------------------------------------------------------------------------------###

# Running 1,000 models of DAPC with the 'a priori' number of clusters (estimated from STRUCTURE): K=6
test1c.k<-list()
for (i in 1:1000){ print(paste("Iteration No.", i, "of 1000"))
  k.native<-find.clusters(native, max.n=19, pca.select="percVar", perc.pca=70, choose.n.clust=F, criterion="goodfit", n.iter=1000, n.start=sample(10:21,1))
  test1c.k[[i]]<-k.native
  rm(k.native) }
rm(i)

# Selection of the most probable K
i=1
x<-length(test1c.k[[i]]$size)
for(i in 1:1000){ x<-c(x, length(test1c.k[[i]]$size))}
rm(i)
table(x)
hist(x)

# Selection of the best model among 1,000 models
i=1
x=test1c.k[[i]]$stat
for(i in 2:1000){ x=c(x, test1c.k[[i]]$stat)}
rm(i)
table(x)
hist(x)

# Models with K=5
index.K5=(which(names(x)=="K=5"))
test1c.k[[index.K5]]
# Get a vector containing the BIC values for K=5
x=c()
for(i in index.K5){ x=c(x, test1c.k[[i]]$stat)}
x
# And then select the best model (lowest BIC value)
which(x==min(x)) # Identify the clusters with the lowest BIC for K=6 -> 2 models: 
# Models with the lowest BIC value for K=5 were explored
# and model 383 was one of the best explicative model

test1c.k[[383]]$stat # lowest BIC for K=6 was here 832.6823
# Model 383 was our chosen model for DAPC in native populations assuming K=5 as the most probable number of genetic clusters
# Don't forget to save our 1,000 iterations for best model research
save(test1c.k,file="Data/DAPC/test1c.k")
load(file="Data/DAPC/test1c.k")

# BUT, perhaps should test for over-fitting

###---------------------------------------------------------------------------------------------------------------------------------------------------###
###  Step 2 in Strategy:  Resolving the Over-fitting Problem
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# cross-validation to estimate the optimal number of retained PCs
# given size of dataset, consider 100 as an upper limit
# maximum "allowed" for model stability (< N/3; 100< 300/3)
mat=tab(native, NA.method="mean")

# 80 PCs to retain to discriminate between the 19 sampling sites, although it seemed that at least 20 is sufficient
xval=xvalDapc(mat, grp=native@pop, n.pca=seq(5,100,5))
# Flat pattern of cross-validation with a slight maximum around 40 PCAs

# 20 was enough to discriminate amongst k=6 clusters for the chosen model n?57, but no evidence of overfitting
xval=xvalDapc(mat, grp=test1c.k[[383]]$grp, n.pca=seq(5,100,5))

#!# NB
# NONE of the plots displayed a nice peaked pattern
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# a-score optimization
# Linear decline of a-score from 0 to 100 (300/3)
# suggests far fewer (<10) PCs to avoid overfitting
dapc.1<-dapc(native, pop=test1c.k[[383]]$grp, pca.select="percVar", perc.pca=70, n.da=5)
scatter(dapc.1)
optim.a.score(dapc.1, n.pca=seq(5,100,5))


dapc.2<-dapc(native, pop=test1c.k[[383]]$grp, pca.select="percVar", n.pca=85, n.da=5)
scatter(dapc.2)
dapc.3<-dapc(native, pop=test1c.k[[383]]$grp, pca.select="percVar", n.pca=10, n.da=5)
scatter(dapc.3)
# scatterplot shows greater degree of separation with 85 PCs than with 10 PCs

###===================================================================================================================================================###
###  Working Clusters -- Last one, and definitive value of K=6
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# Given the a priori from STRUCTURE (clear assessment of K=6) and the rage of probable values of K within 4-9 from DAPC under different criterion
# We retained K=6 and the best model was n?57
k.native=test1c.k[[383]]
k.native$size
save(k.native,file="Data/DAPC/k.native.rda")
load("Data/DAPC/k.native.rda")

xval<-xvalDapc(mat, grp=k.native$grp, n.pca=seq(5,100,5))

dapc.native<-dapc(native, pop=k.native$grp, n.pca=85,n.da=6)
save(dapc.native,file="Data/DAPC/dapc.native.rda")
load("Data/DAPC/dapc.native.rda")

scatter(dapc.native)

compoplot(dapc.native, posi="top", col=c(rainbow(6,s=1,v=1), rainbow(6,s=0.5,v=0.75)))

scatter(dapc.native, xax=1, yax=2)
scatter(dapc.native, xax=1, yax=3)
scatter(dapc.native, xax=1, yax=4)
# Variance explained
round(dapc.native$eig[1]/sum(dapc.native$eig)*100,digits=1)
round(dapc.native$eig[2]/sum(dapc.native$eig)*100,digits=1)
round(dapc.native$eig[3]/sum(dapc.native$eig)*100,digits=1)
round(dapc.native$eig[4]/sum(dapc.native$eig)*100,digits=1)

###---------------------------------------------------------------------------------------------------------------------------------------------------###
### Study of admixture between inferred clusters
###---------------------------------------------------------------------------------------------------------------------------------------------------###

### Distribution of posterior probabilities across clusters
hist(dapc.native$posterior,breaks=50)
# Posterior probabilities distribution show two extreme values, close to 0 and 1
# Suggesting no admixture between clusters

assignplot(dapc.native, subset=1:300)
# No evidence of admixed individuals with DAPC


###---------------------------------------------------------------------------------------------------------------------------------------------------###
### Identifying Members of Each Cluster
###---------------------------------------------------------------------------------------------------------------------------------------------------###
native.cluster.df=data.frame(Cluster=factor(), Ind.Site=factor(), Ind=factor())

for(g in 1:6){ x=sort(names(k.native$grp)[which(k.native$grp==g)])

for(i in 1:length(x)){ temp.df=data.frame(Cluster=paste("Cluster", g, sep="_"), Ind.Site=strsplit(x[i], "_")[[1]][1], Ind=x[i])

native.cluster.df<-rbind(native.cluster.df, temp.df)

rm(i)
rm(temp.df)}
rm(x) }
rm(g)
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# Create a data frame with all populations allocated to each cluster
unique(native.cluster.df[,1:2])

write.csv(unique(native.cluster.df[,1:2]), "Tables/DAPC/native_cluster.csv", row.names=F,quote=F)





#########################################################################################################################################################
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
###===================================================================================================================================================###
###							Invasive Populations
###===================================================================================================================================================###
###!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!###
#########################################################################################################################################################
# will take multiple approaches:

# 1/ First look at genetic structure in invasive populations, to see if there is genetic divergence
# and/or multiple genetic pools (e.g. suggesting multiple introductions, or stepping-stone dispersal for colonization)
# perform a round of k-means clustering to look at relative divergence since colonization
# divergence from ancestral groups
# divergence w/in the invasive range

# 2/ WORLDWIDE ANALYSIS
# individuals from invasive populations have been assigned to the native clusters
# perhaps will provide some insight into hypotheses to be tested via ABC -> source population(s)

x=c("Aus", "Bel", "Bul1", "Bul2", "Hun", "Ira", "Ita", "Pol", "Spa", "Tur", "UK")
i=1
y=which(Pparva$pop==x[i])
# create an index y of all individuals of Asian populations
for (i in 2:length(x)){ y=c(y, which(Pparva$pop==x[i])) }
rm(i)
rm(x)

invasive<-Pparva
# Retain only individuals from Asian populations
invasive$tab=invasive$tab[y,]
invasive$pop=invasive$pop[y]
invasive$ploidy=invasive$ploidy[y]
rm(y)

save(invasive,file="Data/DAPC/invasive.Rda")
load(file="Data/DAPC/invasive.Rda")


###===================================================================================================================================================###
###  Preliminary Clustering of Data -- k Means
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# run interactively to select the number of clusters (based on minimizing BIC) and to estimate variance explained
# number of retained PCs will be set to maximum allowed for model "stability" (< N/3; 56< 168/3)
# 56 retained PCs appears to exlain 75-80% of cumulative variance
# BIC minimum occurs clearly at K=6, then aplateau at K=7 and increases above K=7

k.invasive<-find.clusters(invasive, max.n=11)

# Compare if clusters corresponded to observed groups (i.e. sampled locations)
table.value(table(pop(invasive), k.invasive$grp), col.lab=paste("inf", 1:11),
            row.lab=paste("ori", 1:11))
# Observed groups globally corresponded to inferred groups, with some outliers

###---------------------------------------------------------------------------------------------------------------------------------------------------###
###  Initial Discriminant Analysis of Principal Components (DAPC)
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# again set to retain maximum numbers of PCs & DA's
# again run interactively to examine scree-plot of F-statistics
# number of discriminant axis is set at (nb of clusters - 1)

dapc.invasive<-dapc(invasive, pop=k.invasive$grp, n.pca=56,n.da=5)

###---------------------------------------------------------------------------------------------------------------------------------------------------###
compoplot(dapc.invasive, posi="top", col=c(rainbow(6,s=1,v=1), rainbow(6,s=0.5,v=0.75)))
# Clear clusters, without shared posterior membership probability between clusters (no 'admixture' evidence within individuals)
dapc.invasive$grp
# no individual show evidence of admixture
# Large pan-european cluster of Aus., Bel., UK, Hun. and Pol.
# Hun. showed signs of membership to 2 clusters (large pan-european and isolated cluster) -> sign of 2 populations? divergence?

# Clusters globally seemed to correspond to existing sampled populations


###===================================================================================================================================================###
###  Alternate Strategy for Clustering of Data:  No. PCs Retained to Achieve Set Cumulative Variance Explained
###---------------------------------------------------------------------------------------------------------------------------------------------------###
#!# am I over-fitting?
# perhaps retain fewer PCs initially
# use cumulative variance as a threshold
###---------------------------------------------------------------------------------------------------------------------------------------------------###
### 75% of variance
###---------------------------------------------------------------------------------------------------------------------------------------------------###
k.invasive<-find.clusters(invasive, max.n=11, pca.select="percVar", perc.pca=75)

# run interactively to select the number of clusters
# BIC minimum occured at k=6
# Same pattern as before, no evidence of overfitting
# plateau with K=7

###---------------------------------------------------------------------------------------------------------------------------------------------------###
### 70% of variance
###---------------------------------------------------------------------------------------------------------------------------------------------------###
k.invasive<-find.clusters(invasive, max.n=11, pca.select="percVar", perc.pca=70)

# run interactively to select the number of clusters
# BIC minimum still occured clearly at k=6
# plateau with K=7

###===================================================================================================================================================###
###  Alternate Strategy for Clustering of Data:  Using a "best fit" Selection Criterion with multiple iterations (in addition to cum. var.)
###---------------------------------------------------------------------------------------------------------------------------------------------------###
test1.k<-list()
for (i in 1:100){ print(paste("Iteration No.", i, "of 100"))
  k.invasive<-find.clusters(invasive, max.n=11, pca.select="percVar", perc.pca=70, choose.n.clust=F, criterion="goodfit", n.iter=1000)
  test1.k[[i]]<-k.invasive
  rm(k.invasive) }
rm(i)

i=1
x<-length(test1.k[[i]]$size)
for(i in 2:100){ x<-c(x, length(test1.k[[i]]$size))}
rm(i)
table(x)
hist(x)
# 100 iterations is not enough
# Mixed scores of 58% for K=4 and 41% for K=5
# BUT "Goodfit" criterion tended to favor a small K, that might not be informative in our case
# Lack of resolution in invasive populations (weak genetic structure in data) or no real clustering?

# Same iterative test but with the criterion "min" of BIC
# As BIC~Nb of clusters decreases up to a plateau and stop decreasing,
# in an exponential negative shape like an island model (see details of ?find.clusters())
test1.k<-list()
for (i in 1:100){ print(paste("Iteration No.", i, "of 100"))
  k.invasive<-find.clusters(invasive, max.n=11, pca.select="percVar", perc.pca=70, choose.n.clust=F, criterion="min", n.iter=1000)
  test1.k[[i]]<-k.invasive
  rm(k.invasive) }
rm(i)

i=1
x<-length(test1.k[[i]]$size)
for(i in 2:100){ x<-c(x, length(test1.k[[i]]$size))}
rm(i)
table(x)
hist(x)
# K=6 for 17 iterations, and K=7 for 83 iterations
# Once again, there is a range between 4 and 7 by considering different criterion

###---------------------------------------------------------------------------------------------------------------------------------------------------###

# probably best to run at least 1,000 total iterations though if I want to be somewhat "Structure-esque"
###---------------------------------------------------------------------------------------------------------------------------------------------------###
for (i in 1:1000){ print(paste("Iteration No.", i, "of 1000"))
  k.invasive=find.clusters(invasive, max.n=19, pca.select="percVar", perc.pca=70, choose.n.clust=F, criterion="goodfit", n.iter=1000, n.start=sample(10:21,1))
  test1.k[[i]]=k.invasive
  rm(k.invasive) }
rm(i)

i=1
x=length(test1.k[[i]]$size)
for(i in 1:1000){ x=c(x, length(test1.k[[i]]$size))}
rm(i)
table(x)
hist(x)
# 4 clusters (n=794/1000) or 5 clusters (n=205/1000) 

for (i in 1:1000){ print(paste("Iteration No.", i, "of 1000"))
  k.invasive=find.clusters(invasive, max.n=19, pca.select="percVar", perc.pca=70, choose.n.clust=F, criterion="min", n.iter=1000, n.start=sample(10:21,1))
  test1.k[[i]]=k.invasive
  rm(k.invasive) }
rm(i)

i=1
x=length(test1.k[[i]]$size)
for(i in 1:1000){ x=c(x, length(test1.k[[i]]$size))}
rm(i)
table(x)
hist(x)
# K=6-7 were favored for minimum criterion (n=119 for K=6 & n=882 for K=7) 
# But might be overfitting

# SO, we tested for over-fitting

###---------------------------------------------------------------------------------------------------------------------------------------------------###
###  Step 2 in Strategy:  Resolving the Over-fitting Problem
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# cross-validation to estimate the optimal number of retained PCs
# given size of dataset, consider 56 as an upper limit
# maximum "allowed" for model stability (< N/3; 56< 168/3)
mat=tab(invasive, NA.method="mean")

# It seemed that at least 20 PCA was sufficient to discriminate amongst clusters
xval=xvalDapc(mat, grp=invasive@pop, n.pca=seq(5,55,5))

# 20 seemed enough to discriminate amongst k=7 clusters
xval=xvalDapc(mat, grp=test1.k[[10]]$grp, n.pca=seq(5,55,5))

#!# NB
# NONE of the plots displayed a nice peaked pattern
###---------------------------------------------------------------------------------------------------------------------------------------------------###
# a-score optimization for K=7 -> 'min' criterion that increased the number of clusters compared to 'goodfit' criterion
# Linear decline of a-score from 0 to 55
# suggested far fewer (<10) PCs to avoid overfitting
dapc.1=dapc(invasive, pop=test1.k[[10]]$grp, pca.select="percVar", perc.pca=70, n.da=6)
scatter(dapc.1)
optim.a.score(dapc.1, n.pca=seq(5,55,5))


dapc.2=dapc(invasive, pop=test1.k[[10]]$grp, pca.select="percVar", n.pca=56, n.da=6)
scatter(dapc.2)
dapc.3=dapc(invasive, pop=test1.k[[10]]$grp, pca.select="percVar", n.pca=10, n.da=6)
scatter(dapc.3)
# scatterplot showed greater degree of separation with 55 PCs than with 10 PCs -> we should decide to keep the max. of PCs (i.e. 56)


###===================================================================================================================================================###
###  Working Clusters
###---------------------------------------------------------------------------------------------------------------------------------------------------###
for (k in 4:5) {
  k.invasive<-find.clusters(invasive, max.n=11,n.pca=56,n.clust=k)
  dapc.invasive<-dapc(invasive, pop=k.invasive$grp, n.pca=55,n.da=k-1)
  scatter(dapc.invasive, xax=1, yax=2)
  scatter(dapc.invasive, xax=1, yax=3)
  # Variance explained
  round(dapc.invasive$eig[1]/sum(dapc.invasive$eig)*100,digits=1)
  round(dapc.invasive$eig[2]/sum(dapc.invasive$eig)*100,digits=1)
  round(dapc.invasive$eig[3]/sum(dapc.invasive$eig)*100,digits=1)
  invasive.cluster.df<-data.frame(Cluster=factor(), Ind.Site=factor(), Ind=factor())
  for(g in 1:k){ x<-sort(names(k.invasive$grp)[which(k.invasive$grp==g)])
  for(i in 1:length(x)){ temp.df<-data.frame(Cluster=paste("Cluster", g, sep="_"), Ind.Site=strsplit(x[i], "_")[[1]][1], Ind=x[i])
  invasive.cluster.df<-rbind(invasive.cluster.df, temp.df)
  rm(i)
  rm(temp.df)}
  rm(x) }
  rm(g)
  unique(invasive.cluster.df[,1:2])
  compoplot(dapc.invasive, posi="top")
}

# K=4 was a parcimonious number of K and the lowest BIC value (goodfit criterion)
# BIC for K=5 was close to the BIC value of K=4
k.invasive<-find.clusters(invasive, max.n=11,n.pca=56,n.clust=4)
k.invasive$size
save(k.invasive,file="Data/DAPC/k.invasive.rda")
load(file="Data/DAPC/k.invasive.rda")

xval<-xvalDapc(mat, grp=k.invasive$grp, n.pca=seq(5,55,5))

dapc.invasive<-dapc(invasive, pop=k.invasive$grp, n.pca=55,n.da=3)
save(dapc.invasive,file="Data/DAPC/dapc.invasive.rda")
load("Data/DAPC/dapc.invasive.rda")

scatter(dapc.invasive, xax=1, yax=2)
scatter(dapc.invasive, xax=1, yax=3)
# Variance explained
round(dapc.invasive$eig[1]/sum(dapc.invasive$eig)*100,digits=1)
round(dapc.invasive$eig[2]/sum(dapc.invasive$eig)*100,digits=1)
round(dapc.invasive$eig[3]/sum(dapc.invasive$eig)*100,digits=1)

compoplot(dapc.invasive, posi="top", col=c(rainbow(6,s=1,v=1), rainbow(6,s=0.5,v=0.75)))


###---------------------------------------------------------------------------------------------------------------------------------------------------###
### Identifying Members of Each Cluster
###---------------------------------------------------------------------------------------------------------------------------------------------------###
invasive.cluster.df<-data.frame(Cluster=factor(), Ind.Site=factor(), Ind=factor())

for(g in 1:4){ x<-sort(names(k.invasive$grp)[which(k.invasive$grp==g)])

for(i in 1:length(x)){ temp.df<-data.frame(Cluster=paste("Cluster", g, sep="_"), Ind.Site=strsplit(x[i], "_")[[1]][1], Ind=x[i])

invasive.cluster.df<-rbind(invasive.cluster.df, temp.df)

rm(i)
rm(temp.df)}
rm(x) }
rm(g)


###---------------------------------------------------------------------------------------------------------------------------------------------------###
# Create a data frame with all populations allocated to each cluster
unique(invasive.cluster.df[,1:2])

write.csv(unique(invasive.cluster.df[,1:2]), "Tables/DAPC/invasive_cluster.csv", row.names=F,quote=F)

###---------------------------------------------------------------------------------------------------------------------------------------------------###
# To look for a given cluster

subset((unique(invasive.cluster.df[,1:2])), Cluster=="Cluster_3")
barplot(table(subset(invasive.cluster.df, Cluster=="Cluster_3")$Ind.Site), las=1)

barplot(sort(table(subset(invasive.cluster.df, Cluster=="Cluster_3")$Ind.Site)), las=1)

###---------------------------------------------------------------------------------------------------------------------------------------------------###
# Ou bien pour identifier chaque cluster ou se trouvent des echantillons d'une site en particulier
subset((unique(invasive.cluster.df[,1:2])), Ind.Site=="Aus")

x=c("Aus", "Bel", "Bul1", "Bul2", "Hun", "Ira", "Ita", "Pol", "Spa", "Tur", "UK")
for (i in 1:length(x)){ print(i)
  print(subset((unique(invasive.cluster.df[,1:2])), Ind.Site==x[i])) }
rm(i)
rm(x)



###===================================================================================================================================================###
###===================================================================================================================================================###
###
###         WORLDWIDE ANALYSIS
###
###  2/     Assignments of invasive populations to native Sampling Site
###===================================================================================================================================================###
###===================================================================================================================================================###

###===================================================================================================================================================###
###  Assignments by Sampling Site
###---------------------------------------------------------------------------------------------------------------------------------------------------###

###---------------------------------------------------------------------------------------------------------------------------------------------------###
# Compute assignment probability to native clusters for any invasive population in the vector 'invasive.names'
for (p in invasive.names) {
  cat("Population assignment for", p,"\n")
  temp.genind=Pparva[pop=p]
  
  pred=predict.dapc(dapc.native, newdata=temp.genind)
  table(tail(pred$assign, length(which(Pparva$pop==p))))
  barplot(t((tail(pred$posterior, length(which(Pparva$pop==p))))), legend=T, args.legend=list(x="right", horiz=F),
          col=c(rainbow(6,s=1,v=1), rainbow(3,s=0.5,v=0.75)),main=paste("Assignment of",p,"individuals to native clusters",sep=" "))
  rm(p)
  rm(temp.genind)
}

###---------------------------------------------------------------------------------------------------------------------------------------------------###
### All in a Single Plot
###---------------------------------------------------------------------------------------------------------------------------------------------------###
png("Figures/Pop.assignment.InvasiveToNative.png",width=1600,height = 1050)
par(mfrow=c(4,3))
par(oma=c(1,1,2,1))
par(mar=c(2,3,2,0))

for (p in invasive.names) {
  cat("Population assignment for", p,"\n")
  temp.genind=Pparva[pop=p]
  
  pred=predict.dapc(dapc.native, newdata=temp.genind)
  table(tail(pred$assign, length(which(Pparva$pop==p))))
  barplot(t((tail(pred.Aus$posterior, length(which(Pparva$pop==p))))), legend=F,
          col=c(rainbow(6,s=1,v=1), rainbow(3,s=0.5,v=0.75)),main=paste("Individuals of",p,sep=" "))
  rm(p)
  rm(temp.genind)
}
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("left", legend=(1:7), fill=c(rainbow(6,s=1,v=1), rainbow(3,s=0.5,v=0.75)),cex=2)
mtext("Posterior Probability of Assignment", side=3, outer=T, line=0, font=2)

par(mfrow=c(1,1))
par(oma=c(1,1,1,1))
par(mar=c(1,1,1,1))
dev.off()



#======================================
# THE END
#======================================