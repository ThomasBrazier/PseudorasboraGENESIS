############################################################################
#   INRA - ANR project Pseudorasbora
#
#       FIGURES
#
############################################################################
# by Thomas Brazier
# brazier.thomas@gmail.com
# MSc internship

# Supervisor: Scott McCairns
#             INRA UMR ESE


# DESCRIPTION
# This script is the main workflow for producing figures that will serve in the report (beautiful plots)

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
# SAMPLED SITES
#----------------------------------------------------------

#----------------------------------
# WORLDWIDE MAP OF SAMPLED SITES FOR M&M (ALL IN ONE MAP)
# Download database of river network at scale 50
# rivers50=ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')
# save(rivers50,file="Data/rivers50.Rda")
load("Data/rivers50.Rda")

png("Figures/Worldwide M&M map.png",width=3000,height=1000)
map("world", xlim=c(-15,150), ylim=c(20,59), col="gray90",  fill=TRUE) # Plot the map of European area
# Nomenclature
north.arrow(149,56,2, lab="", lab.pos="above")
scalebar(c(125,20),2000, bg=NA, border=NA, division.cex=3)
points(invasive.coords$X, invasive.coords$Y, pch = 16, cex = 4,col="grey49") # Add points to the map
# Labels of sampled sites
# Jitter the 4rth label
text(invasive.coords$X, y = invasive.coords$Y, invasive.coords$Pop,
     pos = c(rep(4,3), 1, rep(4,7)), offset = 1, cex = 3.6, font = 2)
points(native.coords$X, native.coords$Y, pch = 16, cex = 4,col="black") # Add points to the map
points(native.coords$X[c(11,13,18)], native.coords$Y[c(11,13,18)], pch = 16, cex = 4,col="grey49") #\ Add points to the map
# Labels of sampled sites
# Jitter the 13rd label
text(native.coords$X, y = native.coords$Y, native.coords$Pop,
     pos = c(rep(4,12), 1, rep(4,5)), offset = 1, cex = 3.6, font = 2)
dev.off()


# For maps, Bul1 and Bul2 pie charts are displayed on the same point
# S19 and S20 too
# Jitter these populations coordinates
native.coords[native.coords$Pop=="S19",3]=native.coords[native.coords$Pop=="S19",3]-2
native.coords[native.coords$Pop=="S20",3]=native.coords[native.coords$Pop=="S19",3]+2

invasive.coords[invasive.coords$Pop=="Bul1",3]=invasive.coords[invasive.coords$Pop=="Bul1",3]-1
invasive.coords[invasive.coords$Pop=="Bul2",3]=invasive.coords[invasive.coords$Pop=="Bul2",3]+1




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
# DAPC results & figures (explanations in DAPC.R)
#----------------------------------------------------------

# Number of K and DAPC have been computed in 'DAPC.R'. Results were loaded in .rda objects for figures purposes


#----------------------------------------------------------
# BEAUTIFUL PLOTS
#----------------------------------------------------------

# Load the K value selected in DAPC.R
load("Data/DAPC/native.rda") # Data set genind for DAPC
load("Data/DAPC/k.native.rda") # Clusters
load("Data/DAPC/dapc.native.rda") # Results of the dapc for 6 clusters in native populations
# Set a vector of colors for clusters
# cols.native=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
#               "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9")

cols.native=c("#FF7F00",  "#377EB8", "#4DAF4A", "#FFFF33", "#984EA3")
# Load the K value selected in DAPC.R
# load("Data/DAPC/invasive.rda") # Data set genind for DAPC
# load("Data/DAPC/k.invasive.rda") # Clusters
# load("Data/DAPC/dapc.invasive.rda") # Results of the dapc for 6 clusters in native popualtions
# # Set a vector of colors for the number of clusters
# # cols.invasive=brewer.pal(n = length(levels(dapc.invasive$grp)), name = "Set2")
# cols.invasive=c("#A65628", "#F781BF",
#               "#8DD3C7", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9")
# Custom K value
# K = 6
# k.native<-find.clusters(native, max.n=11,n.pca=56,n.clust=K)
# dapc.native<-dapc(native, pop=k.native$grp, n.pca=56,n.da=(K-1))

#------------
# Fig. Scatter plot of 2 first axis of DAPC with clusters
# With barplot of clusters
png("Figures/DAPC/Final/DAPC_Native_Scatter.png",width=600,height=600)
scatter.dapc(dapc.native,xax=1,yax=2,addaxes=TRUE,xlab="Axe1",ylab="Axe 2",grp=dapc.native$grp,col=cols.native, cex=2, bg="white",cstar=1,cellipse=0.95,
             csub =2, clabel = 2, posi.da="topright",
             sub=paste("PC1 (",round(dapc.native$eig[1]/sum(dapc.native$eig)*100,digits=1),"%) vs PC2 (",round(dapc.native$eig[2]/sum(dapc.native$eig)*100,digits=1),"%)",sep=""))
dev.off()

# Custom barplot (Posterior membership probability as a function of individuals)
# Format a data frame with 4 columns, and as many rows per individual as there is clusters
dapc.results=as.data.frame(dapc.native$posterior)
dapc.results$pop=pop(native)
dapc.results$indNames=rownames(dapc.results)
dapc.results=melt(dapc.results)
colnames(dapc.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")

######### GGPLOT
# Cluster populations together in reference demes before plotting...
dapc.results$Original_Population=factor(dapc.results$Original_Population,
                                        levels = c("S9","S18","S19","S20","Tib","S10","S11",
                                                   "S4","S6",
                                                   "S3","S1","S2","S16","S13","S14",
                                                   "S15","S17", "Jap"))# Plotting the barplot
dapc_barplot=ggplot(data=dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.native) +
  facet_grid(~Original_Population, scales = "free",space="free_x") +
  labs(x="Sampled individual", y="Posterior membership probabilities\n", fill="Assigned\npopulation") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=22, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=22,hjust = 0.5),
        axis.title.x = element_text(color="black", size=22),
        axis.title.y = element_text(color="black", size=22),
        axis.text=element_text(size=22, colour="black"),
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=18, colour="black"),
        # legend.key = element_rect(fill = "white", size = 1),
        # legend.text=element_text(size=22),
        # legend.title=element_text(size=22),
        legend.position = "none")
dapc_barplot

ggsave(paste(figuresdir,"/DAPC/Final/Native_Posterior_membership_probability_barplot_K5.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)

scat.native = as.ggplot(expression(scatter.dapc(dapc.native,xax=1,yax=2,addaxes=TRUE,xlab="Axe1",ylab="Axe 2",grp=dapc.native$grp,col=cols.native, cex=2, bg="white",cstar=1,cellipse=0.95,
                                                csub =2, clabel = 2, posi.da="topright",
                                                sub=paste("PC1 (",round(dapc.native$eig[1]/sum(dapc.native$eig)*100,digits=1),"%) vs PC2 (",round(dapc.native$eig[2]/sum(dapc.native$eig)*100,digits=1),"%)",sep=""))))

ggarrange(scat.native, dapc_barplot, ncol = 2, nrow = 1, widths = c(1,2), labels = "auto",
          font.label = list(size = 30))
ggsave(paste(figuresdir,"/DAPC/Final/Native_Scatter&Barplot.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=20)




#----------------------------------------------------------
# NATIVE AREA
#----------------------------------------------------------

# Load the K value selected in DAPC.R
load("Data/DAPC/native.rda") # Data set genind for DAPC
load("Data/DAPC/k.native.rda") # Clusters
load("Data/DAPC/dapc.native.rda") # Results of the dapc for 6 clusters in native populations

# or run a custom value of K
K = 5
k.native = find.clusters(native, max.n=19,n.pca=100,n.clust=K)
dapc = dapc(native, pop=k.native$grp, n.pca=40,n.da=(K-1))


dapc.native
summary(dapc.native)

# Contribution of alleles to the variance in the data set
loadingplot(dapc.native$var.contr, axis = 2, thres = 0.07, lab.jitter = 1)
# Pattern of contribution seemed neutral (random pattern) with multiple contributors to the variance
sum(dapc.native$var.contr>0.005)
# 236 alleles contributed to more than 0.005 of the explained variance
1/ncol(native$tab)
# null hypothesis of equal contribution was 1/nb of alleles = 0.00024
sum(dapc.native$var.contr>0.00024)
sum(dapc.native$var.contr>0.00024)/ncol(native$tab)
# 2456 alleles/4216 (58%) contributed more than expected by chance

####### Figures #######

### BATCH PLOT FOR K=2-11
# Download database of river network at scale 50
# rivers50=ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')
# save(rivers50,file="Data/rivers50.Rda")
load("Data/rivers50.Rda")

for (p in 2:11) { # p is the value of K to plot, best K canidates
  # Run DAPC for each K
  k=find.clusters(native, max.n=19,n.pca=100,n.clust=p)
  dapc=dapc(native, pop=k$grp, n.pca=85,n.da=(p-1))
  
  # Set a vector of colors for clusters
  cols.native=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
                "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9")
  
  #------------
  # Fig 1. Bar plot of admixture proportions per individual, sorted by populations
  # Custom barplot (Posterior membership probability as a function of individuals)
  # Format a data frame with 4 columns, and as many rows per individual as there is clusters
  dapc.results=as.data.frame(dapc$posterior)
  dapc.results$pop=pop(native)
  dapc.results$indNames=rownames(dapc.results)
  dapc.results=melt(dapc.results)
  colnames(dapc.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
  
  ######### GGPLOT
  # Cluster populations together in reference demes before plotting...
  dapc.results$Original_Population=factor(dapc.results$Original_Population,
                                          levels = c("S9","S18","S19","S20","Tib","S10","S11",
                                                     "S4","S6",
                                                     "S3","S1","S2","S16","S13","S14",
                                                     "S15","S17", "Jap"))
  # Plotting the barplot
  dapc_barplot=ggplot(data=dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
    geom_bar(stat='identity',width=1) +
    scale_fill_manual(values = cols.native) +
    facet_grid(~Original_Population, scales = "free") +
    labs(x="Sampled individual", y="Posterior\nmembership\nprobabilities\n", fill="Assigned\npopulation") +
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
          axis.text.x=element_blank(), # No samples names
          axis.ticks.x=element_blank(), # No x axis
          strip.text=element_text(size=28, colour="black"),
          legend.key = element_rect(fill = "white", size = 1),
          legend.text=element_text(size=28),
          legend.title=element_text(size=28))
  dapc_barplot
  
  ggsave(paste(figuresdir,"/DAPC/Batch2/Native_Posterior_membership_probability_barplot_K",p,".png",sep=""),
         device="png",dpi=320,units="cm",width=60,height=9)
  # ggsave(paste(figuresdir,"/DAPC/Batch/Native_Posterior_membership_probability_barplot_K",p,".jpg",sep=""),
  #        device="jpg",dpi=150,units="cm",width=60,height=9)
  
  #------------
  # Fig 2. Scatter plot of 2 first axis of DAPC with clusters
  png(paste("Figures/DAPC/Batch2/DAPC_scatter_native_pops_K",p,".png",sep=""),width=800,height=500)
  scatter.dapc(dapc,xax=1,yax=2,addaxes=TRUE,xlab="Axe1",ylab="Axe 2",grp=dapc$grp,col=cols.native, cex=2, bg="white",cstar=1,cellipse=0.95,
               sub=paste("PC1 (",round(dapc$eig[1]/sum(dapc$eig)*100,digits=1),"%) vs PC2 (",round(dapc$eig[2]/sum(dapc$eig)*100,digits=1),"%)",sep=""))
  dev.off()
  
  #------------
  # Fig 3. Map of the native area with pies chart of cluster membership on each sampling location
  # Prepare the data set for mapping
  Native.pies=cbind(as.character(pop(native)),as.data.frame(dapc$posterior))
  names(Native.pies)=c("Population", paste("Cluster ",c(1:p)))
  Native.pie_means=sapply(split(Native.pies[2:(1+p)],Native.pies$Population), colMeans) ## This just makes population cluster averages for each cluster
  # Mapping
  png(paste("Figures/DAPC/Batch2/Native_pies_map_K",p,".png",xsep=""),width=2000,height=1100)
  map("worldHires", xlim=c(92, 160), ylim=c(20,48), col="gray90",  fill=TRUE) # Plot the map of Asian area
  sp::plot(rivers50,lwd=2, col = 'blue',add=TRUE)
  # Nomenclature
  north.arrow(155,45,1, lab="N", lab.pos="above")
  scalebar(c(150,20),1000, bg=NA, border=NA, division.cex=2)
  text(x=105, y=34,"China",cex=3,font=2)
  text(x=139, y=39,"Japan",cex=3,font=2)
  points(native.coords$X, native.coords$Y, pch = 16, cex = 0.7,col="red") #\ Add points to the map
  ## Pies ##
  # Draw a pie with proportions of membership for each population
  for (i in 1:ncol(Native.pie_means)) {
    add.pie(Native.pie_means[,i],x=native.coords$X[i], y=native.coords$Y[i],labels=native.coords[i,1],label.dist=1.5,
            radius=1,edges=200,clockwise=T,	col	=	cols.native,cex=2)
    # Add labels to pies
  }
  dev.off()
}

################# INDIVIDUAL PLOTS
# Set a vector of colors for clusters
cols.native=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
              "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9")
# Fig 1. Bar plot of admixture proportions per individual, sorted by populations
png("Figures/DAPC/DAPC_barplot_native_pops.png",width=800,height=400)
compoplot(dapc.native,col=cols.native)
dev.off()

# Custom barplot (Posterior membership probability as a function of individuals)
# Format a data frame with 4 columns, and as many rows per individual as there is clusters
dapc.results=as.data.frame(dapc.native$posterior)
dapc.results$pop=pop(native)
dapc.results$indNames=rownames(dapc.results)
dapc.results=melt(dapc.results)

colnames(dapc.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")

######### GGPLOT
# Cluster populations together in reference demes before plotting...
dapc.results$Original_Population=factor(dapc.results$Original_Population,
                                        levels = c("S9","S18","S19","S20","Tib","S10","S11",
                                                   "S4","S6",
                                                   "S3","S1","S2","S16","S13","S14",
                                                   "S15","S17", "Jap"))
# Plotting the barplot
dapc_barplot=ggplot(data=dapc.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.native) +
  facet_grid(~Original_Population, scales = "free",space="free_x") +
  labs(x="Sampled individual", y="Posterior\nmembership\nprobabilities\n", fill="Assigned\npopulation") +
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
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=20, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
dapc_barplot

ggsave(paste(figuresdir,"/DAPC/Native_Posterior_membership_probability_barplot.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)

# Fig 2. Scatter plot of 2 first axis of DAPC with clusters
png("Figures/DAPC/DAPC_scatter_native_pops.png",width=800,height=500)
scatter.dapc(dapc.native,xax=1,yax=2,addaxes=TRUE,xlab="Axe1",ylab="Axe 2",grp=dapc.native$grp,col=cols.native, cex=2, bg="white",cstar=1,cellipse=0.95,
             sub=paste("PC1 (",round(dapc.native$eig[1]/sum(dapc.native$eig)*100,digits=1),"%) vs PC2 (",round(dapc.native$eig[2]/sum(dapc.native$eig)*100,digits=1),"%)",sep=""))
dev.off()
# Look at other axes of DAPC vs PC1
png("Figures/DAPC/DAPC_scatter_native_pops_all_axes.png",width=1600,height=1500)
par(mfrow=c(2,2))
scatter.dapc(dapc.native,xax=1,yax=2,grp=dapc.native$grp,col=cols.native, cex=2,csub=4, clab=2, bg="white",cstar=1,cellipse=0.95,
             sub=paste("PC1 (",round(dapc.native$eig[1]/sum(dapc.native$eig)*100,digits=1),"%) vs PC2 (",round(dapc.native$eig[2]/sum(dapc.native$eig)*100,digits=1),"%)",sep=""))
scatter.dapc(dapc.native,xax=1,yax=3,grp=dapc.native$grp,col=cols.native, cex=2,csub=4, clab=2, bg="white",cstar=1,cellipse=0.95,
             sub=paste("PC1 (",round(dapc.native$eig[1]/sum(dapc.native$eig)*100,digits=1),"%) vs PC3 (",round(dapc.native$eig[3]/sum(dapc.native$eig)*100,digits=1),"%)",sep=""))
scatter.dapc(dapc.native,xax=1,yax=4,grp=dapc.native$grp,col=cols.native, cex=2,csub=4, clab=2, bg="white",cstar=1,cellipse=0.95,
             sub=paste("PC1 (",round(dapc.native$eig[1]/sum(dapc.native$eig)*100,digits=1),"%) vs PC4 (",round(dapc.native$eig[4]/sum(dapc.native$eig)*100,digits=1),"%)",sep=""))
scatter.dapc(dapc.native,xax=1,yax=5,grp=dapc.native$grp,col=cols.native, cex=2,csub=4, clab=2, bg="white",cstar=1,cellipse=0.95,
             sub=paste("PC1 (",round(dapc.native$eig[1]/sum(dapc.native$eig)*100,digits=1),"%) vs PC5 (",round(dapc.native$eig[5]/sum(dapc.native$eig)*100,digits=1),"%)",sep=""))
par(mfrow=c(1,1))
dev.off()

# Fig 3. Map of the native area with pies chart of cluster membership on each sampling location

# Download database of river network at scale 50
# rivers50=ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')

# Prepare the data set for mapping
Native.pies=cbind( as.character(pop(native)),as.data.frame(dapc.native$posterior))
names(Native.pies)=c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")

Native.pie_means=sapply(split(Native.pies[2:7],Native.pies$Population), colMeans) ## This just makes population cluster averages for each cluster

# Mapping
png("Figures/DAPC/Native_pies_map.png",width=2000,height=1100)
map("worldHires", xlim=c(92, 160), ylim=c(20,48), col="gray90",  fill=TRUE) # Plot the map of Asian area
sp::plot(rivers50,lwd=2, col = 'blue',add=TRUE)
# Nomenclature
north.arrow(155,45,1, lab="N", lab.pos="above")
scalebar(c(150,20),1000, bg=NA, border=NA, division.cex=2)

text(x=105, y=34,"China",cex=3,font=2)
text(x=139, y=39,"Japan",cex=3,font=2)
points(native.coords$X, native.coords$Y, pch = 16, cex = 0.7,col="red") #\ Add points to the map
# map.scale(x=140,y=25,relwidth=0.1,cex=2)

## Pies ## 
# Draw a pie with proportions of membership for each population
for (i in 1:ncol(Native.pie_means)) {
  add.pie(Native.pie_means[,i],x=native.coords$X[i], y=native.coords$Y[i],labels=native.coords[i,1],label.dist=1.5,
          radius=1,edges=200,clockwise=T,	col	=	cols.native,cex=2,font=2)
  # Add labels to pies
}
dev.off()



#----------------------------------------------------------
# INVASIVE AREA
#----------------------------------------------------------

#### Plots of invasive area independent clustering (> Supplementary material)

# Load the K value selected in DAPC.R
load("Data/DAPC/invasive.rda") # Data set genind for DAPC
load("Data/DAPC/k.invasive.rda") # Clusters
load("Data/DAPC/dapc.invasive.rda") # Results of the dapc for 6 clusters in native popualtions

# or run a custom value of K
k.invasive<-find.clusters(invasive, max.n=11,n.pca=56,n.clust=6)
dapc.invasive<-dapc(invasive, pop=k.invasive$grp, n.pca=56,n.da=5)


dapc.invasive
summary(dapc.invasive)

# Contribution of alleles to the variance in the data set
loadingplot(dapc.invasive$var.contr, axis = 2, thres = 0.07, lab.jitter = 1)
# Pattern of contribution seemed neutral (random pattern) with multiple contributors to the variance
sum(dapc.invasive$var.contr>0.005)
# 144 alleles contributed to more than 0.005 of the explained variance
1/ncol(invasive$tab)
# null hypothesis of equal contribution was 1/nb of alleles = 0.00024
sum(dapc.invasive$var.contr>0.00024)
sum(dapc.invasive$var.contr>0.00024)/ncol(invasive$tab)
# 1434 alleles/4216 (34%) contributed more than expected by chance

####### Figures #######
# Set a vector of colors for the number of clusters
cols.invasive=brewer.pal(n = length(levels(dapc.invasive$grp)), name = "Set2")

# Fig 1. Bar plot of admixture proportions per individual, sorted by populations
png("Figures/DAPC/DAPC_barplot_invasive_pops.png",width=800,height=400)
compoplot(dapc.invasive,col=cols.invasive)
dev.off()


# Fig 2. Scatter plot of 2 first axis of DAPC with clusters
png("Figures/DAPC/DAPC_scatter_invasive_pops.png",width=800,height=500)
scatter.dapc(dapc.invasive,xax=1,yax=2,grp=dapc.invasive$grp,col=cols.invasive, cex=2, bg="white",cstar=1,cellipse=0.95,
             sub=paste("PC1 (",round(dapc.invasive$eig[1]/sum(dapc.invasive$eig)*100,digits=1),"%) vs PC2 (",round(dapc.invasive$eig[2]/sum(dapc.invasive$eig)*100,digits=1),"%)",sep=""))
dev.off()
# Look at other axes of DAPC vs PC1
png("Figures/DAPC/DAPC_scatter_invasive_pops_all_axes.png",width=800,height=1000)
par(mfrow=c(2,1))
scatter.dapc(dapc.invasive,xax=1,yax=2,grp=dapc.invasive$grp,col=cols.invasive, cex=2,csub=2, clab=2, bg="white",cstar=1,cellipse=0.95,
             sub=paste("PC1 (",round(dapc.invasive$eig[1]/sum(dapc.invasive$eig)*100,digits=1),"%) vs PC2 (",round(dapc.invasive$eig[2]/sum(dapc.invasive$eig)*100,digits=1),"%)",sep=""))
scatter.dapc(dapc.invasive,xax=1,yax=3,grp=dapc.invasive$grp,col=cols.invasive, cex=2,csub=2, clab=2, bg="white",cstar=1,cellipse=0.95,
             sub=paste("PC1 (",round(dapc.invasive$eig[1]/sum(dapc.invasive$eig)*100,digits=1),"%) vs PC3 (",round(dapc.invasive$eig[3]/sum(dapc.invasive$eig)*100,digits=1),"%)",sep=""))
par(mfrow=c(1,1))
dev.off()

# Fig 3. Map of the native area with pies chart of cluster membership on each sampling location

# Prepare the data set for mapping
Invasive.pies=cbind( as.character(pop(invasive)),as.data.frame(dapc.invasive$posterior))
names(Invasive.pies)=c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4")

Invasive.pie_means=sapply(split(Invasive.pies[2:5],Invasive.pies$Population), colMeans) ## This just makes population cluster averages for each cluster

# Mapping
png("Figures/DAPC/Invasive_pies_map.png",width=2000,height=1200)
map("worldHires", xlim=c(-15,60), ylim=c(30,59), col="gray90",  fill=TRUE) # Plot the map of European area
map(database="worldHires", regions="Caspian Sea", col="white", fill=TRUE, add=TRUE) # Correct for unfilled Caspian sea
map(database="worldHires", regions="Aral Sea", col="white", fill=TRUE, add=TRUE)# Correct for unfilled Aral sea
# map('rivers',add=TRUE,xlim=c(-15,30), ylim=c(30,59),col="blue",lwd=2,resolution=2) # European rivers
# map('rivers',add=TRUE,xlim=c(35,50), ylim=c(30,40),col="blue",lwd=2,resolution=2) # Turkish + Iranian rivers
sp::plot(rivers50,xlim=c(-15,30), ylim=c(30,59),lwd=2, col = 'blue',add=TRUE)
# Nomenclature
north.arrow(56,56,1, lab="N", lab.pos="above")
scalebar(c(49,31),1000, bg=NA, border=NA, division.cex=2)


text(x=105, y=34,"",cex=3,font=2)
text(x=139, y=39,"",cex=3,font=2)
points(invasive.coords$X, invasive.coords$Y, pch = 16, cex = 0.7,col="red") #\ Add points to the map
# map.scale(x=50,y=55,relwidth=0.1,cex=2)

## Pies ## 
# Draw a pie with proportions of membership for each population
for (i in 1:ncol(Invasive.pie_means)) {
  add.pie(Invasive.pie_means[,i],x=invasive.coords$X[i], y=invasive.coords$Y[i],labels=invasive.coords[i,1],label.dist=1.5,
          radius=1,edges=200,clockwise=T,	col	=	cols.invasive,cex=2,font=2)
  # Add labels to pies
}
dev.off()




#----------------------------------------------------------
# WORLDWIDE CLUSTERING
#----------------------------------------------------------

#### Plot of invasive area populations assigned to native area clusters (> publications - Results)
###---------------------------------------------------------------------------------------------------------------------------------------------------###
### All Posterior Probabilities of Assignment in a Single Plot
###---------------------------------------------------------------------------------------------------------------------------------------------------###
png("Figures/DAPC/Pop_assignment_InvasiveToNative.png",width=1600,height = 1050)
par(mfrow=c(4,3))
par(oma=c(2,2,3,2))
par(mar=c(3,3,3,3))

dapc.posterior.membership=data.frame()

for (p in invasive.names) {
  cat("Population assignment for", p,"\n")
  temp.genind=Pparva[pop=p]
  
  pred=predict.dapc(dapc.native, newdata=temp.genind)
  table(tail(pred$assign, length(which(Pparva$pop==p))))
  barplot(t((tail(pred$posterior, length(which(Pparva$pop==p))))), legend=F,
          col=cols.native,main=paste("Individuals of",p,sep=" "),cex.main=3)
  
  dapc.posterior.membership=rbind(dapc.posterior.membership,cbind(rep(as.character(p), nrow(pred$posterior)),pred$posterior))
  
  rm(p)
  rm(temp.genind)
}
plot(NULL ,xaxt='n',yaxt='n',bty='n',ylab='',xlab='', xlim=0:1, ylim=0:1)
legend("left", legend=(1:7), fill=cols.native,cex=2)
mtext("Posterior Probability of Assignment", side=3, outer=T, line=0, font=2,cex=3)

par(mfrow=c(1,1))
par(oma=c(1,1,1,1))
par(mar=c(1,1,1,1))
dev.off()


# Fig 3. Map of the native area with pies chart of cluster membership on each sampling location

# Prepare the data set for mapping
names(dapc.posterior.membership)=c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")
for (i in 2:7) {
  dapc.posterior.membership[,i]=as.numeric(as.character(dapc.posterior.membership[,i]))
}

dapc.posterior.membership_means=sapply(split(dapc.posterior.membership[2:7],dapc.posterior.membership$Population), colMeans) ## This just makes population cluster averages for each cluster

# Mapping
png("Figures/DAPC/Invasive_pies_posteriorMembership_map.png",width=2000,height=1200)
map("worldHires", xlim=c(-15,60), ylim=c(30,59), col="gray90",  fill=TRUE) # Plot the map of European area
map(database="worldHires", regions="Caspian Sea", col="white", fill=TRUE, add=TRUE) # Correct for unfilled Caspian sea
map(database="worldHires", regions="Aral Sea", col="white", fill=TRUE, add=TRUE)# Correct for unfilled Aral sea
# map('rivers',add=TRUE,xlim=c(-15,30), ylim=c(30,59),col="blue",lwd=2,resolution=2) # European rivers
# map('rivers',add=TRUE,xlim=c(35,50), ylim=c(30,40),col="blue",lwd=2,resolution=2) # Turkish + Iranian rivers
sp::plot(rivers50,xlim=c(-15,30), ylim=c(30,59),lwd=2, col = 'blue',add=TRUE)
# Nomenclature
north.arrow(56,56,1, lab="N", lab.pos="above")
scalebar(c(49,31),1000, bg=NA, border=NA, division.cex=2)


text(x=105, y=34,"",cex=3,font=2)
text(x=139, y=39,"",cex=3,font=2)
points(invasive.coords$X, invasive.coords$Y, pch = 16, cex = 0.7,col="red") #\ Add points to the map
# map.scale(x=50,y=55,relwidth=0.1,cex=2)

## Pies ## 
# Draw a pie with proportions of membership for each population
for (i in 1:ncol(dapc.posterior.membership_means)) {
  add.pie(dapc.posterior.membership_means[,i],x=invasive.coords$X[i], y=invasive.coords$Y[i],labels=invasive.coords[i,1],label.dist=1.5,
          radius=1,edges=200,clockwise=T,	col	=	cols.native,cex=2,font=2)
  # Add labels to pies
}
dev.off()




#==========================================================
# Part 3. STRUCTURE analyses
#==========================================================

#-----------------------------------
# NATIVE POPULATIONS


#########################
# All populations at once
# 1/ Graphic evaluation with the four plots of Evanno's method (Evanno et al. 2005)
# Most important plot is DeltaK = mean(|L''(K)|)/sd(L(K)) of L(K) as a function of K 

### Choose the most probable K value (number of cluster inferred), K value the most appropriate for the data
# Load results from Harvester: 'evanno.txt'
# evanno=read.table(paste(STRdir,"/2a_Native/Harvester/evanno.txt",sep=""),header=F,sep="\t")
# evanno=read.table(paste(STRdir,"/2b_Native/Harvester/evanno.txt",sep=""),header=F,sep="\t")
evanno=read.table(paste(STRdir,"/2b_Native_OC/Harvester/evanno.txt",sep=""),header=F,sep="\t")

# With K's assessment on convergent replciates only, Delta K suggested K=17 for 18 populations --> may be overclustering
# Look for more reasonable values under K=17
evanno=evanno[1:16,]

colnames(evanno)=c("K","Reps","Mean_LnP","sd_LnP","LnK","Ln2K","DeltaK")

######## Mean Ln(P) as a function of K
(pLnP=ggplot(data=evanno, aes(x=K, y=Mean_LnP)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    geom_errorbar(aes(ymin=Mean_LnP-sd_LnP, ymax=Mean_LnP+sd_LnP), width=.2) +
    xlab("K") + ylab("Ln(P)") +
    scale_x_discrete(limits=c(2:nrow(evanno))) +
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
          legend.title=element_text(size=14)))


######## Mean Ln(K) as a function of K
(pLnK=ggplot(data=evanno, aes(x=K, y=LnK)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    xlab("K") + ylab("Ln'(K)") +
    scale_x_discrete(limits=c(2:nrow(evanno))) +
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
          legend.title=element_text(size=14)))


######## Mean Ln''(K) as a function of K
(pLn2K=ggplot(data=evanno, aes(x=K, y=Ln2K)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    xlab("K") + ylab("Ln''(K)") +
    scale_x_discrete(limits=c(2:nrow(evanno))) +
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
          legend.title=element_text(size=14)))


######## Mean Ln''(K) as a function of K
(pDeltaK=ggplot(data=evanno, aes(x=K, y=DeltaK)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    xlab("K") + ylab("Delta(K)") +
    scale_x_discrete(limits=c(2:nrow(evanno))) +
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
          legend.title=element_text(size=14)))

#### SAVE TO FIGURE
ggarrange(pLnP,pLnK,pLn2K,pDeltaK,widths=1:1,heights=1:1,labels="auto")
# ggsave(paste(figuresdir,"/STRUCTURE/2a_Native/Diagnostic_plots_STRUCTURE_native.png",sep=""),
#        device="png",dpi=320,units="cm",width=40,height=32)
# ggsave(paste(figuresdir,"/STRUCTURE/2b_Native/Diagnostic_plots_STRUCTURE_native.png",sep=""),
#        device="png",dpi=320,units="cm",width=40,height=32)
ggsave(paste(figuresdir,"/STRUCTURE/2b_Native_OC/Diagnostic_plots_STRUCTURE_native.png",sep=""),
       device="png",dpi=320,units="cm",width=40,height=32)
# Based on convergent replicates (beware, low number, 3 replicates), the delta K method suggested multiple levels of genetic structure (i.e. multiple values of most probable K)
# K=3, 7, 9, 13, 15 and 17 evidences that
# a/ the population could be hierarchically structured
# or
# b/ there is no clear genetic structure, but a large metapopulation with gene flows between sampled populations
# We must look for K= 3, 7, 9, 13, 15 and 17


# 2/ Quantitative evaluation with the G' (or H') parameter computed by CLUMPP to assess convergence of STRUCTURE



### Assess the convergence for parameter estimates ###
# See STRUCTURE.R for convergence assessment...

### Convergence of Q


### Convergence of Alpha (admixture parameter) per individual

#-------------------------------------
### Bar plot representation ###
# Representation of admixture (Q) for each individual FOR K=6 (chosen K)
# Format a data frame with 4 columns, and as many rows per individual as there is clusters
# strK6=read.table(paste(STRdir,"/2a_Native/CLUMPP/K6.outfile",sep=""))
strK6=read.table(paste(STRdir,"/2b_Native/CLUMPP/K6.outfile",sep=""))
# Set a vector of colors for clusters
cols.native=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
              "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9")

# Get the same colors and cluster number as DAPC
strK6.results=as.data.frame(strK6[,(c(1,2,3,5,4,6)+5)]) # population clusters in the good order for colors
colnames(strK6.results)=c(1,2,3,4,5,6)
# strK6.results$ind=as.character(strK6[,1])
strK6.results$pop=as.character(strK6[,4])
for (i in 1:nrow(strK6.results)) {
  strK6.results$pop[i]=as.character(labelPops[which(strK6.results$pop[i] == labelPops$Location_index),1])
}
strK6.results$indNames=names(Native$tab[,1])
strK6.results=melt(strK6.results)

colnames(strK6.results)=c("Population", "Individual","Cluster","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK6.results$Population=factor(strK6.results$Population,
                                         levels = c("S9","S18","S19","S20","Tib","S10","S11",
                                                    "S4","S6",
                                                    "S3","S1","S2","S16","S13","S14",
                                                    "S15","S17", "Jap"))


write.table(strK6.results,"Tables/Figure_STRUCTURE_NativeK6.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)


# Plotting the barplot
strK6_barplot=ggplot(data=strK6.results, aes(x=Population, y=Posterior_membership_probability, fill=Cluster))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.native) +
  facet_grid(~Population, scales = "free",space="free_x") +
  labs(x="Sampled individual", y="Posterior\nadmixture\nproportions\n", fill="Assigned\npopulation") +
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
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=22, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
strK6_barplot
ggsave(paste(figuresdir,"/STRUCTURE/2b_Native/Native_barplot_K6b.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)


#############################
############################
# What if K=3, evidences that K=3 was a probable value
# Representation of admixture (Q) for each individual FOR K=3 (chosen K)
# Format a data frame with 3 columns, and as many rows per individual as there is clusters
# strK3=read.table(paste(STRdir,"/2a_Native/CLUMPP/K3.outfile",sep=""))
strK3=read.table(paste(STRdir,"/2b_Native/CLUMPP/K3.outfile",sep=""))
# Get the same colors and cluster number as DAPC
strK3.results=as.data.frame(strK3[,(c(1,2,3)+5)]) # population clusters in the good order for colors

colnames(strK3.results)=c(1,2,3)
strK3.results$pop=as.character(strK3[,4])
for (i in 1:nrow(strK3.results)) {
  strK3.results$pop[i]=as.character(labelPops[which(strK3.results$pop[i] == labelPops$Location_index),1])
}
strK3.results$indNames=names(Native$tab[,1])
strK3.results=melt(strK3.results)

colnames(strK3.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK3.results$Original_Population=factor(strK3.results$Original_Population,
                                         levels = c("S9","S18","S19","S20","S4","S6",
                                                    "S3","S1","S2","S16","S13","S14",
                                                    "S15","S17","Tib","S11","S10", "Jap"))

# Plotting the barplot
strK3_barplot=ggplot(data=strK3.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.native) +
  facet_grid(~Original_Population, scales = "free",space = "free_x") +
  labs(x="Sampled individual", y="Posterior\nadmixture\nproportions\n", fill="Assigned\npopulation") +
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
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=22, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
strK3_barplot


ggsave(paste(figuresdir,"/STRUCTURE/2b_Native/Native_barplot_K3b.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)



#----------------------------------------------------------
# Fig 3. Map of the native area with pies chart of admixture proportions on each sampling location

# Prepare the data set for mapping
Native.pies.str=cbind(as.character(labelPops$Location_name[strK6[,4]]),as.data.frame(strK6[,(c(1,2,3,5,4,6)+5)])) # Clusters are selected to correspodn to DAPC colors
names(Native.pies.str)=c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")

Native.pie.str_means=sapply(split(Native.pies.str[2:7],Native.pies.str$Population), colMeans) ## This just makes population cluster averages for each cluster

write.table(Native.pie.str_means,"Tables/PieCharts_STRUCTURE_NativeK6.txt",sep="\t",quote=FALSE,row.names=FALSE,col.names=TRUE)

# Mapping
# Download database of river network at scale 50
# rivers50=ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')

# png("Figures/STRUCTURE/2a_Native/STR_Native_pies_map.png",width=2000,height=1100)
png("Figures/STRUCTURE/2b_Native/STR_Native_pies_map_K6.png",width=2000,height=1100)
map("worldHires", xlim=c(92, 160), ylim=c(20,48), col="gray90",  fill=TRUE) # Plot the map of Asian area
sp::plot(rivers50,lwd=2, col = 'blue',add=TRUE)
# Nomenclature
north.arrow(155,45,1, lab="N", lab.pos="above")
scalebar(c(150,20),1000, bg=NA, border=NA, division.cex=2)

text(x=105, y=34,"China",cex=3,font=2)
text(x=139, y=39,"Japan",cex=3,font=2)
points(native.coords$X, native.coords$Y, pch = 16, cex = 0.7,col="red") #\ Add points to the map
# map.scale(x=140,y=25,relwidth=0.1,cex=2)

## Pies ## 
# Draw a pie with admixture proportions for each population
for (i in 1:ncol(Native.pie.str_means)) {
  add.pie(Native.pie.str_means[,i],x=native.coords$X[i], y=native.coords$Y[i],labels=native.coords[i,1],label.dist=1.5,
          radius=1,edges=200,clockwise=T,	col	=	cols.native,cex=2,font=2)
  # Add labels to pies
}
dev.off()


#############################
############################
# What if K=7
# Representation of admixture (Q) for each individual FOR K=7 (chosen K)
# Format a data frame with 7 columns, and as many rows per individual as there is clusters
strK7=read.table(paste(STRdir,"/2b_Native/CLUMPP/K7.outfile",sep=""))
# Get the same colors and cluster number as DAPC
strK7.results=as.data.frame(strK7[,(c(2,4,7,3,1,6,5)+5)]) # population clusters in the good order for colors

colnames(strK7.results)=c(1,2,3,4,5,6,7)
strK7.results$pop=as.character(strK7[,4])
for (i in 1:nrow(strK7.results)) {
  strK7.results$pop[i]=as.character(labelPops[which(strK7.results$pop[i] == labelPops$Location_index),1])
}
strK7.results$indNames=names(Native$tab[,1])
strK7.results=melt(strK7.results)

colnames(strK7.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK7.results$Original_Population=factor(strK7.results$Original_Population, levels = c("Jap","S4","S6","S3","S1","S2","S16","S13","S14","S15","S17","S9","S18","S19","S20","S11","S10","Tib"))

# Plotting the barplot
strK7_barplot=ggplot(data=strK7.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.native) +
  facet_grid(~Original_Population, scales = "free", space = "free_x") +
  labs(x="Sampled individual", y="Posterior\nadmixture\nproportions\n", fill="Assigned\npopulation") +
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
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=22, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
strK7_barplot

ggsave(paste(figuresdir,"/STRUCTURE/2b_Native/Native_barplot_K7b.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)
# ggsave(paste(figuresdir,"/STRUCTURE/2b_Native/Native_barplot_K7b.jpg",sep=""),
# device="jpg",dpi=150,units="cm",width=60,height=9)



############################
# What if K=8
strK8=read.table(paste(STRdir,"/2b_Native/CLUMPP/K8.outfile",sep=""))
# Get the same colors and cluster number as DAPC
strK8.results=as.data.frame(strK8[,(c(2,4,7,3,1,6,5,8)+5)]) # population clusters in the good order for colors

colnames(strK8.results)=c(1,2,3,4,5,6,7,8)
strK8.results$pop=as.character(strK8[,4])
for (i in 1:nrow(strK8.results)) {
  strK8.results$pop[i]=as.character(labelPops[which(strK8.results$pop[i] == labelPops$Location_index),1])
}
strK8.results$indNames=names(Native$tab[,1])
strK8.results=melt(strK8.results)

colnames(strK8.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK8.results$Original_Population=factor(strK8.results$Original_Population, levels = c("Jap","S4","S6","S3","S1","S2","S16","S13","S14","S15","S17","S9","S18","S19","S20","S11","S10","Tib"))

# Plotting the barplot
strK8_barplot=ggplot(data=strK8.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.native) +
  facet_grid(~Original_Population, scales = "free", space = "free_x") +
  labs(x="Sampled individual", y="Posterior\nadmixture\nproportions\n", fill="Assigned\npopulation") +
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
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=22, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
strK8_barplot

ggsave(paste(figuresdir,"/STRUCTURE/2b_Native/Native_barplot_K8b.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)


############################
# What if K=9
strK9=read.table(paste(STRdir,"/2b_Native/CLUMPP/K9.outfile",sep=""))
# Get the same colors and cluster number as DAPC
strK9.results=as.data.frame(strK9[,(c(2,4,7,3,1,6,5,8,9)+5)]) # population clusters in the good order for colors

colnames(strK9.results)=c(1,2,3,4,5,6,7,8,9)
strK9.results$pop=as.character(strK9[,4])
for (i in 1:nrow(strK9.results)) {
  strK9.results$pop[i]=as.character(labelPops[which(strK9.results$pop[i] == labelPops$Location_index),1])
}
strK9.results$indNames=names(Native$tab[,1])
strK9.results=melt(strK9.results)

colnames(strK9.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK9.results$Original_Population=factor(strK9.results$Original_Population, levels = c("Jap","S4","S6","S3","S1","S2","S16","S13","S14","S15","S17","S9","S18","S19","S20","S11","S10","Tib"))

# Plotting the barplot
strK9_barplot=ggplot(data=strK9.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.native) +
  facet_grid(~Original_Population, scales = "free", space = "free_x") +
  labs(x="Sampled individual", y="Posterior\nadmixture\nproportions\n", fill="Assigned\npopulation") +
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
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=22, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
strK9_barplot

ggsave(paste(figuresdir,"/STRUCTURE/2b_Native/Native_barplot_K9b.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)


############################
# What if K=10
strK10=read.table(paste(STRdir,"/2b_Native/CLUMPP/K10.outfile",sep=""))
# Get the same colors and cluster number as DAPC
strK10.results=as.data.frame(strK10[,(c(2,4,7,3,1,6,5,8,9,10)+5)]) # population clusters in the good order for colors

colnames(strK10.results)=c(1,2,3,4,5,6,7,8,9,10)
strK10.results$pop=as.character(strK10[,4])
for (i in 1:nrow(strK10.results)) {
  strK10.results$pop[i]=as.character(labelPops[which(strK10.results$pop[i] == labelPops$Location_index),1])
}
strK10.results$indNames=names(Native$tab[,1])
strK10.results=melt(strK10.results)

colnames(strK10.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK10.results$Original_Population=factor(strK10.results$Original_Population, levels = c("Jap","S4","S6","S3","S1","S2","S16","S13","S14","S15","S17","S9","S18","S19","S20","S11","S10","Tib"))

# Plotting the barplot
strK10_barplot=ggplot(data=strK10.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.native) +
  facet_grid(~Original_Population, scales = "free", space = "free_x") +
  labs(x="Sampled individual", y="Posterior\nadmixture\nproportions\n", fill="Assigned\npopulation") +
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
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=22, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
strK10_barplot

ggsave(paste(figuresdir,"/STRUCTURE/2b_Native/Native_barplot_K10b.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)


#----------------------------------------------------------
# Fig 3. Map of the native area with pies chart of admixture proportions on each sampling location
#----------------------------------------------------------
# Prepare the data set for mapping
Native.pies.str=cbind(as.character(labelPops$Location_name[strK7[,4]]),as.data.frame(strK7[,(c(2,4,7,3,1,6,5)+5)])) # Clusters are selected to correspodn to DAPC colors
names(Native.pies.str)=c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6", "Cluster 7")

Native.pie.str_means=sapply(split(Native.pies.str[2:8],Native.pies.str$Population), colMeans) ## This just makes population cluster averages for each cluster

# Mapping
# Download database of river network at scale 50
# rivers50=ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')

# png("Figures/STRUCTURE/2a_Native/STR_Native_pies_map.png",width=2000,height=1100)
png("Figures/STRUCTURE/2b_Native/STR_Native_pies_map_K7.png",width=2000,height=1100)
map("worldHires", xlim=c(92, 160), ylim=c(20,48), col="gray90",  fill=TRUE) # Plot the map of Asian area
sp::plot(rivers50,lwd=2, col = 'blue',add=TRUE)
# Nomenclature
north.arrow(155,45,1, lab="N", lab.pos="above")
scalebar(c(150,20),1000, bg=NA, border=NA, division.cex=2)

text(x=105, y=34,"China",cex=3,font=2)
text(x=139, y=39,"Japan",cex=3,font=2)
points(native.coords$X, native.coords$Y, pch = 16, cex = 0.7,col="red") #\ Add points to the map
# map.scale(x=140,y=25,relwidth=0.1,cex=2)

## Pies ## 
# Draw a pie with admixture proportions for each population
for (i in 1:ncol(Native.pie.str_means)) {
  add.pie(Native.pie.str_means[,i],x=native.coords$X[i], y=native.coords$Y[i],labels=native.coords[i,1],label.dist=1.5,
          radius=1,edges=200,clockwise=T,	col	=	cols.native,cex=2,font=2)
  # Add labels to pies
}
dev.off()




#----------------------------------------------------------
# Admixed populations: assess which population have less than 0.7 of mean Q (ancestry coefficient)
#----------------------------------------------------------
# A population where the mean of the max Q (the major ancestry coeff) over all individuals is inferior to 0.7 should be considered as an admixed population
# Get max Q for each individuals, and then mean(max Q) for each population, for the chosen K: here K=7
# Qk is the proportion of individual's ancestry from population k
ancestry=strK7[,c(1,4,6:12)]
colnames(ancestry)=c("ind","pop","cluster1","cluster2","cluster3","cluster4","cluster5","cluster6","cluster7")
ancestry$Qmax=apply(ancestry[,3:9],1,max)

# Get the mean (max Q)
pop.admixed=data.frame(pop=unique(ancestry$pop))
i=0
for (p in unique(ancestry$pop)) {
  print(p)
  i=i+1
  pop.admixed$Q[i]=mean(ancestry$Qmax[which(ancestry$pop==p)])
  pop.admixed$pop[i]=as.character(labelPops$Location_name[which(labelPops$Location_index==p)])
}
pop.admixed$Q=round(pop.admixed$Q,digits = 3)
write.table(x=pop.admixed,file="Tables/Major ancestry coeff Native K7.txt",sep="\t",quote = FALSE,row.names = FALSE)




#########################
# HIERARCHICAL ANALYSIS
# STRUCTURE on a subset of populations
# Since some populations showed a high elvel of admixture, we performed a subsequent analysis only on these subset to try to find a sub-structure
# Hence, populations that clustered well previously for K=3 were not considered
# 2b1_Native: 2,112 loci, burnin 100,000 + 100,000 iterations, subset of 13 populations in native area that showed admixture



# 1/ Graphic evaluation with the four plots of Evanno's method (Evanno et al. 2005)
# Most important plot is DeltaK = mean(|L''(K)|)/sd(L(K)) of L(K) as a function of K 

### Choose the most probable K value (number of cluster inferred), K value the most appropriate for the data
# Load results from Harvester: 'evanno.txt'
# evanno1=read.table(paste(STRdir,"/2b1_Native/Harvester/evanno.txt",sep=""),header=F,sep="\t")
evanno1=read.table(paste(STRdir,"/2b1_Native_OC/Harvester/evanno.txt",sep=""),header=F,sep="\t")

colnames(evanno1)=c("K","Reps","Mean_LnP","sd_LnP","LnK","Ln2K","DeltaK")

######## Mean Ln(P) as a function of K
(pLnP=ggplot(data=evanno1, aes(x=K, y=Mean_LnP)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    geom_errorbar(aes(ymin=Mean_LnP-sd_LnP, ymax=Mean_LnP+sd_LnP), width=.2) +
    xlab("K") + ylab("Ln(P)") +
    scale_x_discrete(limits=c(2:nrow(evanno))) +
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
          legend.title=element_text(size=14)))


######## Mean Ln(K) as a function of K
(pLnK=ggplot(data=evanno1, aes(x=K, y=LnK)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    xlab("K") + ylab("Ln'(K)") +
    scale_x_discrete(limits=c(2:nrow(evanno))) +
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
          legend.title=element_text(size=14)))


######## Mean Ln''(K) as a function of K
(pLn2K=ggplot(data=evanno1, aes(x=K, y=Ln2K)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    xlab("K") + ylab("Ln''(K)") +
    scale_x_discrete(limits=c(2:nrow(evanno))) +
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
          legend.title=element_text(size=14)))


######## Mean Ln''(K) as a function of K
(pDeltaK=ggplot(data=evanno1, aes(x=K, y=DeltaK)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    xlab("K") + ylab("Delta(K)") +
    scale_x_discrete(limits=c(2:nrow(evanno))) +
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
          legend.title=element_text(size=14)))

#### SAVE TO FIGURE
ggarrange(pLnP,pLnK,pLn2K,pDeltaK,widths=1:1,heights=1:1,labels="auto")
# ggsave(paste(figuresdir,"/STRUCTURE/2b1_Native/Diagnostic_plots_STRUCTURE_native.png",sep=""),
#        device="png",dpi=320,units="cm",width=40,height=32)
ggsave(paste(figuresdir,"/STRUCTURE/2b1_Native_OC/Diagnostic_plots_STRUCTURE_native.png",sep=""),
       device="png",dpi=320,units="cm",width=40,height=32)

# This hierarchical approach (subset of admiwed populations) did not revealed furthermore a cryptic structure of populations in the admixed region. K=2 with Delta K approach







#----------------------------------------------------------
#         INVASIVE POPULATIONS IN EUROPE
#----------------------------------------------------------
evannoEu=read.table(paste(STRdir,"/3a_Invasive/Harvester/evanno.txt",sep=""),header=F,sep="\t")
# evannoEu=read.table(paste(STRdir,"/3a_Invasive_OC/Harvester/evanno.txt",sep=""),header=F,sep="\t")

colnames(evannoEu)=c("K","Reps","Mean_LnP","sd_LnP","LnK","Ln2K","DeltaK")

######## Mean Ln(P) as a function of K
(pLnP=ggplot(data=evannoEu, aes(x=K, y=Mean_LnP)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    geom_errorbar(aes(ymin=Mean_LnP-sd_LnP, ymax=Mean_LnP+sd_LnP), width=.2) +
    xlab("K") + ylab("Ln(P)") +
    scale_x_discrete(limits=c(2:nrow(evannoEu))) +
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
          legend.title=element_text(size=14)))


######## Mean Ln(K) as a function of K
(pLnK=ggplot(data=evannoEu, aes(x=K, y=LnK)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    xlab("K") + ylab("Ln'(K)") +
    scale_x_discrete(limits=c(2:nrow(evannoEu))) +
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
          legend.title=element_text(size=14)))


######## Mean Ln''(K) as a function of K
(pLn2K=ggplot(data=evannoEu, aes(x=K, y=Ln2K)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    xlab("K") + ylab("Ln''(K)") +
    scale_x_discrete(limits=c(2:nrow(evannoEu))) +
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
          legend.title=element_text(size=14)))


######## Mean Ln''(K) as a function of K
(pDeltaK=ggplot(data=evannoEu, aes(x=K, y=DeltaK)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    xlab("K") + ylab("Delta(K)") +
    scale_x_discrete(limits=c(2:nrow(evannoEu))) +
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
          legend.title=element_text(size=14)))

#### SAVE TO FIGURE
ggarrange(pLnP,pLnK,pLn2K,pDeltaK,widths=1:1,heights=1:1,labels="auto")
ggsave(paste(figuresdir,"/STRUCTURE/3a_Invasive/Diagnostic_plots_STRUCTURE_native.png",sep=""),
       device="png",dpi=320,units="cm",width=40,height=32)

# 2/ Quantitative evaluation with the G' (or H') parameter computed by CLUMPP to assess convergence of STRUCTURE



### Assess the convergence for parameter estimates ###
# See STRUCTURE.R for convergence assessment...

### Convergence of Q


### Convergence of Alpha (admixture parameter) per individual

#-------------------------------------
# Set a vector of colors for clusters
cols.invasive=brewer.pal(n=9,name="Pastel1")

### Bar plot representation for K=2 ###
# Representation of admixture (Q) for each individual FOR K=2 (chosen Delta K)
# Format a data frame with 4 columns, and as many rows per individual as there is clusters
strK2Eu=read.table(paste(STRdir,"/3a_Invasive/CLUMPP/K2.outfile",sep=""))
# Get the same colors and cluster number as DAPC
strK2Eu.results=as.data.frame(strK2Eu[,(c(1,2)+5)]) # population clusters in the good order for colors
colnames(strK2Eu.results)=c(1,2)
strK2Eu.results$pop=as.character(strK2Eu[,4])
for (i in 1:nrow(strK2Eu.results)) {
  strK2Eu.results$pop[i]=as.character(labelPops[which(strK2Eu.results$pop[i] == labelPops$Location_index),1])
}
strK2Eu.results$indNames=names(invasive$tab[,1])
strK2Eu.results=melt(strK2Eu.results)
colnames(strK2Eu.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK2Eu.results$Original_Population=factor(strK2Eu.results$Original_Population, levels = c("Ira","Aus","Bel","UK","Hun","Ita",
                                                                                           "Pol","Spa","Tur","Bul1","Bul2"))
# Plotting the barplot
strK2Eu_barplot=ggplot(data=strK2Eu.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.invasive) +
  facet_grid(~Original_Population, scales = "free",space="free_x") +
  labs(x="Sampled individual", y="Posterior\nadmixture\nproportions\n", fill="Assigned\npopulation") +
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
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=22, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
strK2Eu_barplot

ggsave(paste(figuresdir,"/STRUCTURE/3a_Invasive/Invasive_barplot_K2.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)
# ggsave(paste(figuresdir,"/STRUCTURE/3a_Invasive/Invasive_barplot_K2.jpg",sep=""),
#        device="jpg",dpi=150,units="cm",width=60,height=9)



#-------------------------------------
### Bar plot representation for K=3 ###
# Representation of admixture (Q) for each individual FOR K=3 (chosen Delta K)
# Format a data frame with 4 columns, and as many rows per individual as there is clusters
cols.invasive = c("#ABD9E9", "#FEC44F", "#1B5E20", "#FEC44F")
strK3Eu=read.table(paste(STRdir,"/3a_Invasive/CLUMPP/K3.outfile",sep=""))
# Get the same colors and cluster number as DAPC
strK3Eu.results=as.data.frame(strK3Eu[,(c(1,2,3)+5)]) # population clusters in the good order for colors
colnames(strK3Eu.results)=c(1,2,3)
strK3Eu.results$pop=as.character(strK3Eu[,4])
for (i in 1:nrow(strK3Eu.results)) {
  strK3Eu.results$pop[i]=as.character(labelPops[which(strK3Eu.results$pop[i] == labelPops$Location_index),1])
}
strK3Eu.results$indNames=names(Invasive$tab[,1])
strK3Eu.results=melt(strK3Eu.results)
colnames(strK3Eu.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK3Eu.results$Original_Population=factor(strK3Eu.results$Original_Population, levels = c("Aus","Bel","UK","Hun","Ita",
                                                                                           "Pol","Spa","Tur","Bul1","Bul2", "Ira"))
# Plotting the barplot
strK3Eu_barplot=ggplot(data=strK3Eu.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.invasive) +
  facet_grid(~Original_Population, scales = "free",space="free_x") +
  labs(x="Sampled individual", y="Posterior\nadmixture\nproportions\n", fill="Assigned\npopulation") +
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
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=22, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
strK3Eu_barplot

ggsave(paste(figuresdir,"/STRUCTURE/3a_Invasive/Invasive_barplot_K3.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)


#-------------------------------------
### Bar plot representation for K=2 ###
# Representation of admixture (Q) for each individual FOR K=2 (chosen Delta K)
# Format a data frame with 4 columns, and as many rows per individual as there is clusters
cols.invasive = c("#ABD9E9", "#993404", "#A6DBA0", "#FEC44F")
strK4Eu=read.table(paste(STRdir,"/3a_Invasive/CLUMPP/K4.outfile",sep=""))
# Get the same colors and cluster number as DAPC
strK4Eu.results=as.data.frame(strK4Eu[,(c(1,2,3,4)+5)]) # population clusters in the good order for colors
colnames(strK4Eu.results)=c(1,2,3,4)
strK4Eu.results$pop=as.character(strK4Eu[,4])
for (i in 1:nrow(strK4Eu.results)) {
  strK4Eu.results$pop[i]=as.character(labelPops[which(strK4Eu.results$pop[i] == labelPops$Location_index),1])
}
strK4Eu.results$indNames=names(Invasive$tab[,1])
strK4Eu.results=melt(strK4Eu.results)
colnames(strK4Eu.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK4Eu.results$Original_Population=factor(strK4Eu.results$Original_Population, levels = c("Aus","Bel","UK","Hun","Ita",
                                                                                           "Pol","Spa","Tur","Bul1","Bul2", "Ira"))
# Plotting the barplot
strK4Eu_barplot=ggplot(data=strK4Eu.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.invasive) +
  facet_grid(~Original_Population, scales = "free",space="free_x") +
  labs(x="Sampled individual", y="Posterior\nadmixture\nproportions\n", fill="Assigned\npopulation") +
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
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=22, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
strK4Eu_barplot

ggsave(paste(figuresdir,"/STRUCTURE/3a_Invasive/Invasive_barplot_K4.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)



#-------------------------------------
### Bar plot representation for K=9 ###
# Representation of admixture (Q) for each individual FOR K=9
# Format a data frame with 4 columns, and as many rows per individual as there is clusters

# New color palette
cols.invasive = c("#A6DBA0", "#D9F0D3",  "#ABD9E9", "#993404", "#C2A5CF",
                  "#E0F3F8", "#FE9929", "#FEC44F", "#762A83")
strK9Eu=read.table(paste(STRdir,"/3a_Invasive/CLUMPP/K9.outfile",sep=""))
# Get the same colors and cluster number as DAPC
strK9Eu.results=as.data.frame(strK9Eu[,(c(1,2,3,4,5,6,7,8,9)+5)]) # population clusters in the good order for colors
colnames(strK9Eu.results)=c(1,2,3,4,5,6,7,8,9)
strK9Eu.results$pop=as.character(strK9Eu[,4])
for (i in 1:nrow(strK9Eu.results)) {
  strK9Eu.results$pop[i]=as.character(labelPops[which(strK9Eu.results$pop[i] == labelPops$Location_index),1])
}
strK9Eu.results$indNames=names(Invasive$tab[,1])
strK9Eu.results=melt(strK9Eu.results)
colnames(strK9Eu.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK9Eu.results$Original_Population=factor(strK9Eu.results$Original_Population, levels = c("Aus","Bel","UK","Hun","Ita",
                                                                                           "Pol","Spa","Tur","Bul1","Bul2", "Ira"))
# Plotting the barplot
strK9Eu_barplot=ggplot(data=strK9Eu.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.invasive) +
  facet_grid(~Original_Population, scales = "free", space = "free_x") +
  labs(x="Sampled individual", y="Posterior\nadmixture\nproportions\n", fill="Assigned\npopulation") +
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
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=22, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=22),
        legend.title=element_text(size=22))
strK9Eu_barplot

ggsave(paste(figuresdir,"/STRUCTURE/3a_Invasive/Invasive_barplot_K9_newcols.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)

#----------------------------------------------------------
# Fig 3. Map of the invasive area with pies chart of admixture proportions on each sampling location

# Prepare the data set for mapping
Invasive.pies.str=cbind(as.character(labelPops$Location_name[strK9Eu[,4]]),as.data.frame(strK9Eu[,(c(1,2,3,4,5,6,7,8,9)+5)])) # Clusters are selected to correspodn to DAPC colors
names(Invasive.pies.str)=c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4",
                           "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9")

Invasive.pies.str_means=sapply(split(Invasive.pies.str[2:10],Invasive.pies.str$Population), colMeans) ## This just makes population cluster averages for each cluster

# Mapping
# Download database of river network at scale 50
rivers50=ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')

png("Figures/STRUCTURE/3a_Invasive/STR_Invasive_pies_map_K9.png",width=2000,height=1100)
map("worldHires", xlim=c(-15,60), ylim=c(30,59), col="gray90",  fill=TRUE) # Plot the map of European area
map(database="worldHires", regions="Caspian Sea", col="white", fill=TRUE, add=TRUE) # Correct for unfilled Caspian sea
map(database="worldHires", regions="Aral Sea", col="white", fill=TRUE, add=TRUE)# Correct for unfilled Aral sea
sp::plot(rivers50,xlim=c(-15,30), ylim=c(30,59),lwd=2, col = 'blue',add=TRUE)
# Nomenclature
north.arrow(56,56,1, lab="N", lab.pos="above")
scalebar(c(49,31),1000, bg=NA, border=NA, division.cex=2)


text(x=105, y=34,"",cex=3,font=2)
text(x=139, y=39,"",cex=3,font=2)
points(invasive.coords$X, invasive.coords$Y, pch = 16, cex = 0.7,col="red") #\ Add points to the map

## Pies ## 
# Draw a pie with admixture proportions for each population
for (i in 1:ncol(Invasive.pies.str_means)) {
  add.pie(Invasive.pies.str_means[,i],x=invasive.coords$X[i], y=invasive.coords$Y[i],labels=invasive.coords[i,1],label.dist=1.5,
          radius=1,edges=200,clockwise=T,	col	=	cols.invasive,cex=2,font=2)
  # Add labels to pies
}
dev.off()




#################
# HIERARCHICAL PROCEDURE
# STRUCTURE on a subset of populations
# Since some populations showed a high elvel of admixture, we performed a subsequent analysis only on these subset to try to find a sub-structure
# Hence, populations that clustered well previously for K=3 were not considered
# 3a1_Native: 2,112 loci, burnin 100,000 + 100,000 iterations, subset of populations in invasive area that showed admixture

# evannoEu1=read.table(paste(STRdir,"/3a1_Invasive/Harvester/evanno.txt",sep=""),header=F,sep="\t")
evannoEu1=read.table(paste(STRdir,"/3a1_Invasive_OC/Harvester/evanno.txt",sep=""),header=F,sep="\t")

colnames(evannoEu1)=c("K","Reps","Mean_LnP","sd_LnP","LnK","Ln2K","DeltaK")

######## Mean Ln(P) as a function of K
(pLnP=ggplot(data=evannoEu1, aes(x=K, y=Mean_LnP)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    geom_errorbar(aes(ymin=Mean_LnP-sd_LnP, ymax=Mean_LnP+sd_LnP), width=.2) +
    xlab("K") + ylab("Ln(P)") +
    scale_x_discrete(limits=c(2:nrow(evannoEu1))) +
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
          legend.title=element_text(size=14)))


######## Mean Ln(K) as a function of K
(pLnK=ggplot(data=evannoEu1, aes(x=K, y=LnK)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    xlab("K") + ylab("Ln'(K)") +
    scale_x_discrete(limits=c(2:nrow(evannoEu1))) +
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
          legend.title=element_text(size=14)))


######## Mean Ln''(K) as a function of K
(pLn2K=ggplot(data=evannoEu1, aes(x=K, y=Ln2K)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    xlab("K") + ylab("Ln''(K)") +
    scale_x_discrete(limits=c(2:nrow(evannoEu1))) +
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
          legend.title=element_text(size=14)))


######## Mean Ln''(K) as a function of K
(pDeltaK=ggplot(data=evannoEu1, aes(x=K, y=DeltaK)) +
    geom_point(colour="Black",size=1)+
    geom_line() +
    xlab("K") + ylab("Delta(K)") +
    scale_x_discrete(limits=c(2:nrow(evannoEu1))) +
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
          legend.title=element_text(size=14)))

#### SAVE TO FIGURE
ggarrange(pLnP,pLnK,pLn2K,pDeltaK,widths=1:1,heights=1:1,labels="auto")
ggsave(paste(figuresdir,"/STRUCTURE/3a1_Invasive/Diagnostic_plots_STRUCTURE_native.png",sep=""),
       device="png",dpi=320,units="cm",width=40,height=32)
# ggsave(paste(figuresdir,"/STRUCTURE/3a1_Invasive_OC/Diagnostic_plots_STRUCTURE_native.png",sep=""),
#        device="png",dpi=320,units="cm",width=40,height=32)

# This new hierarchical analysis on a subset of Western populations showed a genetic structure with a most probable K of 3 or 4 (most probably 3 for DeltaK)
# Hence, we displayed the clustering for K=3 and K=4 in western Europe
# Selecting only the convergent replicates did not changed the inferenc of K --> K= 3 or 4
# We showed that in Western Europe, populations were close genetically, but there were nonetheless an underlying genetic structure
# A hierarchical genetic structure of the western Europe could indicate a stepping stone model for range expansion with successive founding events and bottlenecks,
# OR
# multiple introductions from the same source populations but with stochasticity of sampling of invading individuals
# See Thibault et al. 2009

#-------------------------------------
### Bar plot representation for K=5 ###
# Representation of admixture (Q) for each individual FOR K=3
# Format a data frame with 4 columns, and as many rows per individual as there is clusters
strK5Eu1=read.table(paste(STRdir,"/3a1_Invasive/CLUMPP/K5.outfile",sep=""))
# Set a vector of colors for clusters
# cols.invasive=brewer.pal(n=8,name="Dark2")
cols.invasive = c("#993404", "#762A83",  "#ABD9E9", "#A6DBA0", "#C2A5CF")
# cols.invasive = c("#A6DBA0", "#D9F0D3",  "#ABD9E9", "#993404", "#C2A5CF",
#                   "#E0F3F8", "#FE9929", "#FEC44F", "#762A83")
# Get the same colors and cluster number as DAPC
strK5Eu1.results=as.data.frame(strK5Eu1[,(c(1,2,3,4,5)+5)]) # population clusters in the good order for colors
colnames(strK5Eu1.results)=c(1,2,3,4,5)
strK5Eu1.results$pop=as.character(strK5Eu1[,4])
for (i in 1:nrow(strK5Eu1.results)) {
  strK5Eu1.results$pop[i]=as.character(labelPops[which(strK5Eu1.results$pop[i] == labelPops$Location_index),1])
}
strK5Eu1.results$indNames=names(Invasive$tab[which(Invasive@pop %in% c("Aus","Bel","Hun","Ita","Pol","Spa","UK")),1])
strK5Eu1.results=melt(strK5Eu1.results)
colnames(strK5Eu1.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK5Eu1.results$Original_Population=factor(strK5Eu1.results$Original_Population, levels = c("Aus","Bel","UK","Hun","Ita",
                                                                                             "Pol","Spa"))
# Plotting the barplot
strK5Eu1_barplot=ggplot(data=strK5Eu1.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.invasive) +
  facet_grid(~Original_Population, scales = "free", space = "free_x") +
  labs(x="Sampled individual", y="Posterior\nadmixture\nproportions\n", fill="Assigned\npopulation") +
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
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=22, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
strK5Eu1_barplot

ggsave(paste(figuresdir,"/STRUCTURE/3a1_Invasive/Invasive_barplot_K5.png",sep=""),
       device="png",dpi=320,units="cm",width=60,height=9)


#----------------------------------------------------------
# Admixed populations: assess which population have less than 0.7 of mean Q (ancestry coefficient)
#----------------------------------------------------------
# A population where the mean of the max Q (the major ancestry coeff) over all individuals is inferior to 0.7 should be considered as an admixed population
# Get max Q for each individuals, and then mean(max Q) for each population, for the chosen K: here K=7
# Qk is the proportion of individual's ancestry from population k

# Hierarchical structure: all Europe
ancestry=strK4Eu[,c(1,4,6:9)]
colnames(ancestry)=c("ind","pop","cluster1","cluster2","cluster3","Cluster4")
ancestry$Qmax=apply(ancestry[,3:6],1,max)

# Get the mean (max Q)
pop.admixed=data.frame(pop=unique(ancestry$pop))
i=0
for (p in unique(ancestry$pop)) {
  print(p)
  i=i+1
  pop.admixed$Q[i]=mean(ancestry$Qmax[which(ancestry$pop==p)])
  pop.admixed$pop[i]=as.character(labelPops$Location_name[which(labelPops$Location_index==p)])
}
pop.admixed$Q=round(pop.admixed$Q,digits = 3)
write.table(x=pop.admixed,file="Tables/Major ancestry coeff Invasive K4.txt",sep="\t",quote = FALSE,row.names = FALSE)


# Hierarchical structure: only Western Europe
ancestry=strK3Eu1[,c(1,4,6:8)]
colnames(ancestry)=c("ind","pop","cluster1","cluster2","cluster3")
ancestry$Qmax=apply(ancestry[,3:5],1,max)

# Get the mean (max Q)
pop.admixed=data.frame(pop=unique(ancestry$pop))
i=0
for (p in unique(ancestry$pop)) {
  print(p)
  i=i+1
  pop.admixed$Q[i]=mean(ancestry$Qmax[which(ancestry$pop==p)])
  pop.admixed$pop[i]=as.character(labelPops$Location_name[which(labelPops$Location_index==p)])
}
pop.admixed$Q=round(pop.admixed$Q,digits = 3)
write.table(x=pop.admixed,file="Tables/Major ancestry coeff Invasive K3b.txt",sep="\t",quote = FALSE,row.names = FALSE)






# #----------------------------------------------------------
# #         INVASIVE POPULATIONS IN EUROPE
# #----------------------------------------------------------
# 
# # Population clustering on worldwide popualtions, without any a priori on populations
# evannoWorld=read.table(paste(STRdir,"/1a_Allpops/Harvester/evanno.txt",sep=""),header=F,sep="\t")
# # Population clustering on worldwide popualtions, with populations in Asia fixed
# evannoWorld=read.table(paste(STRdir,"/1b_Allpops/Harvester/evanno.txt",sep=""),header=F,sep="\t")
# 
# colnames(evannoWorld)=c("K","Reps","Mean_LnP","sd_LnP","LnK","Ln2K","DeltaK")
# 
# ######## Mean Ln(P) as a function of K
# (pLnP=ggplot(data=evannoWorld, aes(x=K, y=Mean_LnP)) +
#     geom_point(colour="Black",size=1)+
#     geom_line() +
#     geom_errorbar(aes(ymin=Mean_LnP-sd_LnP, ymax=Mean_LnP+sd_LnP), width=.2) +
#     xlab("K") + ylab("Ln(P)") +
#     scale_x_discrete(limits=c(2:nrow(evannoWorld))) +
#     theme(axis.line = element_line(colour = "black"),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(),
#           panel.background = element_blank(),
#           plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
#           plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
#           axis.title.x = element_text(color="black", size=14),
#           axis.title.y = element_text(color="black", size=14),
#           axis.text=element_text(size=14, colour="black"),
#           legend.key = element_rect(fill = "white", size = 1),
#           legend.text=element_text(size=14),
#           legend.title=element_text(size=14)))
# 
# 
# ######## Mean Ln(K) as a function of K
# (pLnK=ggplot(data=evannoWorld, aes(x=K, y=LnK)) +
#     geom_point(colour="Black",size=1)+
#     geom_line() +
#     xlab("K") + ylab("Ln'(K)") +
#     scale_x_discrete(limits=c(2:nrow(evannoWorld))) +
#     theme(axis.line = element_line(colour = "black"),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(),
#           panel.background = element_blank(),
#           plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
#           plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
#           axis.title.x = element_text(color="black", size=14),
#           axis.title.y = element_text(color="black", size=14),
#           axis.text=element_text(size=14, colour="black"),
#           legend.key = element_rect(fill = "white", size = 1),
#           legend.text=element_text(size=14),
#           legend.title=element_text(size=14)))
# 
# 
# ######## Mean Ln''(K) as a function of K
# (pLn2K=ggplot(data=evannoWorld, aes(x=K, y=Ln2K)) +
#     geom_point(colour="Black",size=1)+
#     geom_line() +
#     xlab("K") + ylab("Ln''(K)") +
#     scale_x_discrete(limits=c(2:nrow(evannoWorld))) +
#     theme(axis.line = element_line(colour = "black"),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(),
#           panel.background = element_blank(),
#           plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
#           plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
#           axis.title.x = element_text(color="black", size=14),
#           axis.title.y = element_text(color="black", size=14),
#           axis.text=element_text(size=14, colour="black"),
#           legend.key = element_rect(fill = "white", size = 1),
#           
#           legend.text=element_text(size=14),
#           legend.title=element_text(size=14)))
# 
# 
# ######## Mean Ln''(K) as a function of K
# (pDeltaK=ggplot(data=evannoWorld, aes(x=K, y=DeltaK)) +
#     geom_point(colour="Black",size=1)+
#     geom_line() +
#     xlab("K") + ylab("Delta(K)") +
#     scale_x_discrete(limits=c(2:nrow(evannoWorld))) +
#     theme(axis.line = element_line(colour = "black"),
#           panel.grid.major = element_blank(),
#           panel.grid.minor = element_blank(),
#           panel.border = element_blank(),
#           panel.background = element_blank(),
#           plot.title = element_text(color="black", size=14, face="bold.italic",hjust = 0.5),
#           plot.subtitle = element_text(color="black",size=14,hjust = 0.5),
#           axis.title.x = element_text(color="black", size=14),
#           axis.title.y = element_text(color="black", size=14),
#           axis.text=element_text(size=14, colour="black"),
#           legend.key = element_rect(fill = "white", size = 1),
#           
#           legend.text=element_text(size=14),
#           legend.title=element_text(size=14)))
# 
# #### SAVE TO FIGURE
# ggarrange(pLnP,pLnK,pLn2K,pDeltaK,widths=1:1,heights=1:1,labels="auto")
# ggsave(paste(figuresdir,"/STRUCTURE/1b_Allpops/Diagnostic_plots_STRUCTURE_native.png",sep=""),
#        device="png",dpi=320,units="cm",width=40,height=32)

################################
# For 1a_Allpops: burnin 100,000 + 100,000, K=1-32, 20 replicates, 2,112 loci


# For 1b_Allpops: burnin 100,000 + 100,000, K=1-21, 20 replicates, 2,112 loci, PFROMPOPFLAGONLY
# As for K=1-122, using Asian population as fixed for population assignment didn't seemed to work at all. K=2 and linear increase of likelihood with a lot of variance up to K=11
# Wait and see results when K>12, but no reason to be optimist


# For 1c_Allpops: burnin 100,000 + 100,000, K=1-21, 20 replicates, 2,112 loci, PFROMPOPFLAGONLY, use demes instead of pops in Asia



#==========================================================
# Population Assignment with AssignPOP - Maps
#==========================================================

# See 'AssignPOP.R' for script... Here are only presented the codelines for a plot
# of the invasive area with pie charts of predicted populations for sampled sites

#----------------------------------------------------------
# Fig 3. Map of the invasive area with pies chart of admixture proportions on each sampling location

# Prepare the data set for mapping
res=read.table(("AssignPOP/assignPOP.MC.svm.demes/TrimAssignmentResult.txt"), header=TRUE)
res.table = data.frame(Population = row.names(table(res$Ind.ID, res$pred.pop)),
                       Coastal_Admixed_China = table(res$Ind.ID, res$pred.pop)[,1],
                       Continental_Admixed_China = table(res$Ind.ID, res$pred.pop)[,2],
                       Japan = table(res$Ind.ID, res$pred.pop)[,3],
                       North_Central_China = table(res$Ind.ID, res$pred.pop)[,4],
                       North_China = table(res$Ind.ID, res$pred.pop)[,5],
                       South_China = table(res$Ind.ID, res$pred.pop)[,6]) # predicted populations as a function of invasive population

Invasive.pies.means=sapply(split(res.table[2:7],res.table$Population), colMeans) ## This just makes population cluster averages for each cluster

# Mapping
# Set a vector of colors for clusters
cols.assign=c("#E41A1C", "#4DAF4A", "#377EB8", "#FF7F00", "#984EA3",  "#FFFF33")
# Download database of river network at scale 50
rivers50=ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')
# cols.invasive=brewer.pal(n=8,name="Pastel1")
png("Figures/AssignPOP/Map predicted populations.png",width=2000,height=1100)
map("world", xlim=c(-15,60), ylim=c(30,59), col="gray90",  fill=TRUE) # Plot the map of European area
# map(database="worldHires", regions="Caspian Sea", col="white", fill=TRUE, add=TRUE) # Correct for unfilled Caspian sea
# map(database="worldHires", regions="Aral Sea", col="white", fill=TRUE, add=TRUE)# Correct for unfilled Aral sea
sp::plot(rivers50,xlim=c(-15,30), ylim=c(30,59),lwd=2, col = 'blue',add=TRUE)
# Nomenclature
north.arrow(56,56,1, lab="N", lab.pos="above")
scalebar(c(49,31),1000, bg=NA, border=NA, division.cex=2)
text(x=105, y=34,"",cex=3,font=2)
text(x=139, y=39,"",cex=3,font=2)
points(invasive.coords$X, invasive.coords$Y, pch = 16, cex = 0.7,col="red") #\ Add points to the map
## Pies ## 
# Draw a pie with admixture proportions for each population
for (i in 1:ncol(Invasive.pies.means)) {
  add.pie(Invasive.pies.means[,i],x=invasive.coords$X[i], y=invasive.coords$Y[i],labels=invasive.coords[i,1],label.dist=1.5,
          radius=1,edges=200,clockwise=T,	col	=	cols.assign,cex=2,font=2)
  # Add labels to pies
}
dev.off()



#==========================================================
# FINAL MAP FOR REPORT
#==========================================================
# On native side, structure pie charts
# On invasive side, population assignment with AssignPOP

# Prepare the data set for mapping
strK6=read.table(paste(STRdir,"/2b_Native/CLUMPP/K6.outfile",sep=""))
Native.pies.str=cbind(as.character(labelPops$Location_name[strK6[,4]]),as.data.frame(strK6[,(c(1,2,3,5,4,6)+5)])) # Clusters are selected to correspodn to DAPC colors
names(Native.pies.str)=c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")
Native.pie.str_means=sapply(split(Native.pies.str[2:7],Native.pies.str$Population), colMeans) ## This just makes population cluster averages for each cluster

# Prepare the data set for mapping
res=read.table(("AssignPOP/assignPOP.MC.svm.demes/TrimAssignmentResult.txt"), header=TRUE)
res.table = data.frame(Population = row.names(table(res$Ind.ID, res$pred.pop)),
                       Coastal_Admixed_China = table(res$Ind.ID, res$pred.pop)[,1],
                       Continental_Admixed_China = table(res$Ind.ID, res$pred.pop)[,2],
                       Japan = table(res$Ind.ID, res$pred.pop)[,3],
                       North_Central_China = table(res$Ind.ID, res$pred.pop)[,4],
                       North_China = table(res$Ind.ID, res$pred.pop)[,5],
                       South_China = table(res$Ind.ID, res$pred.pop)[,6]) # predicted populations as a function of invasive population
Invasive.pies.means=sapply(split(res.table[2:7],res.table$Population), colMeans) ## This just makes population cluster averages for each cluster

# Download database of river network at scale 50
# rivers50=ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')
# save(rivers50,file="Data/rivers50.Rda")
load("Data/rivers50.Rda")
rivers50.native = crop(rivers50, extent(-15, 30, 30, 59))
rivers50.invasive = crop(rivers50, extent(92, 160, 20, 48))

png("Figures/Worldwide map Report.png",width=3000,height=1000)
map("world", xlim=c(-15,150), ylim=c(20,59), col="gray90",  fill=TRUE) # Plot the map of European area
# Nomenclature
north.arrow(149,56,2, lab="", lab.pos="above")
scalebar(c(128,20),2000, bg=NA, border=NA, division.cex=3)
sp::plot(rivers50.native, lwd=2, col = 'blue',add=TRUE) # invasive range river network
sp::plot(rivers50.invasive, lwd=2, col = 'blue',add=TRUE) # native range river network
# INVASIVE
points(invasive.coords$X, invasive.coords$Y, pch = 16, cex = 4,col="grey49") # Add points to the map
## Pies ## 
# Draw a pie with assignment proportions for each population
# Set a vector of colors for clusters
cols.assign=c("#E41A1C", "#4DAF4A", "#377EB8", "#FF7F00", "#984EA3",  "#FFFF33")
for (i in 1:ncol(Invasive.pies.means)) {
  add.pie(Invasive.pies.means[,i],x=invasive.coords$X[i], y=invasive.coords$Y[i],labels="",
          radius=1.3,edges=200,clockwise=T,	col	=	cols.assign,cex=2,font=2)
  # Add labels to pies
}
# NATIVE
# Labels of sampled sites
text(invasive.coords$X, y = invasive.coords$Y, invasive.coords$Pop,
     pos = c(rep(4,10), 2), offset =2.5, cex = 3, font = 2)
points(native.coords$X, native.coords$Y, pch = 16, cex = 4,col="black") # Add points to the map
points(native.coords$X[c(11,13,18)], native.coords$Y[c(11,13,18)], pch = 16, cex = 4,col="grey49") #\ Add points to the map
# Labels of sampled sites
# Jitter the 13rd label
text(native.coords$X, y = native.coords$Y, native.coords$Pop,
     pos = c(4, 4, 2, 2, 4, 4, 4, 2, 2, 2, 4, 4, 4, 4, 4, 4, 2, 4), offset = 2.5, cex = 3, font = 2)
## Pies ## 
# Draw a pie with admixture proportions for each population
# Set a vector of colors for clusters
cols.native=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
              "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9")
for (i in 1:ncol(Native.pie.str_means)) {
  add.pie(Native.pie.str_means[,i],x=native.coords$X[i], y=native.coords$Y[i], labels = "",
          radius=1.3,edges=200,clockwise=T,	col	=	cols.native,cex=2,font=2)
  # Add labels to pies
}
dev.off()




#==========================================================
# FINAL MAP FOR ARTICLE
#==========================================================
# On native side, structure pie charts
# On invasive side, population assignment with AssignPOP
# Underneath the map a Barplot STRUCTURE + ASSIGNPOP

# Prepare the data set for mapping
strK6=read.table(paste(STRdir,"/2b_Native/CLUMPP/K6.outfile",sep=""))
Native.pies.str=cbind(as.character(labelPops$Location_name[strK6[,4]]),as.data.frame(strK6[,(c(1,2,3,5,4,6)+5)])) # Clusters are selected to correspodn to DAPC colors
names(Native.pies.str)=c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4", "Cluster 5", "Cluster 6")
Native.pie.str_means=sapply(split(Native.pies.str[2:7],Native.pies.str$Population), colMeans) ## This just makes population cluster averages for each cluster

# Prepare the data set for mapping
res=read.table(("AssignPOP/assignPOP.MC.svm.demes/TrimAssignmentResult.txt"), header=TRUE)
res.table = data.frame(Population = row.names(table(res$Ind.ID, res$pred.pop)),
                       Coastal_Admixed_China = table(res$Ind.ID, res$pred.pop)[,1],
                       Continental_Admixed_China = table(res$Ind.ID, res$pred.pop)[,2],
                       Japan = table(res$Ind.ID, res$pred.pop)[,3],
                       North_Central_China = table(res$Ind.ID, res$pred.pop)[,4],
                       North_China = table(res$Ind.ID, res$pred.pop)[,5],
                       South_China = table(res$Ind.ID, res$pred.pop)[,6]) # predicted populations as a function of invasive population
Invasive.pies.means=sapply(split(res.table[2:7],res.table$Population), colMeans) ## This just makes population cluster averages for each cluster

# Download database of river network at scale 50
# rivers50=ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')
# save(rivers50,file="Data/rivers50.Rda")
load("Data/rivers50.Rda")
rivers50.native = crop(rivers50, extent(-15, 30, 30, 59))
rivers50.invasive = crop(rivers50, extent(92, 160, 20, 48))

svg("Figures/Worldwide map Report.svg",width=60,height=20)
map("world", xlim=c(-15,150), ylim=c(20,59), col="gray90",  fill=TRUE) # Plot the map of European area
# Nomenclature
north.arrow(149,56,2, lab="", lab.pos="above")
scalebar(c(128,20),2000, bg=NA, border=NA, division.cex=3)
sp::plot(rivers50.native, lwd=2, col = 'blue',add=TRUE) # invasive range river network
sp::plot(rivers50.invasive, lwd=2, col = 'blue',add=TRUE) # native range river network
# INVASIVE
points(invasive.coords$X, invasive.coords$Y, pch = 16, cex = 4,col="grey49") # Add points to the map
## Pies ## 
# Draw a pie with assignment proportions for each population
# Set a vector of colors for clusters
cols.assign=c("#E41A1C", "#4DAF4A", "#377EB8", "#FF7F00", "#984EA3",  "#FFFF33")
for (i in 1:ncol(Invasive.pies.means)) {
  add.pie(Invasive.pies.means[,i],x=invasive.coords$X[i], y=invasive.coords$Y[i],labels="",
          radius=1.3,edges=200,clockwise=T,	col	=	cols.assign,cex=2,font=2)
  # Add labels to pies
}
# NATIVE
# Jitter S19 & S20

# Labels of sampled sites
text(invasive.coords$X, y = invasive.coords$Y, invasive.coords$Pop,
     pos = c(rep(4,10), 2), offset = 3.5, cex = 3, font = 2)
points(native.coords$X, native.coords$Y, pch = 16, cex = 4,col="black") # Add points to the map
points(native.coords$X[c(11,13,18)], native.coords$Y[c(11,13,18)], pch = 16, cex = 4,col="grey49") #\ Add points to the map
# Labels of sampled sites
# Jitter the 13rd label
text(native.coords$X, y = native.coords$Y, native.coords$Pop,
     pos = c(4, 4, 2, 2, 4, 4, 4, 2, 2, 2, 4, 4, 4, 4, 4, 4, 2, 4), offset = 3.5, cex = 3, font = 2)
## Pies ## 
# Draw a pie with admixture proportions for each population
# Set a vector of colors for clusters
cols.native=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
              "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9")
for (i in 1:ncol(Native.pie.str_means)) {
  add.pie(Native.pie.str_means[,i],x=native.coords$X[i], y=native.coords$Y[i], labels = "",
          radius=1.3,edges=200,clockwise=T,	col	=	cols.native,cex=2,font=2)
  # Add labels to pies
}
dev.off()


# A barplot with AssignPop Invasive on the left and STRUCTURE Native on the right
#-------------------------------------
### Bar plot representation ###
# Representation of admixture (Q) for each individual FOR K=6 (chosen K)
# Format a data frame with 4 columns, and as many rows per individual as there is clusters
# strK6=read.table(paste(STRdir,"/2a_Native/CLUMPP/K6.outfile",sep=""))
strK6=read.table(paste(STRdir,"/2b_Native/CLUMPP/K6.outfile",sep=""))
# Set a vector of colors for clusters
cols.native=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
              "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9")

# Get the same colors and cluster number as DAPC
strK6.results=as.data.frame(strK6[,(c(1,2,3,5,4,6)+5)]) # population clusters in the good order for colors
colnames(strK6.results)=c(1,2,3,4,5,6)
# strK6.results$ind=as.character(strK6[,1])
strK6.results$pop=as.character(strK6[,4])
for (i in 1:nrow(strK6.results)) {
  strK6.results$pop[i]=as.character(labelPops[which(strK6.results$pop[i] == labelPops$Location_index),1])
}
strK6.results$indNames=names(Native$tab[,1])
strK6.results=melt(strK6.results)

colnames(strK6.results)=c("Population", "Individual","Cluster","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK6.results$Population=factor(strK6.results$Population,
                                levels = c("S9","S18","S19","S20","Tib","S10","S11",
                                           "S4","S6",
                                           "S3","S1","S2","S16","S13","S14",
                                           "S15","S17", "Jap"))

res.trim = read.table(file = "AssignPOP/assignPOP.MC.svm.demes/TrimAssignmentResult.txt", header = TRUE)

# cols.demes = c("#4DAF4A", "#377EB8", "#FF7F00") # Blue, Green, Orange
# cols.invasive=brewer.pal(n = 3, name = "Set2")
# cols.invasive = c("#8DA0CB", "#FC8D62", "#66C2A5")

# Replace pop name by deme name, for legend
res.trim$Invasive.deme = as.character(res.trim$Ind.ID)
res.trim$Invasive.deme[res.trim$Invasive.deme %in% c("Ira")] = "Iran"
res.trim$Invasive.deme[res.trim$Invasive.deme %in% c("Bul1", "Bul2", "Tur")] = "Eastern Europe"
res.trim$Invasive.deme[res.trim$Invasive.deme %in% c("Aus", "Bel", "Hun", "Ita", "Pol", "Spa", "UK")] = "Western Europe"

# Demes are displayed in this order...
res.trim$pred.pop=factor(res.trim$pred.pop, levels = c("Japan","South_China","Coastal_Admixed_China",
                                                       "Continental_Admixed_China","North_Central_China",
                                                       "North_China"))
# res.trim$Ind.ID=factor(res.trim$Ind.ID, levels = c("UK", "Bel", "Spa", "Ita", "Aus", "Pol", "Hun", "Bul1", "Bul2", "Tur", "Ira"))
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
df.melt$Pop = factor(df.melt$Pop, levels = c("UK", "Bel", "Spa", "Ita", "Aus", "Pol", "Hun", "Bul", "Tur", "Ira"),
                     labels = c("UK", "Bel", "Spa", "Ita", "Aus", "Pol", "Hun", "Bul", "Tur", "Ira"))



# Plotting the barplot
png("Figures/Worldwide map Report STRUCTURE barplot.png",width=3000,height=500)
txtsize = 30
#---------------------------------
# Beautiful plot of assignment proba
PredSource.prop.plot = ggplot(data=df.melt, aes(x=Ind, y = Posterior_assignment_probability, fill=source.demes))+
  geom_col(width=1) +
  # facet_grid(~inv.demes, scales = "free",space="free_x") +
  facet_grid(~Pop, scales = "free",space="free_x") +
  scale_fill_manual(values = cols.demes, labels = c("Japan" = "Japan",
                                                    "South_China" = "South China",
                                                    "Coastal_Admixed_China" = "Coastal China",
                                                    "Continental_Admixed_China" = "Continental China",
                                                    "North_Central_China" = "North Central China",
                                                    "North_China" = "North China")) +
  labs(x="Sampled individual", y="Posterior assignment probability", fill="Native deme") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=txtsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=txtsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=txtsize),
        axis.title.y = element_text(color="black", size=txtsize),
        axis.text=element_text(size=txtsize, colour="black"),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=26, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=txtsize),
        legend.title=element_text(size=txtsize))

strK6_barplot = ggplot(data=strK6.results, aes(x=Individual, y=Posterior_membership_probability, fill=Cluster))+
  geom_col(width=1) +
  # geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.native) +
  facet_grid(~Population, scales = "free",space="free_x") +
  labs(x="Sampled individual", y="Posterior admixture proportions\n", fill="Genetic cluster") +
  theme(axis.line = element_blank(),
        # axis.line.x = element_blank(), # No x axis
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        plot.title = element_text(color="black", size=txtsize, face="bold.italic",hjust = 0.5),
        plot.subtitle = element_text(color="black",size=txtsize,hjust = 0.5),
        axis.title.x = element_text(color="black", size=txtsize),
        axis.title.y = element_text(color="black", size=txtsize),
        axis.text=element_text(size=txtsize, colour="black"),
        axis.text.x=element_blank(), # No samples names
        axis.ticks.x=element_blank(), # No x axis
        strip.text=element_text(size=26, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=txtsize),
        legend.title=element_text(size=txtsize))

library(egg)
ggarrange(PredSource.prop.plot, strK6_barplot, nrow = 1, ncol = 2, widths = c(2,3))

# ggsave(paste(figuresdir,"/STRUCTURE/2b_Native/Native_barplot_K6b.png",sep=""),
#        device="png",dpi=320,units="cm",width=60,height=9)
dev.off()











#==========================================================
# FINAL MAP FOR ORAL DEFENSE
#==========================================================
# On both sides, structure pie charts

# strK9Eu=read.table(paste(STRdir,"/3a_Invasive/CLUMPP/K9.outfile",sep=""))
# Invasive.pies.str=cbind(as.character(labelPops$Location_name[strK9Eu[,4]]),as.data.frame(strK9Eu[,(c(1,2,3,4,5,6,7,8,9)+5)])) # Clusters are selected to correspodn to DAPC colors
# names(Invasive.pies.str)=c("Population", "Cluster 1", "Cluster 2", "Cluster 3", "Cluster 4",
#                            "Cluster 5", "Cluster 6", "Cluster 7", "Cluster 8", "Cluster 9")
# 
# Invasive.pies.str_means=sapply(split(Invasive.pies.str[2:10],Invasive.pies.str$Population), colMeans) ## This just makes population cluster averages for each cluster

strK3Eu=read.table(paste(STRdir,"/3a_Invasive/CLUMPP/K3.outfile",sep=""))
Invasive.pies.str=cbind(as.character(labelPops$Location_name[strK3Eu[,4]]),as.data.frame(strK3Eu[,(c(1,2,3)+5)])) # Clusters are selected to correspodn to DAPC colors
names(Invasive.pies.str)=c("Population", "Cluster 1", "Cluster 2", "Cluster 3")

Invasive.pies.str_means=sapply(split(Invasive.pies.str[2:4],Invasive.pies.str$Population), colMeans) ## This just makes population cluster averages for each cluster

#----------------------------------
# WORLDWIDE MAP OF SAMPLED SITES FOR M&M (ALL IN ONE MAP)
# Download database of river network at scale 50
# rivers50=ne_download(scale = 50, type = 'rivers_lake_centerlines', category = 'physical')
# save(rivers50,file="Data/rivers50.Rda")
load("Data/rivers50.Rda")
rivers50.native = crop(rivers50, extent(-15, 30, 30, 59))
rivers50.invasive = crop(rivers50, extent(92, 160, 20, 48))

png("Figures/Worldwide map Oral.png",width=3000,height=1000)
map("world", xlim=c(-15,150), ylim=c(20,59), col="gray90",  fill=TRUE) # Plot the map of European area
# Nomenclature
north.arrow(149,56,2, lab="", lab.pos="above")
scalebar(c(125,20),2000, bg=NA, border=NA, division.cex=3)
# sp::plot(rivers50.native, lwd=2, col = 'blue',add=TRUE) # invasive range river network
# sp::plot(rivers50.invasive, lwd=2, col = 'blue',add=TRUE) # native range river network
# INVASIVE
points(invasive.coords$X, invasive.coords$Y, pch = 16, cex = 4,col="grey49") # Add points to the map
## Pies ## 
# Draw a pie with admixture proportions for each population
# Set a vector of colors for clusters
cols.invasive = c("#ABD9E9", "#FEC44F", "#1B5E20", "#A6DBA0")
for (i in 1:ncol(Invasive.pies.str_means)) {
  add.pie(Invasive.pies.str_means[,i],x=invasive.coords$X[i], y=invasive.coords$Y[i],labels="",
          radius=1.3,edges=200,clockwise=T,	col	=	cols.invasive,cex=2,font=2)
  # Add labels to pies
}
# NATIVE
# Labels of sampled sites
text(invasive.coords$X, y = invasive.coords$Y, invasive.coords$Pop,
     pos = c(rep(4,10), 2), offset =2.5, cex = 3, font = 2)
points(native.coords$X, native.coords$Y, pch = 16, cex = 4,col="black") # Add points to the map
points(native.coords$X[c(11,13,18)], native.coords$Y[c(11,13,18)], pch = 16, cex = 4,col="grey49") #\ Add points to the map
# Labels of sampled sites
# Jitter the 13rd label
text(native.coords$X, y = native.coords$Y, native.coords$Pop,
     pos = c(4, 4, 2, 2, 4, 4, 4, 2, 2, 2, 4, 4, 4, 4, 4, 4, 2, 4), offset = 2.5, cex = 3, font = 2)
### Pies ## 
# Draw a pie with admixture proportions for each population
# Set a vector of colors for clusters
cols.native=c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
              "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9")
for (i in 1:ncol(Native.pie.str_means)) {
  add.pie(Native.pie.str_means[,i],x=native.coords$X[i], y=native.coords$Y[i], labels = "",
          radius=1.3,edges=200,clockwise=T,	col	=	cols.native,cex=2,font=2)
  # Add labels to pies
}
dev.off()













#==========================================================
# Report cover
#----------------------------------------------------------

#----------------------------------
# WORLDWIDE MAP
# European countries marked as invaded
europe <- c("Austria","Belgium","Bulgaria","Croatia", "Serbia","Bosnia",
            "Montenegro","Albania","Moldova", "Kosovo", "Macedonia",
                   "Czech Rep.","Denmark","Estonia","France",
                   "Germany","Greece","Hungary","Italy","Latvia",
                   "Lithuania","Luxembourg","Netherlands","Poland",
                   "Romania","Slovakia","Slovenia","Spain",
                   "United Kingdom", "UK", "Ukrain")

png("Figures/Report cover.png",width=2000,height=1000)
par(mai=c(0,0,0,0), mar=c(0,0,0,0))
map("world", col="gray45", border = "gray15", bg = "gray15", fill=TRUE, mar = c(0,0,0,0)) # Plot the worldwide dark grey map with darker seas
map("world", regions = c("China", "North Korea", "South Korea", "Japan"),
    col="darkred", border="darkred",  fill=TRUE, add = TRUE) # Plot the native areas
map("world", regions = c(europe, "Europe", "Morocco", "Algeria",
                         "Turkey", "Iran"), col="red2", border="red2",  fill=TRUE, add = TRUE) # Plot the invasive areas
dev.off()


map("world", col="gray45", bg = "gray15", fill=TRUE, mar = c(0,0,0,0), namesonly = TRUE) # Plot the worldwide dark grey map with darker seas



#==========================================================
# Presentation invasion history
#----------------------------------------------------------
# 2 maps of countries known as invaded (1960-69 and 2000-2009)

#----------------------------------
# European countries marked as invaded in 1960
europe <- c("Hungary","Lithuania","Romania")

png("Figures/Presentation invasion 1960.png",width=2000,height=1000)
par(mai=c(0,0,0,0), mar=c(0,0,0,0))
map("world", xlim=c(-15,65), ylim=c(25,59), col="gray45", border = "gray15", bg = "gray15", fill=TRUE, mar = c(0,0,0,0)) # Plot the worldwide dark grey map with darker seas
map("world", regions = c(europe), col="red2", border="red2",  fill=TRUE, add = TRUE) # Plot the invasive areas
dev.off()


#----------------------------------
# European countries marked as invaded in 2000
europe <- c("Austria","Belgium","Bulgaria","Croatia", "Serbia","Bosnia",
            "Montenegro","Albania","Moldova", "Kosovo", "Macedonia",
            "Czech Rep.","Denmark","Estonia","France",
            "Germany","Greece","Hungary","Italy","Latvia",
            "Lithuania","Luxembourg","Netherlands","Poland",
            "Romania","Slovakia","Slovenia","Spain",
            "United Kingdom", "UK", "Ukrain")

png("Figures/Presentation invasion 2000.png",width=2000,height=1000)
par(mai=c(0,0,0,0), mar=c(0,0,0,0))
map("world", xlim=c(-15,65), ylim=c(25,59), col="gray45", border = "gray15", bg = "gray15", fill=TRUE, mar = c(0,0,0,0)) # Plot the worldwide dark grey map with darker seas
map("world", regions = c(europe, "Morocco", "Algeria",
                         "Turkey", "Iran"), col="red2", border="red2",  fill=TRUE, add = TRUE) # Plot the invasive areas
dev.off()







#==========================================================
# THE END
#==========================================================