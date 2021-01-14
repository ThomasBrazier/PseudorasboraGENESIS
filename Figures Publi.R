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
# FINAL STRUCTURE BARPLOT FOR ARTICLE
#==========================================================
#-------------------------------------
### Bar plot representation ###
# Representation of admixture (Q) for each individual FOR K=6 (chosen K)
# Format a data frame with 4 columns, and as many rows per individual as there is clusters
# strK6=read.table(paste(STRdir,"/2a_Native/CLUMPP/K6.outfile",sep=""))
strK6 = read.table(paste(STRdir,"/2b_Native/CLUMPP/K6.outfile",sep=""))
# Set a vector of colors for clusters
cols.native = c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3", "#FF7F00", "#FFFF33", "#A65628", "#F781BF",
              "#8DD3C7", "#FFFFB3", "#BEBADA", "#FB8072", "#80B1D3", "#FDB462", "#B3DE69", "#FCCDE5", "#D9D9D9")

# Get the same colors and cluster number as DAPC
strK6.results=as.data.frame(strK6[,(c(1,2,3,5,4,6)+5)]) # population clusters in the good order for colors
colnames(strK6.results) = c(1,2,3,4,5,6)
# strK6.results$ind=as.character(strK6[,1])
strK6.results$pop = as.character(strK6[,4])
strK6.results$deme = as.character(strK6[,4])
labelPops$deme = c(rep("Europe", 7), "Japan", "Europe", "North Central China", "Continental", "Continental", "North China", "North China", "North China",
                   "North Central China", "North China", "South China", "South China", "North Central China", "South China", " ", "Coastal", "Coastal",
                   "South China", "Europe", "Tibet", "Europe", "Europe")
for (i in 1:nrow(strK6.results)) {
  strK6.results$pop[i] = as.character(labelPops[which(strK6.results$pop[i] == labelPops$Location_index),1])
  strK6.results$deme[i] = as.character(labelPops$deme[which(strK6.results$deme[i] == labelPops$Location_index)])
}
strK6.results$indNames=names(Native$tab[,1])
strK6.results=melt(strK6.results)

colnames(strK6.results)=c("Population", "Deme", "Individual","Cluster","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK6.results$Population=factor(strK6.results$Population,
                                levels = c("S9","S18","S19","S20","Tib","S10","S11",
                                           "S4","S6",
                                           "S3","S1","S2","S16","S13","S14",
                                           "S15","S17", "Jap"))
strK6.results$Deme=factor(strK6.results$Deme,
                                levels = c("South China", "Tibet", "Continental", " ", "Coastal", "North Central China", "North China", "Japan"))
# EXPORT SVG FIGURE 1A
# Plotting the STRUCTURE barplot
# source("sources/facet_nested.R")
# devtools::install_github("teunbrand/ggh4x")
library(ggh4x)
txtsize = 30
strK6_barplot = ggplot(data = strK6.results, aes(x=Individual, y=Posterior_membership_probability, fill=Cluster))+
  geom_col(width=1) +
  # geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.native) +
  # facet_grid(~ Population, scales = "free",space="free_x") +
  facet_nested(~ Deme + Population, scales = "free",space="free_x") +
  labs(x="Sampled individual", y="Admixture proportions\n", fill="Genetic cluster") +
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
        strip.background = element_rect(colour="white", fill="white"),
        legend.position = "none")
# strK6_barplot
# ggsave("Figures/ARTICLE STRUCTURE barplot.svg", plot = strK6_barplot, device = "svg", units = "cm", width = 60, height = 15)

# Background color of populations facets grid
japan_col = "#377EB8"
coastal_col = "#E41A1C"
continental_col = "#4DAF4A"
north_central_col = "#FF7F00"
north_col = "#984EA3"
south_col = "#FFFF33"
na_col = "darkgrey"
color_demes = c(south_col, south_col, south_col, south_col,
                na_col, continental_col, continental_col,
                na_col, coastal_col, coastal_col,
                north_central_col, north_central_col, north_central_col,
                north_col, north_col, north_col, north_col,
                japan_col)
g = ggplot_gtable(ggplot_build(strK6_barplot))
strip = which(grepl('strip-t-[0-9]*-2', g$layout$name))
for (i in 1:length(strip)) {
  j = which(grepl('rect', g$grobs[[strip[i]]]$grobs[[1]]$childrenOrder))
  g$grobs[[strip[i]]]$grobs[[1]]$children[[j]]$gp$fill <- color_demes[i]
}

ggsave("Figures/ARTICLE STRUCTURE barplot.svg", plot = g, device = "svg", units = "cm", width = 60, height = 15)



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

# Jitter S19 & S20 pie charts
native.coords$Y[native.coords$Pop == "S19"] = native.coords$Y[native.coords$Pop == "S19"] - 1
native.coords$Y[native.coords$Pop == "S20"] = native.coords$Y[native.coords$Pop == "S20"] + 1

# Jitter Bul1 & Bul2 pie charts
invasive.coords$Y[invasive.coords$Pop == "Bul2"] = invasive.coords$Y[invasive.coords$Pop == "Bul2"] - 1
invasive.coords$Y[invasive.coords$Pop == "Bul1"] = invasive.coords$Y[invasive.coords$Pop == "Bul1"] + 1



res.trim = read.table(file = "AssignPOP/assignPOP.MC.svm.demes/TrimAssignmentResult.txt", header = TRUE)
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
df.melt$Pop = as.character(df.melt$Pop)
df.melt$Pop[df.melt$Pop %in% c("Bul1", "Bul2")] = "Bul"
df.melt$Pop = as.factor(df.melt$Pop)
df.melt$Pop = factor(df.melt$Pop, levels = c("UK", "Bel", "Spa", "Ita", "Aus", "Pol", "Hun", "Bul", "Tur", "Ira"),
                     labels = c("UK", "Bel", "Spa", "Ita", "Aus", "Pol", "Hun", "Bul", "Tur", "Ira"))


# EXPORT SVG FIGURE 2A
svg("Figures/ARTICLE Worldwide map Piecharts.svg",width=60,height=20)
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
     pos = c(rep(4,10), 2), offset = 3.5, cex = 4, font = 2)
points(native.coords$X, native.coords$Y, pch = 16, cex = 4,col="black") # Add points to the map
points(native.coords$X[c(11,13,18)], native.coords$Y[c(11,13,18)], pch = 16, cex = 4,col="grey49") #\ Add points to the map
# Labels of sampled sites
# Jitter the 13rd label
text(native.coords$X, y = native.coords$Y, native.coords$Pop,
     pos = c(4, 4, 2, 2, 4, 4, 4, 2, 2, 2, 4, 4, 4, 4, 4, 4, 2, 4), offset = 3.5, cex = 4, font = 2)
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

# EXPORT SVG FIGURE 2B
# Plotting the AssignPop barplot
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
ggsave("Figures/ARTICLE Worldwide map AssignPop.svg", plot = PredSource.prop.plot, device = "svg", units = "cm", width = 60, height = 20)

#png(output_file, width=800, height=400)
# Cairo(800,400,file=paste(output_file, ".svg", sep=""),type="svg",bg="transparent",pointsize=8, units="px",dpi=400)
# svg("Figures/ARTICLE Worldwide map AssignPop.svg", width = 60, height = 10)
# gt <- ggplot_gtable(ggplot_build(p))
# gt$layout$clip[gt$layout$name=="panel"] <- "off"
# grid.draw(PredSource.prop.plot)
# dev.off()


#==========================================================
# THE END
#==========================================================