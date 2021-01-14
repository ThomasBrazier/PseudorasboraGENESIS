# 1/ Graphic evaluation with the four plots of Evanno's method (Evanno et al. 2005)
# Most important plot is DeltaK = mean(|L''(K)|)/sd(L(K)) of L(K) as a function of K 

### Choose the most probable K value (number of cluster inferred), K value the most appropriate for the data
# Load results from Harvester: 'evanno.txt'
setwd("~/Documents/2b_Native_tmp/Harvester")
evanno=read.table("evanno.txt")
setwd("~/Documents/2b_Native_tmp/")

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
# ggsave(paste(figuresdir,"/STRUCTURE/b_Diagnostic_plots_STRUCTURE_native.tiff",sep=""),
#        device="tiff",dpi=320,units="cm",width=40,height=32)

# 2/ Quantitative evaluation with the G' (or H') parameter computed by CLUMPP to assess convergence of STRUCTURE



### Assess the convergence for parameter estimates ###
# See STRUCTURE.R for convergence assessment...

### Convergence of Q


### Convergence of Alpha (admixture parameter) per individual

#-------------------------------------
### Bar plot representation ###
# Representation of admixture (Q) for each individual FOR K=6 (chosen K)
# Format a data frame with 4 columns, and as many rows per individual as there is clusters
setwd("~/Documents/2b_Native_tmp/Harvester")
strK6=read.table("K6.indfile")
setwd("~/Documents/2b_Native_tmp/")
# Get the same colors and cluster number as DAPC
strK6.results=as.data.frame(strK6[1:300,(c(1,2,3,4,5,6)+5)]) # population clusters in the good order for colors

colnames(strK6.results)=c(1,2,3,4,5,6)
strK6.results$pop=as.character(strK6[1:300,4])
for (i in 1:nrow(strK6.results)) {
  strK6.results$pop[i]=as.character(labelPops[which(strK6.results$pop[i] == labelPops$Location_index),1])
}
strK6.results$indNames=names(native$tab[,1])
strK6.results=melt(strK6.results)

colnames(strK6.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK6.results$Original_Population=factor(strK6.results$Original_Population, levels = c("Jap","S4","S6","S3","S1","S2","S16","S13","S14","S15","S17","S9","S18","S19","S20","S11","S10","Tib"))

# Plotting the barplot
strK6_barplot=ggplot(data=strK6.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  scale_fill_manual(values = cols.native) +
  facet_grid(~Original_Population, scales = "free") +
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
        strip.text=element_text(size=28, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
strK6_barplot

ggsave("Native_barplot_K6b.tiff",
       device="tiff",dpi=320,units="cm",width=60,height=9)

###################################
# K=13
setwd("~/Documents/2b_Native_tmp/Harvester")
strK13=read.table("K13.indfile")
setwd("~/Documents/2b_Native_tmp/")
# Get the same colors and cluster number as DAPC
strK13.results=as.data.frame(strK13[1:300,(c(1:13)+5)]) # population clusters in the good order for colors

colnames(strK13.results)=c(1:13)
strK13.results$pop=as.character(strK13[1:300,4])
for (i in 1:nrow(strK13.results)) {
  strK13.results$pop[i]=as.character(labelPops[which(strK13.results$pop[i] == labelPops$Location_index),1])
}
strK13.results$indNames=names(native$tab[,1])
strK13.results=melt(strK13.results)

colnames(strK13.results)=c("Original_Population","Sample","Assigned_Population","Posterior_membership_probability")
######### GGPLOT
# Cluster populations together in reference demes before plotting...
strK13.results$Original_Population=factor(strK13.results$Original_Population, levels = c("Jap","S4","S6","S3","S1","S2","S16","S13","S14","S15","S17","S9","S18","S19","S20","S11","S10","Tib"))

# Plotting the barplot
strK13_barplot=ggplot(data=strK13.results, aes(x=Sample, y=Posterior_membership_probability, fill=Assigned_Population))+
  geom_bar(stat='identity',width=1) +
  facet_grid(~Original_Population, scales = "free") +
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
        strip.text=element_text(size=28, colour="black"),
        legend.key = element_rect(fill = "white", size = 1),
        legend.text=element_text(size=28),
        legend.title=element_text(size=28))
strK13_barplot

ggsave("Native_barplot_K13b.tiff",
       device="tiff",dpi=320,units="cm",width=60,height=9)
