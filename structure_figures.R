# LIBRARIES
if (!require("ggplot2", quietly = TRUE)) install.packages("ggplot2")
if (!require("reshape2", quietly = TRUE)) install.packages("reshape2")
if (!require("ggsci", quietly = TRUE)) install.packages("ggsci")
library(ggplot2)
library(reshape2)
library(ggsci)

library(tidyverse)
library(dplyr)
library(stringr)

################################################################################
# ADMIXTURE - cross-validation error
################################################################################

DATASET="chr22_pop_dist_pruned.admixture.cv.error"
DATASET="chr22_pop_finer_pruned.admixture.cv.error"

CV <- read.table(DATASET)
colnames(CV) <- c("K","CV")
graph_cv <- ggplot(CV,aes(x=K,y=CV)) + 
  geom_line() + geom_point() +
  scale_x_continuous(breaks = c(1:max(CV$K))) +
  labs(title="Cross-validation plot", x="K", y="Cross-validation error") +
  theme(axis.text.y = element_text(colour="black",size=12),
        axis.text.x = element_text(colour="black",size=12),
        axis.title = element_text(size=16,colour="black",face="bold"),
        plot.title = element_text(size=18,colour="black",face="bold")) 
ggsave(paste0(DATASET,".jpg"), graph_cv, width=4.2, height=4.5, dpi=600)

################################################################################
# ADMIXTURE - admixture analysis
################################################################################

DATASET="chr22_pop_dist_pruned.4.Q"
POPLIST="poplist_human_pop_dist.txt"
  
DATASET="chr22_pop_finer_pruned.4.Q"
POPLIST="poplist_human_pop_finer.txt"
  
# Read data
admixture <- read.table(DATASET)

# Table with individual names and population
id <- read.table(POPLIST,header=FALSE)
colnames(id) <- c("V1","Individual","Population")
admixture <- cbind(id[-1],admixture)
colnames(admixture) <- gsub("V", "K", colnames(admixture))
admixture <- data.frame(POP_IND=paste(admixture$Population,admixture$Individual,sep="_"),admixture)

# Transform the admixture object into a long format
admixture_long <- melt(admixture,id.vars=c("POP_IND","Individual","Population"),variable.name="ANCESTRY",value.name="PERC")

# Plot
admixture_long$POP_IND_s <- factor(admixture_long$POP_IND, levels = unique(admixture_long$POP_IND))
graph_admixture <- ggplot(admixture_long,aes(x=POP_IND_s,y=PERC,fill=ANCESTRY)) +
  geom_bar(stat="identity",width=0.95) +
  scale_y_continuous(expand = c(0,0)) +
  scale_x_discrete(expand = c(0,0)) +
  labs(x="Sample", y="Ancestry") +
  scale_fill_npg() + 
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.title.x = element_text(colour="black",size=14,face="bold"),
        axis.text.x = element_text(colour="black",size=6,face="bold", angle=90, hjust=0, vjust=0.5),
        axis.title.y = element_text(colour="black",size=14,face="bold"),
        axis.text.y = element_text(colour="black",size=14),
        legend.position = "none") 
ggsave(paste0(DATASET,".jpg"), graph_admixture, width=14, height=4, dpi=600)










################################
# PCA

# Read data
pca <- read_table("chr22_pop_dist_pruned_pca.eigenvec")
eigenval <- scan("chr22_pop_dist_pruned_pca.eigenval")

# Clean data
pca <- pca[,-1]
names(pca)[1] <- "Individual"

# Assign population info
id <- read.table("poplist_human_pop_dist.txt",header=FALSE)
colnames(id) <- c("V1","Individual","Population")
pca <- cbind(id[3],pca)

# Plot - percentage of variance each principal component explains
pve <- data.frame(PC = 1:ncol(pca[-c(1,2)]), pve = eigenval/sum(eigenval)*100)
pve_plot <- ggplot(pve, aes(PC, pve)) + 
  geom_bar(stat = "identity") + 
  labs(title="Variance explained by each PC",y="Percentage") + 
  scale_x_continuous(breaks = c(0:ncol(pca[-c(1,2)]))) +
  theme(axis.text.y = element_text(colour="black",size=12),
        axis.text.x = element_text(colour="black",size=12),
        axis.title = element_text(size=16,colour="black",face="bold"),
        plot.title = element_text(size=18,colour="black",face="bold"))

ggsave("Percentage-variance-explained.jpg",pve_plot,width=5.5,height=4.5,dpi=600)

cumsum(pve$pve)

# Plot - PCA
pop <- c("YRI","LWK","GBR","GIH","CHB","PEL")
pop <- c("GBR","IBS","TSI","GIH")
pca$Population <- factor(pca$Population, levels = pop)
pca_plot <- ggplot(pca, aes(x=PC1, y=PC2, fill = Population)) + 
  geom_hline(yintercept=0,size=0.2,linetype=5,color="gray20") +
  geom_vline(xintercept=0,size=0.2,linetype=5,color="gray20") +
  geom_point(size = 3.5, alpha=0.7, shape=21, color="gray40",stroke=0.2) +
  scale_fill_npg() +
  xlab(paste0("PC1 (", signif(pve$pve[1], 3), "%)")) + 
  ylab(paste0("PC2 (", signif(pve$pve[2], 3), "%)")) +
  theme_bw() +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        axis.text.y = element_text(colour="black",size=12),
        axis.text.x = element_text(colour="black",size=12),
        axis.title = element_text(size=16,colour="black",face="bold"),
        legend.text = element_text(colour="black",size=10),
        legend.title = element_text(colour="black",size=12, face="bold"),
        plot.title = element_text(size=18,colour="black",face="bold"))
ggsave("PCA.jpg",pca_plot,width=5.5,height=4.5,dpi=600)




