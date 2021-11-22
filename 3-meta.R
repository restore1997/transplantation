rm(list = ls())
library(ggplot2)
library(ggpubr)

load("1-meta.rda")
pred <- meta[,c("data","status","seurat_clusters")]
SAMPLE <- as.data.frame(unique(pred$data))


#################SAMPLE
##L1
L1 <- pred[pred$data=="L1",]
L1 <- as.data.frame(table(L1$seurat_clusters))
L1$percentage <- L1$Freq/sum(L1$Freq)
L1$status <- "transplantation"
L1$data <- "L1"
L1 <- L1[,-2]
colnames(L1) <- c("cluster","percentage","status","data")

##L2
L2 <- pred[pred$data=="L2",]
L2 <- as.data.frame(table(L2$seurat_clusters))
L2$percentage <- L2$Freq/sum(L2$Freq)
L2$status <- "transplantation"
L2$data <- "L2"
L2 <- L2[,-2]
colnames(L2) <- c("cluster","percentage","status","data")


##L4
L4 <- pred[pred$data=="L4",]
L4 <- as.data.frame(table(L4$seurat_clusters))
L4$percentage <- L4$Freq/sum(L4$Freq)
L4$status <- "transplantation"
L4$data <- "L4"
L4 <- L4[,-2]
colnames(L4) <- c("cluster","percentage","status","data")


##L5
L5 <- pred[pred$data=="L5",]
L5 <- as.data.frame(table(L5$seurat_clusters))
L5$percentage <- L5$Freq/sum(L5$Freq)
L5$status <- "transplantation"
L5$data <- "L5"
L5 <- L5[,-2]
colnames(L5) <- c("cluster","percentage","status","data")

##N1
N1 <- pred[pred$data=="N1",]
N1 <- as.data.frame(table(N1$seurat_clusters))
N1$percentage <- N1$Freq/sum(N1$Freq)
N1$status <- "health"
N1$data <- "N1"
N1 <- N1[,-2]
colnames(N1) <- c("cluster","percentage","status","data")

##N2
N2 <- pred[pred$data=="N2",]
N2 <- as.data.frame(table(N2$seurat_clusters))
N2$percentage <- N2$Freq/sum(N2$Freq)
N2$status <- "health"
N2$data <- "N2"
N2 <- N2[,-2]
colnames(N2) <- c("cluster","percentage","status","data")

##N3
N3 <- pred[pred$data=="N3",]
N3 <- as.data.frame(table(N3$seurat_clusters))
N3$percentage <- N3$Freq/sum(N3$Freq)
N3$status <- "health"
N3$data <- "N3"
N3 <- N3[,-2]
colnames(N3) <- c("cluster","percentage","status","data")


proportion_total <- rbind(
  L1,L2,L4,L5,N1,N2,N3               
)


ggplot(proportion_total, aes(fill=status, y=percentage, x=cluster)) + 
  geom_bar(position="dodge",  stat="summary",fun="mean")+
  scale_fill_manual(values = c("health"="#7A7269","transplantation"="#b3003c"))+
  geom_point(aes(shape=status),position = position_dodge(0.9))+
  theme_classic()+    #theme_bw()
  stat_compare_means(aes(group = status),label = "p.format",method = "wilcox.test") ##wilcox.test
