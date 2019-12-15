
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/") 

#install.packages(c("Rtsne", "ggplot2"))
library("Rtsne")
library("ggplot2")

load("./3_DataProcessing/Assembler_Panca_Normalized_filtered/Pancancer/Pancancer_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered_for_Rosalyn.Rdata")
load("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Clustering/5.3.Hallmark_and_ICR_cluster_assignment_allcancers_SkippedDLBC.Rdata")
load("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Signature_Enrichment/ssGSEA_Selected.pathways_ES.Rdata")

RNASeq_QN = t(filtered.norm.RNAseqData.all)
RNASeq_QN_LOG2 = log(RNASeq_QN +1, 2)


tsne_model_1 = Rtsne(t(ES.all), check_duplicates=FALSE, pca=TRUE, perplexity=15, theta=0.5, dims=2)
tsne_model_2 = Rtsne(RNASeq_QN_LOG2, check_duplicates=FALSE, pca=TRUE, perplexity=15, theta=0.5, dims=2)

load("~/Desktop/tsne_model2.Rdata")

df = data.frame(t(ES.all))
df$Cancer = Hallmark_and_ICR_cluster_assignment_allcancers$Cancer[match(rownames(df), rownames(Hallmark_and_ICR_cluster_assignment_allcancers))]
#df$Type = "Other cancer types"
df$Type = as.character(df$Cancer)
df$Type[which(df$Cancer == "SKCM" & substring(rownames(df), 14, 15) == "06")] = "SKCM_Meta"
df$Type[which(df$Cancer == "SKCM" & substring(rownames(df), 14, 15) == "01")] = "SKCM_Primary"
df$Type[which(df$Cancer == "UVM")] = "UVM"
df = df[rownames(RNASeq_QN_LOG2),]

colors = c("grey", "orange", "blue", "black")
names(colors) = unique(df$Type)


plot(tsne_model_1$Y, t='n', main="tSNE", xlab="tSNE dimension 1", ylab="tSNE dimension 2", "cex.main"=2, "cex.lab"=1.5)
text(tsne_model_1$Y, labels=df$Type, col=colors[df$Type])

d_tsne_1 = as.data.frame(tsne_model_2$Y) 
d_tsne_2 <- d_tsne_1

dir.create("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/2_Figures", showWarnings = FALSE)
png("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/2_Figures/tSNE_6.png", height = 3.8, width = 4.2,
    units = "in", res = 600)
plot= ggplot(d_tsne_2, aes(x=V1, y=V2,color=df$Type)) +
  geom_point(aes(x=V1, y=V2,fill=df$Type),size=0.2) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  #scale_colour_brewer(palette = "Set2")+
  scale_color_manual(values=c("grey", "orange", "blue", "darkgreen"))+
  xlab("tSNE dimension 1") + ylab("tSNE dimension 2") +
  #ggtitle("Perplexity=15") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none")
plot(plot)
dev.off()

load("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Cancer_color_table.Rdata")
rownames(Cancer_color_table) = Cancer_color_table$Group.1
Cancer_color_table = Cancer_color_table[unique(df$Type),]
Cancer_color_table[24,1] = "SKCM_Primary"
Cancer_color_table[24,2] = "orange"
Cancer_color_table[25,1] = "SKCM_Meta"
Cancer_color_table[25,2] = "blue"

png("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/2_Figures/tSNE_all_cancers_colored2.png", height = 3.8, width = 4.2,
    units = "in", res = 600)
plot= ggplot(d_tsne_2, aes(x=V1, y=V2,color=df$Type)) +
  geom_point(aes(x=V1, y=V2,fill=df$Type),size=0.2) +
  guides(colour=guide_legend(override.aes=list(size=6))) +
  #scale_colour_brewer(palette = "Set2")+
  scale_color_manual(values=c(Cancer_color_table$color))+
  xlab("tSNE dimension 1") + ylab("tSNE dimension 2") +
  #ggtitle("Perplexity=15") +
  theme_light(base_size=20) +
  theme(axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        legend.position = "none")
plot(plot)
dev.off()
