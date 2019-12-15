
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
immune_variable = "LF" #"mean_centered_ICRscore"
ICR_pancancer= "ICR_pancancer"

# Load data
#mutations = read.csv("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/Revision 2 submission/CTNNB1_EXON3 mutations_V2.csv",
                    # stringsAsFactors = FALSE)
load("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/4. PTEN and CTTNB1/PTEN_CTNNB1_MUT.RData")
all_samples = read.csv("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/Revision 2 submission/Additional materials related t CTNNB1/all_samples_CTNNB1_analyzed.csv",
                       stringsAsFactors = FALSE)
load("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers_Skipped_DLBC.Rdata")
load("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Clustering/5.3.1.and.5.8.Hallmark_and_ICR_cluster_assignment_allcancers_ICRbenefit_cluster.Rdata")
Immunity = read.csv("./3_DataProcessing/External/mmc2-PanImmune_MS.csv", stringsAsFactors = FALSE)
Immunity$LF_MC = NA

cancer.types = as.character(unique(Immunity$TCGA.Study))

for (j in cancer.types) {
  cancer.scaled = as.numeric(scale(Immunity$Leukocyte.Fraction[Immunity$TCGA.Study == j],center = TRUE))
  Immunity$LF_MC[Immunity$TCGA.Study==j] = cancer.scaled
}

#load("./3_DataProcessing/Assembler_Panca_Normalized_filtered/Pancancer/Pancancer_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered_for_Rosalyn.Rdata")
mean_centered_ICR = read.csv("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/4. PTEN and CTTNB1/Mean_centered_per_cancer_ICRscore.csv",
                             stringsAsFactors = FALSE)

#RNASeqLOG2 = log(filtered.norm.RNAseqData.all + 1, 2)
#t_RNASeqLOG2 = as.data.frame(t(RNASeqLOG2))

# Analysis
mutations = subset_mut[which(subset_mut$Hugo_Symbol == "CTNNB1"),]

df_plot = data.frame(Sample = rownames(ICR_cluster_assignment_allcancers), Cancer = NA, mutation = NA)
#df_plot = df_plot[which(substring(df_plot$Sample, 1, 12) %in% copy_number_alterations$Patient.ID),]

df_plot$Cancer = ICR_cluster_assignment_allcancers$Cancer[match(substring(df_plot$Sample, 1, 12), substring(rownames(ICR_cluster_assignment_allcancers), 1, 12))]
table(df_plot$Cancer)
df_plot$ICRscore = ICR_cluster_assignment_allcancers$ICRscore[match(substring(df_plot$Sample, 1, 12), substring(rownames(ICR_cluster_assignment_allcancers), 1, 12))]
df_plot$mean_centered_ICRscore = mean_centered_ICR$ICRscore.MC[match(substring(df_plot$Sample, 1, 12), 
                                                                     substring(mean_centered_ICR$TCGA_Sample_Barcode, 1, 12))]
df_plot$LF= Immunity$Leukocyte.Fraction[match(substring(df_plot$Sample, 1, 12), Immunity$TCGA.Participant.Barcode)]
df_plot$LF_MC = Immunity$LF_MC[match(substring(df_plot$Sample, 1, 12), Immunity$TCGA.Participant.Barcode)]

if(ICR_pancancer == "ICR_pancancer"){
  df_plot$ICR_cluster = Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster[match(substring(df_plot$Sample, 1, 12), 
                                                                                         substring(rownames(Hallmark_and_ICR_cluster_assignment_allcancers), 1, 12))]
}else{
  df_plot$ICR_cluster = ICR_cluster_assignment_allcancers$HML_cluster[match(substring(df_plot$Sample, 1, 12), substring(rownames(ICR_cluster_assignment_allcancers), 1, 12))]
}


#df_plot$IFNG = t_RNASeqLOG2$IFNG[match(substring(rownames(df_plot), 1, 12), substring(rownames(t_RNASeqLOG2), 1, 12))]
#df_plot$GZMB = t_RNASeqLOG2$GZMB[match(substring(rownames(df_plot), 1, 12), substring(rownames(t_RNASeqLOG2), 1, 12))]

df_plot$mutation = "non-mutated"
df_plot$mutation[which(substring(df_plot$Sample, 1, 12) %in% substring(mutations$Tumor_Sample_Barcode, 1, 12))] = "mutated"

df_plot$mutation = factor(df_plot$mutation, levels = c("non-mutated", "mutated"))
table(df_plot$mutation, exclude = NULL)

df_plot = df_plot[which(substring(df_plot$Sample, 1, 12) %in% all_samples$Patient.ID),]
  
if(immune_variable %in% c("ICRscore", "mean_centered_ICRscore", "LF", "LF_MC")){
  immune_lab = gsub("_", " ", immune_variable)
}else{
  immune_lab = paste0("log2 transformed ", immune_variable, " (offset +1)")
}

plot = ggplot(df_plot, aes(x= mutation, y=get(immune_variable))) +
  geom_boxplot(outlier.stroke = 0.2, aes(color = mutation)) +
  scale_color_manual(values = c("non-mutated" = "#01C03C", "mutated" = "#FF3800")) +
  theme_bw() +
  facet_grid(cols = vars(Cancer)) +
  theme(axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 15),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", size = 13)) +
  xlab(paste0("CTNNB1 mutation", " ")) +
  stat_compare_means(method = "t.test", comparisons = list(c("non-mutated", "mutated"))) +
  ylab(immune_lab)

dir.create("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/4. PTEN and CTTNB1/Figures", 
           showWarnings = FALSE)
png(filename = paste0("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/4. PTEN and CTTNB1/Figures/v2_Mut_CTNNB1_",
                       "_", immune_variable, "_split_by_cancer.png"), res = 600, width =20, height = 5, units = "in")
plot(plot)
dev.off()


plot = ggplot(df_plot, aes(x= mutation, y=get(immune_variable))) +
  geom_boxplot(outlier.stroke = 0.2, aes(color = mutation)) +
  scale_color_manual(values = c("non-mutated" = "#01C03C", "mutated" = "#FF3800")) +
  theme_bw() +
  #facet_grid(cols = vars(Cancer)) +
  theme(axis.title.x = element_text(colour = "black", size = 15),
        axis.title.y = element_text(colour = "black", size = 15),
        axis.text.x = element_blank(),
        axis.text.y = element_text(colour = "black", size = 15),
        legend.title = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(colour = "black", size = 13)) +
  xlab(paste0("CTNNB1 mutation", " ")) +
  stat_compare_means(method = "t.test", comparisons = list(c("non-mutated", "mutated"))) +
  ylab(immune_lab)

dir.create("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/4. PTEN and CTTNB1/Figures", 
           showWarnings = FALSE)
png(filename = paste0("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/4. PTEN and CTTNB1/Figures/v2_Mut_CTNNB1_",
                      "_", immune_variable, "_pancancer.png"), res = 600, width =3.2, height = 3.8, units = "in")
plot(plot)
dev.off()

# investigation

test = data.frame(table(df_plot$Cancer, df_plot$mutation, df_plot$ICR_cluster))
test = data.frame(table(df_plot$mutation, df_plot$ICR_cluster))
excluded_mut = mutations[-which(substring(mutations$Sample.ID, 1, 12) %in% substring(rownames(ICR_cluster_assignment_allcancers), 1, 12)),]
#write.csv(test, file = "~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/Revision 2 submission/Additional materials related t CTNNB1/mut_per_cancer_overview.csv")
write.csv(test, file = paste0("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/Revision 2 submission/Additional materials related t CTNNB1/v2_mut_pancancer_overview_", ICR_pancancer, ".csv"))

