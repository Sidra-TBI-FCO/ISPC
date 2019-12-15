# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("ggplot2")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = "DLBC"                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized_filtered"                                                                                      # Specify download method "TCGA_Assembler", "Pancancer_matrix" (this information to be used when saving the file)
assay.platform = "gene_RNAseq"
ICR_classification_k = "HML_classification"                                                                                                         # Subset can be "ICR High", "ICR Low", "All"
order = "Enabled_disabled"
Source_surv_data = "Cell_paper"
Cutoff_HR = 1

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.
load("./4_Analysis/Pancancer_matrix/Pan_Cancer/Clustering/Cancer_color_table.Rdata")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

load(paste0("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Clustering/3.17.Hallmark_and_ICR_cluster_assignment_allcancers_SkippedDLBC.Rdata"))
dir.create("./5_Figures", showWarnings = FALSE)
dir.create("./5_Figures/ICR_distribution_plots", showWarnings = FALSE)
dir.create(paste0("./5_Figures/ICR_distribution_plots/ICR_histogram_PER_CANCER_ICR_Clusters_across_cancers/"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/ICR_distribution_plots/ICR_histogram_PER_CANCER_ICR_Clusters_across_cancers/", download.method), showWarnings = FALSE)

Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster = as.factor(Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster)
Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster = relevel(Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster, "ICR High")
Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster = relevel(Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster, "ICR Medium")
Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster = relevel(Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster, "ICR Low")

load(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis/", Source_surv_data, "_Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))

ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])

if(order == "Enabled_disabled"){
  #All_survival_analysis_data = All_survival_analysis_data[order(All_survival_analysis_data$HR, decreasing = TRUE),]
  #Cancer_order = as.character(All_survival_analysis_data$Cancertype[-which(All_survival_analysis_data$Cancertype == "LAML")])
  ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
  ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
  ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])
  All_survival_analysis_data = All_survival_analysis_data[order(All_survival_analysis_data$HR, decreasing = TRUE),]
  Enabled = All_survival_analysis_data[which(All_survival_analysis_data$Cancer %in% ICR_enabled_cancers),]
  Disabled = All_survival_analysis_data[which(All_survival_analysis_data$Cancer %in% ICR_disabled_cancers),]
  Neutral = All_survival_analysis_data[which(All_survival_analysis_data$Cancer %in% ICR_neutral_cancers),]
  
  Cancer_order = c(as.character(Enabled$Cancertype), as.character(Neutral$Cancertype), as.character(Disabled$Cancertype))
  Hallmark_and_ICR_cluster_assignment_allcancers$Cancer = factor(Hallmark_and_ICR_cluster_assignment_allcancers$Cancer, levels = Cancer_order)
}

png(paste0("./5_Figures/ICR_distribution_plots/ICR_histogram_PER_CANCER_ICR_Clusters_across_cancers/", download.method,
           "/Barplot_", order, "_ICR_clusters_Pancancer.png"),res=600,height=8,width=30,unit="in")

plot = ggplot(Hallmark_and_ICR_cluster_assignment_allcancers, aes(Cancer)) + geom_bar(aes(fill=HML_cluster)) + labs(fill = NULL) +
  scale_fill_manual("ICR clusters", values = c("blue", "green", "red")) +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), panel.background = element_rect(fill = "white", colour = "grey"), 
        axis.text = element_text(size = 17, colour = "black"),
        axis.title = element_text(size = 22, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 22, colour = "black")) 

plot(plot)
dev.off()

library(dplyr)
H2 <- Hallmark_and_ICR_cluster_assignment_allcancers %>% 
  group_by(Cancer,HML_cluster) %>% 
  summarise(count=n()) %>% 
  mutate(perc=count/sum(count))

png(paste0("./5_Figures/ICR_distribution_plots/ICR_histogram_PER_CANCER_ICR_Clusters_across_cancers/", download.method,
           "/Barplot_v2_", order, "_ICR_clusters_Pancancer.png"),res=600,height=8,width=30,unit="in")
plot = ggplot(H2, aes(x = Cancer, y = perc*100, fill = HML_cluster)) + geom_bar(stat="identity") +
  scale_fill_manual("ICR clusters", values = c("blue", "green", "red")) +
  labs(x = "Cancer", y = "percent", fill = "ICR clusters") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), panel.background = element_rect(fill = "white", colour = "grey"), 
        axis.text = element_text(size = 17, colour = "black"),
        axis.title = element_text(size = 22, colour = "black"),
        legend.title = element_text(size = 10, colour = "black"),
        legend.text = element_text(size = 22, colour = "black")) 
plot(plot)
dev.off()
