
## MSI investigation

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path,"R tools/heatmap.3.R"))

required.packages <- c("gtools", "circlize", "dendsort", "ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
download.method = "Assembler_Panca_Normalized_filtered"
Cancer = "COAD"

# Load data
MSI_data = read.csv("./2_Data/MSI_data_NatComms/ncomms_MSI.csv", stringsAsFactors = FALSE)
load("./3_DataProcessing/External/BinaryMatrix.From.Michele.RData")
load(paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".Assembler_Panca_Normalized_filtered.EDASeq.ICR.reps5000/", 
            Cancer, "_ICR_cluster_assignment_k2-6.Rdata"))
Michele_genes = read.csv("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Mutational_load/Michele_genes_v2_ICRscore.csv",
                         stringsAsFactors = FALSE, header = FALSE)
Michele_genes = c(Michele_genes[,1])
Michele_genes_neg_coefficient = Michele_genes[1:35]

# Analysis
patients = substring(rownames(table_cluster_assignment), 1, 16)
binary_matrix = binary_matrix[which(substring(rownames(binary_matrix), 1, 16) %in% patients), which(colnames(binary_matrix) %in% Michele_genes_neg_coefficient)]
plot_df = data.frame(patient_barcode = substring(rownames(binary_matrix), 1, 12), MSI_status = NA, number_mut_genes = NA, ICR_score = NA)
plot_df$patient_barcode = as.character(plot_df$patient_barcode)
plot_df$number_mutated_neg_coef_genes = rowSums(binary_matrix)
plot_df$MSI_status = MSI_data$MSI_category_nb_from_TCGA_consortium[match(plot_df$patient_barcode, MSI_data$Barcode)]
plot_df$ICR_score = table_cluster_assignment$ICRscore[match(plot_df$patient_barcode, substring(rownames(table_cluster_assignment), 1, 12))]

plot_df = plot_df[-which(is.na(plot_df$MSI_status)),]

#plot = ggplot(plot_df, aes(x = MSI_status, y = ICR_score)) +
  #geom_boxplot(outlier.shape = NA) +
  #theme_bw()

my_comparisons = list(c("msi-h", "msi-l"), c("msi-h", "mss"), c("msi-l", "mss"))

plot = ggboxplot(plot_df, x = "MSI_status", y = "ICR_score", add = "jitter", outlier.shape = NA, 
                 order = c("msi-h", "msi-l", "mss")) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  ggtitle(Cancer) +
  theme(plot.title = element_text(hjust = 0.5))

dir.create(paste0("./4_Analysis/", download.method, "/", Cancer, "/Mutations"), showWarnings = FALSE)
dir.create(paste0("./4_Analysis/", download.method, "/", Cancer, "/Mutations/MSI"), showWarnings = FALSE)
png(paste0("./4_Analysis/", download.method, "/", Cancer, "/Mutations/MSI/Boxplot_ICR_score_by_MSI_status.png"), res = 600,
    width = 6, height = 6, units = "in")
plot(plot)
dev.off()

plot = ggboxplot(plot_df, x = "MSI_status", y = "number_mutated_neg_coef_genes", add = "jitter", outlier.shape = NA, 
                 order = c("msi-h", "msi-l", "mss")) +
  stat_compare_means(comparisons = my_comparisons, method = "t.test") +
  ggtitle(Cancer) +
  theme(plot.title = element_text(hjust = 0.5))

png(paste0("./4_Analysis/", download.method, "/", Cancer, "/Mutations/MSI/Boxplot_number_muts_neg_coefficient_by_MSI_status.png"), res = 600,
    width = 6, height = 6, units = "in")
plot(plot)
dev.off()

