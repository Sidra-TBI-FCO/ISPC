#################################################################
###
###
### Data input:
###
### Output :
### 
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages <- c("ggplot2")
ipak(required.packages)

# Set Parameters
download.method = "Assembler_Panca_Normalized_filtered"                                                                                       # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"
signature = c("TIS")


# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers_Skipped_DLBC.Rdata"))
load(paste0("./3_DataProcessing/", download.method,"/Pancancer/Pancancer_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata"))

if(signature == "TIS"){
  load(paste0(code_path, "Datalists/18_gene_TIS.Rdata"))
  selected_genes = TIS_genes
}else{selected_genes = signature}

subset_RNASeq = filtered.norm.RNAseqData.all[selected_genes,rownames(ICR_cluster_assignment_allcancers)] # Remove DLBC patients
subset_RNASeq = t(subset_RNASeq)
selected_subset_RNAseq_log2 = log(subset_RNASeq +1, 2)
ICR_cluster_assignment_allcancers$TIS_score = rowMeans(selected_subset_RNAseq_log2) 

dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Signature_Correlation_Plots/"), showWarnings = FALSE)
png(paste0("./5_Figures/Pancancer_plots/", download.method,"/Signature_Correlation_Plots/", "TIS","_ICR_Correlation_plot.png"),
    res = 600, height = 4, width = 4, units = "in")
ggplot(ICR_cluster_assignment_allcancers, aes(x = ICRscore, y = `TIS`)) +
  geom_point(size = 0.5, shape = 1) +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  geom_smooth(method=lm, color="blue")

dev.off()

#dev.new()
#plot(ICR_cluster_assignment_allcancers$ICRscore, ICR_cluster_assignment_allcancers$TISscore)
cor.test(x = ICR_cluster_assignment_allcancers$ICRscore, y = ICR_cluster_assignment_allcancers$`TIS`, method = "pearson")

