#################################################################
###
###
### Data input:
### "./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
### "_Bindea_xCell_HallmarkPathways.Rdata")
### Output :
### "./5_Figures/Correlation_plots/ICR_Correlation_plots/", download.method, 
### "/ICR_Correlation_plot_",Cancer,".png"
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path, "R tools/heatmap.3.R"))

required.packages <- c("corrplot", "stringr")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
CancerTYPES = "ALL"                                                                                                      # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("DLBC")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "Assembler_Panca_Normalized_filtered"                                                                                       # Specify download method (this information to be used when saving the file)
my.palette = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 297)                                                 # Specify which genes will be correlated
Log_file = paste0("./1_Log_Files/", download.method, "/5.1_Pancancer_Correlation_matrix_Signatures/5.1_Pancancer_Correlation_matrix_signatures", 
                  "_Bindea_xCell_Hallmark", "_Log_File_", gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"
test = "pearson"
display_correlations = "irrespective_of_significance"                                                                    # Can either be "only_significant" or "irrespective_of_significance"
IPA_excluded = "all"
Source_surv_data = "Cell_paper"
Cutoff_HR = 1

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.
# Define parameters (based on loaded data)

if (CancerTYPES == "ALL") { 
  CancerTYPES = c("Pancancer", TCGA.cancersets$cancerType)
}

# Create folders
dir.create("./5_Figures/",showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method), showWarnings = FALSE)

dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method, "/5.1_Pancancer_Correlation_matrix_Signatures"), showWarnings = FALSE)
cat("This is a log file for creating correlation plots",
    "__________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Parameters Used :",
    paste0("CancerTYPES = ", CancerTYPES),
    paste0("Cancer_skip = ", Cancer_skip),
    paste0("download.method = ", download.method),
    "",
    "Scripts output :",
    file = Log_file,
    append = FALSE, sep= "\n")

N.sets = length(CancerTYPES)

start.time <- Sys.time ()
load(file = paste0("./4_Analysis/", download.method, "/ACC/Correlation/Correlation_matrix Selected.pathways_pearson_ACC.Rdata")) 

TCGA.cancersets = rbind(TCGA.cancersets, c("Pancancer", "Pancancer"))
row.names(TCGA.cancersets) = TCGA.cancersets$cancerType
pancancer_Geneset_cor_table = t(TCGA.cancersets)[-c(1,2),]

pancancer_Geneset_cor_table = rbind(pancancer_Geneset_cor_table,matrix(nrow = nrow(Geneset_cor),ncol=N.sets))
rownames(pancancer_Geneset_cor_table) = rownames(Geneset_cor)

for (i in 1:N.sets) {
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(Cancer == "LAML") 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}
  
  load(paste0("./4_Analysis/",download.method, "/", Cancer, "/Correlation/", "Correlation_matrix Selected.pathways_pearson_", Cancer, ".Rdata"))
  if(display_correlations == "only_significant"){
    N.columns = ncol(Geneset_cor)
    N.rows = nrow(Geneset_cor)
    for (i in 1:N.columns){
      for (j in 1:N.rows){
        Geneset_cor[i,j] = ifelse(Geneset_cor_sign[[1]][i,j] <0.05 | Geneset_cor_sign[[1]][i,j] >0.95, Geneset_cor[i,j], 0)
      }
    }
  }
  pancancer_Geneset_cor_table[, Cancer] = as.numeric(Geneset_cor[,"ICR_score"])
}
  
# convert to numeric matrix
mode(pancancer_Geneset_cor_table) = "numeric"

# Remove LAML and DLBC from all correlation tables
pancancer_Geneset_cor_table = pancancer_Geneset_cor_table[,-c(7, 14)]

if(IPA_excluded == "IPA_excluded"){
  pancancer_Geneset_cor_table = pancancer_Geneset_cor_table[grep(pattern = "IPA]", rownames(pancancer_Geneset_cor_table), invert = TRUE),]
}

dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Oncogenic_pathways_ICR_correlation/"), showWarnings = FALSE)

rownames_order = rownames(pancancer_Geneset_cor_table)[order(rowMeans(pancancer_Geneset_cor_table),decreasing = TRUE)]

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/", Source_surv_data ,"_Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype %in% c("LAML", "DLBC")),]
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])
All_survival_analysis_data = All_survival_analysis_data[order(All_survival_analysis_data$HR, decreasing = TRUE),]
Enabled = All_survival_analysis_data[which(All_survival_analysis_data$Cancer %in% ICR_enabled_cancers),]
Disabled = All_survival_analysis_data[which(All_survival_analysis_data$Cancer %in% ICR_disabled_cancers),]
Neutral = All_survival_analysis_data[which(All_survival_analysis_data$Cancer %in% ICR_neutral_cancers),]
Cancer_order = c("Pancancer", as.character(Enabled$Cancertype), as.character(Neutral$Cancertype), as.character(Disabled$Cancertype))

pancancer_Geneset_cor_table = pancancer_Geneset_cor_table[rownames_order,Cancer_order]
pancancer_Geneset_cor_table = pancancer_Geneset_cor_table[-which(rownames(pancancer_Geneset_cor_table) == "ICR_score"),]

rownames(pancancer_Geneset_cor_table) = gsub(".*] ", "",rownames(pancancer_Geneset_cor_table))
rownames(pancancer_Geneset_cor_table)[which(rownames(pancancer_Geneset_cor_table) == "Phopholipase")] = "Phospholipase"

## Correlation heatmap
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Oncogenic_pathways_ICR_correlation/v4_5.1.Oncogenic_pathways_ICR_correlation_", IPA_excluded, "_", display_correlations, "_for_paper.png"),res=600,height= 19,width=25,unit="in")
heatmap.3 (pancancer_Geneset_cor_table,
           col= my.palette,
           main = paste0("Pan-Cancer ", test, " correlation \n between ICR and oncogenic pathway gene signatures "),
           cex.main = 1,
           cexCol = 1.6,
           cexRow = 1.5,
           Rowv = NULL,
           Colv = NULL,
           margins=c(14, 30))

#title(sub = list(paste0("Correlations between ICR and gene signatures were calculated with R package corrplot. \n", 
                        #gsub("_", " ", str_to_title(display_correlations)), "."), cex = 1), outer = FALSE, line = -1)
dev.off()

# Color Key
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Oncogenic_pathways_ICR_correlation/5.1.Color_Key_Oncogenic_pathways_ICR_correlation_", IPA_excluded, "_", display_correlations, ".png"),res=600,height= 7,width=7,unit="in")
heatmap.3 (pancancer_Geneset_cor_table,
           col= my.palette,
           main = paste0("Pan-Cancer ", test, " correlation \n between ICR and oncogenic pathway gene signatures "),
           cex.main = 1,
           cexCol = 1.6,
           cexRow = 1.5,
           Rowv = FALSE,
           margins=c(14, 30))

title(sub = list(paste0("Correlations between ICR and gene signatures were calculated with R package corrplot. \n", 
                        gsub("_", " ", str_to_title(display_correlations)), "."), cex = 1), outer = FALSE, line = -1)
dev.off()

dir.create(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Correlation"), showWarnings = FALSE)
save(pancancer_Geneset_cor_table,
     file = paste0("./4_Analysis/", download.method, "/Pan_Cancer/Correlation/", "Correlation_Oncogenic_pathways_and_ICR_", display_correlations, ".Rdata"))
end.time <- Sys.time ()
time <- end.time - start.time
print (paste0("Between start script and completion script: ", time))
