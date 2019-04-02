#################################################################
###
###
### Data input:
### "./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
### "_Bindea_xCell_HallmarkPathways.Rdata")
### Output :
### "./5_Figures/Correlation_plots/ICR_Correlation_plots/", download.method, 
### "/ICR_Correlation_plot_",Cancer,".png"
### Manual adjustment of min and max
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path,"R tools/heatmap.3.R"))

required.packages <- c("gtools", "circlize")
ibiopak("ComplexHeatmap")
ipak(required.packages)

# Set Parameters
Pathways = "ALL"
CancerTYPES = "ALL"
Pathway_skip = ""                                                                                                        
download.method = "Assembler_Panca_Normalized_filtered"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
Log_file = paste0("./1_Log_Files/", download.method ,
                  "/5.3.Pancancer.Clustering.Oncogenic.Pathways/5.3.Pancancer.Clustering.Oncogenic.Pathways_",          # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
Legend = c("ICR Low","ICR Med","ICR High")
Gene.set = "Selected.pathways"

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/ssGSEA_", Gene.set,"_ES.Rdata"))
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers_Skipped_DLBC.Rdata"))

if (Pathways == "ALL") { 
  Pathways = rownames(ES.all)
}
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

# Create folders
dir.create("./5_Figures/",showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Oncogenic_pathway_Categorization"), showWarnings = FALSE)

dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method, "/5.3.Pancancer.Clustering.Oncogenic.Pathways"), showWarnings = FALSE)

cat("This is a log file for clustering Pancancer samples by expression of oncogenic pathways",
    "__________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Parameters Used :",
    paste0("Pathways = ", Pathways),
    paste0("Pathway_skip = ", Pathway_skip),
    paste0("download.method = ", download.method),
    "",
    "Scripts output :",
    file = Log_file,
    append = FALSE, sep= "\n")

N.pathways = length(Pathways)

Hallmark_and_ICR_cluster_assignment_allcancers = ICR_cluster_assignment_allcancers

i=24
for (i in 1:N.pathways){
  Pathway = Pathways[i]
  # Pancancer Oncogenic pathway classification
  Enrichment_score_all = t(ES.all[which(rownames(ES.all) == Pathway), , drop = FALSE])
  Enrichment_score_all = data.frame(Enrichment_score_all)
  colnames(Enrichment_score_all) = "Pathway"
  Enrichment_score_all$Pathway = as.numeric(Enrichment_score_all$Pathway)
  Enrichment_score_all$Category = NA
  Enrichment_score_all$Category[which(Enrichment_score_all$Pathway < median(Enrichment_score_all$Pathway))] = paste0(Pathway, " Low")
  Enrichment_score_all$Category[which(Enrichment_score_all$Pathway >= median(Enrichment_score_all$Pathway))] = paste0(Pathway, " High")
  
  Hallmark_and_ICR_cluster_assignment_allcancers[, paste0(Pathway, "_cluster_Pancancer")] = Enrichment_score_all$Category[match(rownames(Hallmark_and_ICR_cluster_assignment_allcancers),
                                                                                                                                  rownames(Enrichment_score_all))]
}

save(Cancer_color_table, Hallmark_and_ICR_cluster_assignment_allcancers, file = paste0("./4_Analysis/", download.method,
                                                                                       "/Pan_Cancer/Clustering/5.3.Hallmark_and_ICR_cluster_assignment_allcancers_SkippedDLBC.Rdata"))

