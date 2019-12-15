#################################################################
###
###
### Data input:
###
### Output :
### 
#################################################################


## Remark: The forest plot function has issues with zero and/or inf number values.
## Beware of this for troubleshooting!

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path, "R tools/heatmap.3.R"))

required.packages <- c("RColorBrewer", "forestplot")
ipak(required.packages)

## Set Parameters
CancerTYPES = c("ALL")                                                                                                   # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "Assembler_Panca_Normalized_filtered"                                                                                       # Specify download method (this information to be used when saving the file)
my.palette = colorRampPalette(c("#152B7E", "white", "#1B7E09"))(n = 297)                                                 # Specify which genes will be correlated
Log_file = paste0("./1_Log_Files/", download.method, "/5.1_Pancancer_Correlation_matrix_Signatures/5.1_Pancancer_Correlation_matrix_signatures", 
                  "_Bindea_xCell_Hallmark", "_Log_File_", gsub(":",".",gsub(" ","_",date())),".txt")
assay.platform = "gene_RNAseq"
Source_surv_data = "Cell_paper"
Outcome = "OS"

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations)

if(Outcome == ""){
  load(paste0("./4_Analysis/", download.method,"/Pan_Cancer/Survival_Analysis/", Source_surv_data, "_Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))
}else{
  load(paste0("./4_Analysis/", download.method,"/Pan_Cancer/Survival_Analysis/",Outcome, "_", Source_surv_data, "_Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))
}
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers_Skipped_DLBC.Rdata"))

i =1
# Add number of patients per cancer to All_survival_analysis data
for (i in 1:length(TCGA.cancersets$cancerType)){
  Cancer = TCGA.cancersets$cancerType[i]
  if(Source_surv_data == "TCGA_Assembler"){
    Survival_data = read.csv(paste0("./3_DataProcessing/TCGA_Assembler/",Cancer,"/SurvivalData/updatedsurvivaldata.csv"))
  }
  if(Source_surv_data == "Cell_paper"){
    Survival_data = read.csv("./2_Data/TCGA cell 2018 clinical/TCGA_CLINICAL_DATA_CELL_2018_S1.csv",
                             stringsAsFactors = FALSE)
    Survival_data = Survival_data[which(Survival_data$type == Cancer),]
  }
  Survival_data = Survival_data[which(Survival_data$bcr_patient_barcode %in% 
                                  substring(rownames(ICR_cluster_assignment_allcancers), 1, 12)),]
  
  All_survival_analysis_data$N[All_survival_analysis_data$Cancertype == Cancer] = nrow(Survival_data)
  
  table_cluster_assignment = ICR_cluster_assignment_allcancers[which(ICR_cluster_assignment_allcancers$Cancer == Cancer),]
  nICR_High = nrow(table_cluster_assignment[which(table_cluster_assignment$HML_cluster == "ICR High"),])
  nICR_Low = nrow(table_cluster_assignment[which(table_cluster_assignment$HML_cluster == "ICR Low"),])
  All_survival_analysis_data$ICR_High[All_survival_analysis_data$Cancertype == Cancer] = nICR_High
  All_survival_analysis_data$ICR_Low[All_survival_analysis_data$Cancertype == Cancer] = nICR_Low
}

if (CancerTYPES == "ALL") { 
  CancerTYPES = TCGA.cancersets$cancerType
}

# Remove LAML + DLBC 
# and TGCT + PCPG (have zero and inf values)
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype == "DLBC"),]
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype == "LAML"),]
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype == "TGCT"),]
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype == "THYM"),]
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype == "PCPG"),]

#All_survival_analysis_data[which(All_survival_analysis_data$Cancertype == "TGCT"), "CI_lower"] = 0.001
#All_survival_analysis_data[which(All_survival_analysis_data$Cancertype == "TGCT"), "CI_upper"] = 16
#All_survival_analysis_data[which(All_survival_analysis_data$Cancertype == "TGCT"), "HR"] = 1
#All_survival_analysis_data[which(All_survival_analysis_data$Cancertype == "THYM"), "CI_lower"] = 0.001
#All_survival_analysis_data[which(All_survival_analysis_data$Cancertype == "THYM"), "CI_upper"] = 16
#All_survival_analysis_data[which(All_survival_analysis_data$Cancertype == "THYM"), "HR"] = 1
#All_survival_analysis_data[which(All_survival_analysis_data$Cancertype == "PCPG"), "CI_lower"] = 0.001
#All_survival_analysis_data[which(All_survival_analysis_data$Cancertype == "PCPG"), "CI_upper"] = 16
#All_survival_analysis_data[which(All_survival_analysis_data$Cancertype == "PCPG"), "HR"] = 1

All_survival_analysis_data$p_value = signif(All_survival_analysis_data$p_value, digits = 3)
HR.table = All_survival_analysis_data

#HR.table <- HR.table[-which(is.infinite(HR.table$Upper)),]
#HR.table <- HR.table[-which(HR.table$Lower==0),]

n.cancers = length(All_survival_analysis_data$Cancertype)
x = n.cancers + 2

HR.table = HR.table[order(HR.table$HR),]

# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta =
  structure(list(
    mean  = as.numeric(c(NA,HR.table$HR[1:n.cancers]), NA),
    lower = c(NA,HR.table$CI_lower[c(1:n.cancers)], NA),
    upper = c(NA,HR.table$CI_upper[c(1:n.cancers)], NA)),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -x),
    class = "data.frame")

tabletext<-cbind(
  c("Cancer", as.character(HR.table$Cancertype)[c(1:n.cancers)]),
  c("N", HR.table$N[c(1:n.cancers)]),
  c("N events", HR.table$N_events[c(1:n.cancers)]),
  c("N ICR High", HR.table$ICR_High[c(1:n.cancers)]),
  c("N ICR Low", HR.table$ICR_Low[c(1:n.cancers)]),
  c("p-value", HR.table$p_value[c(1:n.cancers)]),
  c("HR",      HR.table$HR[c(1:n.cancers)]))

#dev.new()
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Forest_plots"), showWarnings = FALSE)

if(Outcome == ""){
  pdf(file = paste0("./5_Figures/Pancancer_plots/", download.method, "/Forest_plots/3.16.", Source_surv_data, "_forest_plot_Skipped_DBLC_LAML_TGCT_THYM_PCPG.pdf"),
      height = 6.5, width = 13)
}else{
  pdf(file = paste0("./5_Figures/Pancancer_plots/", download.method, "/Forest_plots/3.16.Final_Revised_", Outcome, "_", Source_surv_data, "_forest_plot_Skipped_DBLC_LAML_TGCT_THYM_PCPG.pdf"),
      height = 6.5, width = 13)
}

forestplot(tabletext,
           cochrane_from_rmeta,new_page = TRUE,
           is.summary=c(TRUE,rep(FALSE,n.cancers),TRUE,rep(FALSE,n.cancers),TRUE,FALSE),
           #clip=c(0,8),
           xlog=TRUE,
           #xlim = c(0.010, 15),
           boxsize = .25,
           vertices = TRUE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 13), xlab = gpar(fontsize = 20), cex = 1))
dev.off()
