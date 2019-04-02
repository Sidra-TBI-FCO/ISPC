

rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                  # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                        # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R")) 

required.packages = c("base64enc", "HGNChelper","RCurl","httr","stringr","digest","bitops",
                      "rjson")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                 # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                     # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Michele_matrix"
assay.platform = "gene_RNAseq"
file_to_combine = "_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered" #"_gene_RNAseq_normalized_TP_filtered"   "_gene_RNAseq_normalized_TP_filtered"    

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 

# Create folders
dir.create("./3_DataProcessing/",showWarnings = FALSE)                                                                  # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./3_DataProcessing/",download.method),showWarnings = FALSE)
dir.create(paste0("./3_DataProcessing/",download.method, "/Pancancer"), showWarnings = FALSE)

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

N.sets = length(CancerTYPES)
load(paste0("./3_DataProcessing/", download.method, "/ACC/RNASeqData/ACC", file_to_combine,".Rdata"))
filtered.norm.RNAseqData.all = filtered.norm.RNAseqData

i=1
for (i in 2:N.sets){
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  #if (Cancer %in% c("LAML", "DLBC")){next}
  if(Cancer == "SKCM"){
    load(paste0("./3_DataProcessing/", download.method, "/", Cancer, "/RNASeqData/", Cancer, file_to_combine, ".Rdata"))
  }else{
    load(paste0("./3_DataProcessing/", download.method, "/", Cancer, "/RNASeqData/", Cancer, file_to_combine,".Rdata"))
  }
  filtered.norm.RNAseqData.all = cbind(filtered.norm.RNAseqData.all, filtered.norm.RNAseqData)
}
  
save(filtered.norm.RNAseqData.all, file = paste0("./3_DataProcessing/", download.method, "/Pancancer/Pancancer", file_to_combine, ".Rdata"))
  
  
  