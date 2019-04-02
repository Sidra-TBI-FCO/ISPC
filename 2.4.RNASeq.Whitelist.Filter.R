
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

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
white_list_patients = read.csv("./2_Data/TCGA_Whitelist/Whitelist.20.03.17.csv",
                               stringsAsFactors = FALSE)
white_list_samples = read.csv("./2_Data/TCGA_Whitelist/merged_sample_quality_annotations.TCGA.PC.Synapse.v11.tsv", 
                              sep = "\t", stringsAsFactors = FALSE)


# Create folders
dir.create("./3_DataProcessing/",showWarnings = FALSE)                                                                  # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./3_DataProcessing/",download.method),showWarnings = FALSE)

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

N.sets = length(CancerTYPES)

i=1
for (i in 1:N.sets){
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(Cancer == "SKCM"){
    load(paste0("./3_DataProcessing/", download.method, "/", Cancer, "/RNASeqData/", Cancer, "_gene_RNAseq_normalized_TPandTM_filtered.Rdata"))
  }else{
    load(paste0("./3_DataProcessing/", download.method, "/", Cancer, "/RNASeqData/", Cancer, "_gene_RNAseq_normalized_TP_filtered.Rdata"))
  }
 wl_filtered.norm.RNAseqData = filtered.norm.RNAseqData[,which(substring(colnames(filtered.norm.RNAseqData), 1, 12) %in% white_list_patients$patient_barcode)]
 samples_to_include = white_list_samples$aliquot_barcode[which(white_list_samples$Do_not_use == "False")]
 wl2_filtered.norm.RNAseqData = wl_filtered.norm.RNAseqData[,which(colnames(wl_filtered.norm.RNAseqData) %in% samples_to_include)]
 rm("filtered.norm.RNAseqData")
 filtered.norm.RNAseqData = wl2_filtered.norm.RNAseqData
 save(filtered.norm.RNAseqData, geneInfo, file = paste0("./3_DataProcessing/", download.method, "/", Cancer,"/RNASeqData/", Cancer,
                                                        "_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata"))
}