#################################################################
###
### This script downloads ANY-CANCER RNASeq Data 
### from the TCGA database
### It will download and process the data. 
### Downloaded raw data is saved in:
### "./2_Data/",download.method,"/",Cancer,"/RNASeqData/"
### Processed data is saved in:
###("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData/")
###
#################################################################

##Download the RNASeq data using TCGA Assembler

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                   # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/" 

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/Module_A.R"))
source(paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/Module_B.R"))

required.packages = c("RCurl","httr", "rjson", "stringr", "HGNChelper")

ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                     # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "TCGA_Assembler"                                                                                      # Specify download method (this information to be used when saving the file)
Log_file = paste0("./1_Log_Files/1.5_RNASeq_Download/RNASeq_Download_Log_File_",                                        # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
GenomeFileHg18 = paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/SupportingFiles/Hg18GenePosition.txt")                # Specify RefGenomeFiles from Assembler_v2.0.3
GenomeFileHg19 = paste0(code_path, "R tools/TCGA-Assembler_v2.0.3/SupportingFiles/Hg19GenePosition.txt")                # Specify RefGenomeFiles from Assembler_v2.0.3

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
                                                                                                                        # in the Manual of Assembler v2.0.3 and was saved as csv file
# Create folders
dir.create(paste0("./2_Data/"),showWarnings = FALSE)                                                                    # Create folder to save downloaded raw data (by Assembler module A)
dir.create(paste0("./2_Data/",download.method), showWarnings = FALSE)

dir.create("./3_DataProcessing/",showWarnings = FALSE)                                                                  # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./3_DataProcessing/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/1.5_RNASeq_Download/"), showWarnings = FALSE)
cat("This is a log file for Downloading and Processing RNASeq data",                                                    # Set-up logfile
    "__________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Script Running Date :",
    capture.output(Sys.time()),
    "",
    "Parameters Used :",
    paste0("CancerTYPES = ", CancerTYPES),                                                          
    paste0("Cancer_skip = ", Cancer_skip),
    paste0("download.method = ", download.method),
    "",
    "Scripts output :",
    "",
    "Downloads",
    file = Log_file,
    append = FALSE, sep= "\n")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

N.sets = length(CancerTYPES)

start.time.script <- Sys.time()

#Download data
start.time.download.all <- Sys.time ()

for (i in 1:N.sets) {
  start.time.download.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  print(paste0 ("Downloading ",Cancer,"."))
  Cancer_path = paste0 ("./2_Data/",download.method,"/",Cancer,"/RNASeqData/")
  DownloadRNASeqData(cancerType = CancerTYPES[i],
                     assayPlatform = "gene_RNAseq",
                     tissueType = NULL,
                     saveFolderName = Cancer_path,
                     outputFileName = "",
                     inputPatientIDs = NULL)
  end.time.download.cancer = Sys.time()
  time = substring(as.character(capture.output(round(end.time.download.cancer - start.time.download.cancer, 2))),20,100)
  msg = paste0("Downloading time for ", Cancer, ": ",time, ".", "\n", "outputfile is downloaded to ", Cancer_path)
  cat(msg)
  cat(msg, file= Log_file,sep = "\n",append=TRUE)
}
end.time.download.all = Sys.time ()
time = substring(as.character(capture.output(round(end.time.download.all - start.time.download.all, 2))),20,100)
msg = paste0("\n", "Downloading time for all cancertypes: ",time, " min", "\n", "---------------------------------------------------------------")
cat(msg)   
cat(msg, file= Log_file,sep = "\n",append=TRUE)


## Process data
start.time.process.all = Sys.time()
msg = paste0("Processing", "\n")
cat(msg, file= Log_file,sep = "\n",append=TRUE)


for (j in 1:N.sets) {
  start.time.loop.cancer = Sys.time()
  Cancer = CancerTYPES[j]
  if (Cancer %in% Cancer_skip) {next}
  print (paste0 ("Processing ",Cancer,"."))
  Cancer_path = paste0 ("./2_Data/",download.method,"/",Cancer,"/RNASeqData/")
  file.list.all <- list.files(Cancer_path, full.names = TRUE)
  
  dir.create(paste0("./3_DataProcessing/",download.method,"/",Cancer),showWarnings = FALSE)
  dir.create(paste0("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData/"),showWarnings = FALSE)
  
  folder = paste0("./3_DataProcessing/",download.method,"/",Cancer,"/RNASeqData/")
  
  N.files = length(file.list.all)
  if (N.files>0) {
    for(k in 1:N.files){
      start.time.process.cancer = Sys.time()
      file = file.list.all[k]
      print(paste0("Processing ",file,"."))
      if(length(grep("gene_RNAseq",file))){assay.platform = "gene_RNAseq"}
      outputname = paste0(Cancer, "_", assay.platform, "_", "Processed")
      ProcessRNASeqData(inputFilePath = file,
                        outputFileName = outputname,
                        outputFileFolder = folder,
                        dataType = "geneExp",
                        verType = "RNASeqV2")
      outputfiles <- list.files(folder, full.names = TRUE)
      if(length(grep(".txt",outputfiles))){
        file.remove(outputfiles[grep(".txt",outputfiles)])
      }
      if(length(grep(".rda",outputfiles))){
        old.name <- outputfiles[grep(".rda",outputfiles)]
        new.name <- gsub(".rda", ".Rdata", old.name)
        file.rename(old.name, new.name)
      }
      end.time.process.cancer <- Sys.time ()
      time = substring(as.character(capture.output(round(end.time.process.cancer - start.time.process.cancer, 2))),20,100)
      msg = paste0("Processing time for ", file, ": ",time, ".", "\n", "Outputfile is ", outputname, ".Rdata.")
      cat(msg)
      cat(msg, file = Log_file, sep = "\n",append=TRUE)
    }
  }
}
end.time.process.all <- Sys.time ()
time <- substring(as.character(capture.output(round(end.time.process.all - start.time.process.all, 2))),20,100)
msg = paste0("\n","Processing time for all cancertypes: ",time, " min.", "\n", "---------------------------------------------------------------")
cat(msg)
cat(msg, file = Log_file,"",sep = "\n",append=TRUE)




