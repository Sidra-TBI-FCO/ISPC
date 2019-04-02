
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"                                                                # Setwd to location were output files have to be saved.

source(paste0(code_path,"R tools/ipak.function.R"))
source(paste0(code_path,"R tools/heatmap.3.R"))

required.bioconductor.packages = c("GSVA","heatmap3", "gclus")                                                                   
ibiopak(required.bioconductor.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized_filtered"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
pw_selection_version = "3.4"
Log_file = paste0("./1_Log_Files/", download.method ,"/3.8_Deconvolution_Bindea/3.8_Deconvolution_Bindea_Log_File_",                          # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
Gene.set = "Bindea_ORIG"
Clustering = "No_Clustering"
exclude_Tgd = "exclude_Tgd"

# Load data and R scripts
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0(code_path, "Datalists/Selected.pathways.",pw_selection_version,".Rdata"))
load(paste0(code_path, "Datalists/immune.gene.lists.v3.Rdata"))

# Create folders and log file
dir.create(paste0("./5_Figures/"),showWarnings = FALSE)
dir.create(paste0("./5_Figures/Heatmaps/"),showWarnings = FALSE)
dir.create(paste0("./5_Figures/Heatmaps/", Gene.set, "_Heatmaps"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Heatmaps/", Gene.set, "_Heatmaps/", download.method), showWarnings = FALSE)
dir.create("./4_Analysis/",showWarnings = FALSE)                                                                        # Create folders to save Rdata.files
dir.create(paste0("./4_Analysis/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/", download.method), showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/", download.method, "/3.8_Deconvolution_Bindea"), showWarnings = FALSE)
cat("This is a log file for Deconvolution using Bindeas gene signatures, xCell and Hallmark pathways on RNASeq data",   # Set-up logfile
    "_________________________________________________________________________________",
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
    paste0("pathway selection version used =", pw_selection_version),
    "",
    "Scripts output :",
    "",
    "Calculating Deconvolution scores",
    file = Log_file,
    append = FALSE, sep= "\n")

# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}
N.sets = length(CancerTYPES)

start.time.process.all = Sys.time()
msg = paste0("Generating heatmaps", "\n")
cat(msg)

i = 1
for (i in 1:N.sets) {
  start.time.process.cancer = Sys.time()
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(Cancer == "LAML") 
  {cat(paste0("For ", Cancer, ", a normalization file does not exist, file is skipped.", 
              "\n",
              "-----------------------------------------------------------------------------------------------------------",
              "\n"), file = Log_file, sep = "\n", append = TRUE)
    next}
  load(paste0("./3_DataProcessing/", download.method, "/", Cancer, "/RNASeqData/", Cancer, "_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata"))
  cat (paste0 ("Generating Heatmap ",Cancer,"."))
  load(paste0("./4_Analysis/", download.method, "/", Cancer, "/Signature_Enrichment/GSEA_",
              Cancer, "_", Gene.set, ".Rdata"))
  if(Gene.set == "Bindea_ORIG" & exclude_Tgd == "exclude_Tgd"){
    ESz = ESz[-which(rownames(ESz) == "Tgd"),]
  }
  Cluster_file = paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/",
                        Cancer, "_ICR_cluster_assignment_k2-6.Rdata")
  load(Cluster_file)
  
  if(download.method == "TCGA_Assembler" | download.method == "Assembler_Panca_Normalized" | download.method == "Michele_matrix" | download.method == "Assembler_Panca_Normalized_filtered"){
    Expression.data = log(filtered.norm.RNAseqData +1, 2)
  }
  
  ## Annotation for plotting
  annotation <- table_cluster_assignment[,"HML_cluster",drop=FALSE]
  annotation$HML_cluster.col[annotation$HML_cluster=="ICR High"] = "red"
  annotation$HML_cluster.col[annotation$HML_cluster=="ICR Medium"] = "green"
  annotation$HML_cluster.col[annotation$HML_cluster=="ICR Low"] = "blue"
  
  ICR_order = c("ICR Low","ICR Medium","ICR High")
  annotation = annotation[order(match(annotation$HML_cluster,ICR_order)),]
  ColsideLabels = c("HML ICR clusters")
  Legend = c("ICR Low","ICR Med","ICR High")
  Legend_colors = c("blue","green","red")
  Colv = NULL
  annotation.blot = as.matrix(annotation[,c("HML_cluster.col"), drop= FALSE])
  ESz = ESz[,rownames(annotation.blot)]
  side.height.fraction = 0.15
  
  if(Clustering == "Clustering"){
    sHc = hclust(ddist <- dist(t(ESz)), method = "ward.D2")
    plot(sHc,labels=FALSE)
    annotation$Geneset_cluster = cutree(sHc,k = 2)[match(rownames(annotation),names(cutree(sHc,k = 2)))]
    annotation = annotation[colnames(ES),]
    annotation$Geneset_score = colMeans(ES)
    cluster_means = aggregate(Geneset_score ~ Geneset_cluster,data=annotation,FUN=mean)
    cluster_means = cluster_means[order(cluster_means$Geneset_score),]
    cluster_means$cluster_name = c(paste0(Gene.set," Low"), paste0(Gene.set, " High"))
    annotation$cluster_name = cluster_means$cluster_name[match(annotation$Geneset_cluster, cluster_means$Geneset_cluster)]
    annotation$cluster.col[annotation$cluster_name == paste0(Gene.set, " High")] = "purple"
    annotation$cluster.col[annotation$cluster_name == paste0(Gene.set, " Low")] = "pink"
    cluster_order = c(paste0(Gene.set, " Low"), paste0(Gene.set, " High"))
    annotation = annotation[, -which(colnames(annotation) %in% c("Geneset_cluster", "Geneset_score"))]
    annotation = annotation[order(match(annotation$HML_cluster,ICR_order), match(annotation$cluster_name, cluster_order)),]
    annotation.blot = as.matrix(annotation[,c("HML_cluster.col","cluster.col"), drop= FALSE])
    ColsideLabels = c("HML ICR clusters", paste0(Gene.set, " clusters"))
    Legend = c("ICR Low","ICR Med","ICR High", paste0(Gene.set, " Low"), paste0(Gene.set, " High"))
    Legend_colors = c("blue","green","red", "pink", "purple")
    Colv = as.dendrogram(sHc)
    annotation.blot = annotation.blot[colnames(ESz),]
    side.height.fraction = 0.3
  }
  
  ### Bindea plotting
  png(paste0("./5_Figures/Heatmaps/", Gene.set, "_Heatmaps/", download.method, "/Bindea_Heatmap.3_RNASeq_",Cancer, "_", Clustering, "_", exclude_Tgd, 
             ".png"),res=600,height=9,width=9,unit="in")
  heatmap.3((as.matrix(ESz)),
            main= paste0(Cancer, "\nssGSEA/", Gene.set," signatures"),
            col=my.palette,
            ColSideColors=annotation.blot,
            font_size_col_Labs = 1.5,
            cex.main = 10,
            ColSideLabs = ColsideLabels,
            Colv = Colv,
            Rowv = NA,
            keysize = 1,
            labCol=NA,
            side.height.fraction = side.height.fraction,
            cexRow = 1.3,
            margins = c(13, 11))
  
  title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, ". \n",
                          Gene.set, " enrichment z-scores were used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.65)
  legend("topright",legend = Legend,
         col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
  
  dev.off()
}
  
  
  
  
  
  