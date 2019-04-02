#################################################################
###
###
### Data input:
### "./4_Analysis/",download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, 
### "_Bindea_xCell_HallmarkPathways.Rdata")
### Output :
### "./5_Figures/Correlation_plots/ICR_Correlation_plots/", download.method, 
### "/ICR_Correlation_plot_",Cancer,".png"
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
source(paste0(code_path,"R tools/heatmap.3.R"))

required.packages <- c("gtools", "circlize")
ibiopak("ComplexHeatmap")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized_filtered"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
Log_file = paste0("./1_Log_Files/", download.method ,"/5.5.1.Pancancer_Bindea_Heatmap/5.5.1.Pancancer.Bindea.Heatmap_", # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
ColsideLabels = c("ICR clusters", "Cancers", "ICR Enabled/Neutral/Disabled")
Legend = c("ICR Low","ICR Med","ICR High")
Gene.set = "Bindea_ORIG"                                                                 # "Bindea.enrichment.score", "bindea_patrick.enrichment.score"
Cutoff_HR = 1
Source_surv_data = "Cell_paper"
exclude_Tgd = "exclude_Tgd"

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
# in the Manual of Assembler v2.0.3 and was saved as csv file.
# Define parameters (based on loaded data)
if (CancerTYPES == "ALL") { 
  CancerTYPES = TCGA.cancersets$cancerType
}

# Create folders
dir.create("./5_Figures/",showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method,"/ES_Heatmaps"), showWarnings = FALSE)

load(paste0("./4_Analysis/", download.method, "/ACC/Signature_Enrichment/GSEA_ACC", "_", Gene.set,".Rdata"))
ES.all = ES

N.sets = length(CancerTYPES)

for (i in 2:N.sets){
  Cancer = CancerTYPES[i]
  if(Cancer %in% c("LAML", "DLBC")){next}
  rm(ES)
  load(paste0("./4_Analysis/", download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, "_", Gene.set,".Rdata"))
  ES.all = cbind(ES.all, ES)
}

if(Gene.set == "Bindea_ORIG" & exclude_Tgd == "exclude_Tgd"){
  ES.all = ES.all[-which(rownames(ES.all) == "Tgd"),]
}

dir.create(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment"), showWarnings = FALSE)
save(ES.all, file = paste0("./4_Analysis/", download.method,"/Pan_Cancer/Signature_Enrichment/ssGSEA_", Gene.set, "_ES.Rdata"))

## t-test
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers_Skipped_DLBC.Rdata"))
ES.test = as.data.frame(t(ES.all))

t.tests <- data.frame (cell.types= colnames(ES.test),p.value=0,p.value.fdr=0,mean.ICR.Low=0,mean.ICR.High=0,stringsAsFactors = FALSE)
t.tests[t.tests == 0] = NA

N.cells = length(colnames(ES.test))
cell.types = colnames(ES.test)

i=21
for (i in 1:N.cells){
  cell = cell.types[i]
  df = ES.test[,c(cell), drop=FALSE]
  df$ICR = ICR_cluster_assignment_allcancers$HML_cluster[match(rownames(df), rownames(ICR_cluster_assignment_allcancers))]
  colnames(df) = c("cell_type", "ICR")
  subset = df[which(df$ICR %in% c("ICR High", "ICR Low")),]
  test = t.test(cell_type~ICR, data = subset, paired = FALSE, var.equal = FALSE)
  mean_ICR_Low = unname(test$estimate["mean in group ICR Low"])
  mean_ICR_High = unname(test$estimate["mean in group ICR High"])
  p.value = test$p.value
  p.value.fdr = signif(p.adjust(p = p.value, method = "fdr", n = ncol(ES.test)), 4)
  t.tests[t.tests$cell.types==cell,c(2:ncol(t.tests))] <- c(p.value,p.value.fdr, mean_ICR_Low, mean_ICR_High)
}

t.tests$sig = NA
t.tests$sig[which(t.tests$p.value.fdr < 0.05)] = "*"
#t.tests$inverse_log = -log(t.tests$p.value)

## Plotting
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers_Skipped_DLBC.Rdata"))

annotation = ICR_cluster_assignment_allcancers[,c("HML_cluster", "Cancer")]

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/", Source_surv_data ,"_Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype %in% c("LAML", "DLBC")),]
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])

annotation$ICR_ED = NA
annotation$ICR_ED[which(annotation$Cancer %in% ICR_enabled_cancers)] = "ICR_enabled"
annotation$ICR_ED[which(annotation$Cancer %in% ICR_neutral_cancers)] = "ICR_neutral"
annotation$ICR_ED[which(annotation$Cancer %in% ICR_disabled_cancers)] = "ICR_disabled"
annotation$ICR_ED.col[annotation$ICR_ED == "ICR_enabled"] = "orange"
annotation$ICR_ED.col[annotation$ICR_ED == "ICR_neutral"] = "grey"
annotation$ICR_ED.col[annotation$ICR_ED == "ICR_disabled"] = "purple"

annotation$HML_cluster.col[annotation$HML_cluster=="ICR High"] = "red"
annotation$HML_cluster.col[annotation$HML_cluster=="ICR Medium"] = "green"
annotation$HML_cluster.col[annotation$HML_cluster=="ICR Low"] = "blue"

annotation$Cancer.col = Cancer_color_table$color[match(annotation$Cancer, Cancer_color_table$Group.1)]

ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])
All_survival_analysis_data = All_survival_analysis_data[order(All_survival_analysis_data$HR, decreasing = TRUE),]
Enabled = All_survival_analysis_data[which(All_survival_analysis_data$Cancer %in% ICR_enabled_cancers),]
Disabled = All_survival_analysis_data[which(All_survival_analysis_data$Cancer %in% ICR_disabled_cancers),]
Neutral = All_survival_analysis_data[which(All_survival_analysis_data$Cancer %in% ICR_neutral_cancers),]
Cancer_order = c(as.character(Enabled$Cancertype), as.character(Neutral$Cancertype), as.character(Disabled$Cancertype))

ICR_order = c("ICR Low","ICR Medium","ICR High")

annotation = annotation[order(match(annotation$HML_cluster,ICR_order), match(annotation$Cancer, Cancer_order)),]

annotation.blot = as.matrix(annotation[,c("HML_cluster.col","Cancer.col", "ICR_ED.col"), drop = FALSE])
#annotation.blot = annotation.blot[colnames(Expression.data),]                                                                                        # The sample order in annotation.blot needs to be the same as in Expression.data
#Expression.data = Expression.data[colnames(annotation.blot),]

ESz.all = ES.all 
for(j in 1: nrow(ESz.all))  {
  ESz.all[j,] = (ES.all[j,]-mean(ES.all[j,]))/sd(ES.all[j,]) # z-score the enrichment matrix
}
ESz.all = ESz.all[,rownames(annotation.blot)]

## Determine order rows/signatures
to_aggregate = t(ESz.all)
to_aggregate = as.data.frame(to_aggregate)
to_aggregate$ICR_cluster = ICR_cluster_assignment_allcancers$HML_cluster[match(rownames(to_aggregate), rownames(ICR_cluster_assignment_allcancers))]
mean_ES_all = aggregate(.~ ICR_cluster, data = to_aggregate, FUN = mean)
mean_ES_all = t(mean_ES_all)
colnames(mean_ES_all) = mean_ES_all[1,]
mean_ES_all = mean_ES_all[-1,]
mode(mean_ES_all) = "numeric"
mean_ES_all = as.data.frame(mean_ES_all)
mean_ES_all$DeltaHL = c(mean_ES_all$`ICR High` - mean_ES_all$`ICR Low`)
mean_ES_all = mean_ES_all[order(mean_ES_all$DeltaHL, decreasing = TRUE),]
cell_order = rownames(mean_ES_all)

save(cell_order, file = paste0("./5_Figures/Pancancer_plots/", download.method, "/ES_Heatmaps/", Gene.set, "_Cell.order.Rdata"))
write.csv(mean_ES_all, file = paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/mean_ES_all_", Gene.set,".csv"))

## order t-tests by cell_order
t.tests = t.tests[match(cell_order, t.tests$cell.types),]

## Plot prep
ESz.all = ESz.all[c(cell_order),]
Legend = Cancer_order
Cancer_color_table = Cancer_color_table[-which(Cancer_color_table$Group.1 == "DLBC"),]
Legend_colors = c(Cancer_color_table$color[order(match(Cancer_color_table$Group.1, Cancer_order))])

Legend2 = c("ICR Low","ICR Med","ICR High")
Legend_colors2 = c("blue","green","red")

Legend3 = c("ICR enabled", "ICR neutral", "ICR disabled")
Legend_colors3 = c("orange", "grey", "purple")

### Plotting

png(paste0("./5_Figures/Pancancer_plots/", download.method, "/ES_Heatmaps/5.5.1.", Gene.set, "_", exclude_Tgd, "_Enabled_Disabled_Neutral_Cancers_reordered_Heatmap_RNASeq_Pancancer_HML.png"), res = 600, height = 10, width = 15, unit = "in")
heatmap.3((as.matrix(ESz.all)),
          main= paste0("Pancancer enrichment scores \nssGSEA/ ", Gene.set, " signature"),
          col= my.palette,
          ColSideColors=annotation.blot,
          font_size_col_Labs = 1.5,
          cex.main = 10,
          ColSideLabs = ColsideLabels,
          #Colv= as.dendrogram(sHc),
          Colv = NULL,
          Rowv = NULL,
          labCol=NA,
          side.height.fraction = 0.3,
          cexRow = 1.3,
          margins = c(13, 30))

title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, ". \n",
                        "Immune geneset enrichment z-scores were used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.60)
legend("bottomleft",legend = Legend,
       col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
legend("topright", legend = Legend2,
       col = Legend_colors2, lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
legend("top", legend = Legend3,
       col = Legend_colors3, lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
dev.off()

# Color key
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/ES_Heatmaps/", Gene.set,"_Color_Key.png"), res = 600, height = 7, width = 7, unit = "in")
heatmap.3((as.matrix(ESz.all)),
          main= paste0("Pancancer enrichment scores \nssGSEA/ ", Gene.set, " signature"),
          col=my.palette,
          ColSideColors=annotation.blot,
          font_size_col_Labs = 1.5,
          cex.main = 10,
          ColSideLabs = ColsideLabels,
          #Colv= as.dendrogram(sHc),
          Colv = NULL,
          Rowv = NULL,
          labCol=NA,
          side.height.fraction = 0.3,
          cexRow = 1.3,
          margins = c(13, 30))
dev.off()


