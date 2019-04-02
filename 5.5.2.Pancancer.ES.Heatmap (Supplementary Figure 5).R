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

required.packages <- c("gtools", "circlize")
ibiopak("ComplexHeatmap")
ipak(required.packages)                                                                                                  # Install and load required packages

# Set Parameters
CancerTYPES = c("ALL")                                                                                                   # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "Assembler_Panca_Normalized_filtered"                                                                       # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"
Gene.set = "Bindea_ORIG"                                                                                           # c("Bindea", "bindea_patrick", "xCell")
Plot_type = "Mean_all"                                                                                                  # "Mean_all" or "ICR_High_vs_Low"     
enrichment.score.type = ".enrichment.score"
scaling = "z_score"
Source_surv_data = "Cell_paper"
Cutoff_HR = 1
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
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/DotHeatmaps"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method,"/DotHeatmaps/", Gene.set, "_enrichment_pancancer_Dot_Heatmap"), showWarnings = FALSE)

N.sets = length(CancerTYPES)

start.time = Sys.time ()

load(paste0("./4_Analysis/", download.method, "/ACC/Signature_Enrichment/GSEA_ACC_", Gene.set,".Rdata"))

Cancers_included = CancerTYPES[-which(CancerTYPES %in% c("LAML", "DLBC"))]
Deconvolution_score_per_ICR_df = data.frame(matrix(nrow = length(rownames(ES)), ncol = length(Cancers_included) * 4), row.names = rownames(ES))
ICR_Highs = paste(Cancers_included, "_ICR_High", sep = "")
ICR_Lows = paste(Cancers_included, "_ICR_Low", sep = "")
ICR_delta = paste(Cancers_included, "_ICR_High_vs_Low", sep = "")
Mean_all = paste(Cancers_included, "_Mean_all", sep = "")
colnames(Deconvolution_score_per_ICR_df) = c(ICR_Highs, ICR_Lows, ICR_delta, Mean_all)

i=1
for (i in 1:N.sets){
  Cancer = CancerTYPES[i]
  if(Cancer %in% c("LAML", "DLBC")){next}
  load(paste0("./4_Analysis/", download.method, "/", Cancer, "/Signature_Enrichment/GSEA_", Cancer, "_", Gene.set,".Rdata"))
  load(paste0("./4_Analysis/", download.method, "/", Cancer, "/Clustering/", Cancer, ".", download.method, ".EDASeq.ICR.reps5000/", 
              Cancer, "_ICR_cluster_assignment_k2-6.Rdata"))
  
  ICR_high_patients = rownames(table_cluster_assignment)[which(table_cluster_assignment$HML_cluster == "ICR High")]
  ICR_low_patients = rownames(table_cluster_assignment)[which(table_cluster_assignment$HML_cluster == "ICR Low")]
  ES = as.data.frame(ES)
  ES$Mean_all = rowMeans(ES)
  ES$ICR_High_mean = rowMeans(ES[, which(colnames(ES) %in% ICR_high_patients)])
  ES$ICR_Low_mean = rowMeans(ES[, which(colnames(ES) %in% ICR_low_patients)])
  
  ICR_High_Cancer = paste0(Cancer, "_ICR_High")
  Deconvolution_score_per_ICR_df[,which(colnames(Deconvolution_score_per_ICR_df) == ICR_High_Cancer)] = ES$ICR_High_mean[match(rownames(ES),
                                                                                                                               rownames(Deconvolution_score_per_ICR_df))]
                                                                                                                                                                
  ICR_Low_Cancer = paste0(Cancer, "_ICR_Low")
  Deconvolution_score_per_ICR_df[,which(colnames(Deconvolution_score_per_ICR_df) == ICR_Low_Cancer)] = ES$ICR_Low_mean[match(rownames(ES),
                                                                                                                             rownames(Deconvolution_score_per_ICR_df))]
                                                                                                                                                              
  ICR_Delta_Cancer = paste0(Cancer, "_ICR_High_vs_Low")
  Deconvolution_score_per_ICR_df[,which(colnames(Deconvolution_score_per_ICR_df) == ICR_Delta_Cancer)] = Deconvolution_score_per_ICR_df[,which(colnames(Deconvolution_score_per_ICR_df) == ICR_High_Cancer)] -
    Deconvolution_score_per_ICR_df[,which(colnames(Deconvolution_score_per_ICR_df) == ICR_Low_Cancer)]
  
  Mean_all_Cancer = paste0(Cancer, "_Mean_all")
  Deconvolution_score_per_ICR_df[,which(colnames(Deconvolution_score_per_ICR_df) == Mean_all_Cancer)] = ES$Mean_all[match(rownames(ES),
                                                                                                                    rownames(Deconvolution_score_per_ICR_df))]
                                                                                                                                                           
}

if(Gene.set == "Bindea_ORIG" & exclude_Tgd == "exclude_Tgd"){
  Deconvolution_score_per_ICR_df = Deconvolution_score_per_ICR_df[-which(rownames(Deconvolution_score_per_ICR_df) == "Tgd"),]
}

if(Plot_type == "ICR_High_vs_Low"){
  Deconvolution_score_for_plot = Deconvolution_score_per_ICR_df[, grep(pattern = "_ICR_High_vs_Low", colnames(Deconvolution_score_per_ICR_df))]
  colnames(Deconvolution_score_for_plot) = gsub("_ICR_High_vs_Low","", colnames(Deconvolution_score_for_plot))
}

if(Plot_type == "Mean_all"){
  Deconvolution_score_for_plot = Deconvolution_score_per_ICR_df[, grep(pattern = "_Mean_all", colnames(Deconvolution_score_per_ICR_df))]
  colnames(Deconvolution_score_for_plot) = gsub("_Mean_all","", colnames(Deconvolution_score_for_plot))
}

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/", Source_surv_data ,"_Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])
All_survival_analysis_data = All_survival_analysis_data[order(All_survival_analysis_data$HR, decreasing = TRUE),]
Enabled = All_survival_analysis_data[which(All_survival_analysis_data$Cancer %in% ICR_enabled_cancers),]
Disabled = All_survival_analysis_data[which(All_survival_analysis_data$Cancer %in% ICR_disabled_cancers),]
Neutral = All_survival_analysis_data[which(All_survival_analysis_data$Cancer %in% ICR_neutral_cancers),]
Cancer_order = c(as.character(Enabled$Cancertype), as.character(Neutral$Cancertype), as.character(Disabled$Cancertype))
Deconvolution_score_for_plot = Deconvolution_score_for_plot[,order(match(colnames(Deconvolution_score_for_plot), Cancer_order))]

#load("./4_Analysis/Assembler_Panca_Normalized/Pan_Cancer/Signature_Enrichment/cell_order_Trajanoski_heatmap.Rdata")

load(paste0("./5_Figures/Pancancer_plots/", download.method, "/ES_Heatmaps/", Gene.set, "_Cell.order.Rdata"))
Deconvolution_score_for_plot = Deconvolution_score_for_plot[order(match(rownames(Deconvolution_score_for_plot), cell_order)),]

#Deconvolution_score_for_plot = Deconvolution_score_for_plot[order(rowMeans(Deconvolution_score_for_plot), decreasing = TRUE),]

Deconvolution_score_for_plot = data.matrix(Deconvolution_score_for_plot)

Deconvolution_score_for_plot_z_scored = Deconvolution_score_for_plot
for(j in 1: nrow(Deconvolution_score_for_plot))  {
  Deconvolution_score_for_plot_z_scored[j,] = (Deconvolution_score_for_plot[j,]-mean(Deconvolution_score_for_plot[j,]))/sd(Deconvolution_score_for_plot[j,]) # z-score the enrichment matrix
}

if(scaling == "z_score"){
  Deconvolution_score_for_plot = Deconvolution_score_for_plot_z_scored
}

min = -4.4
max = 4.4

col_fun = circlize::colorRamp2(c(min, 0, max), c("blue", "white", "red"))

#DOT HEATMAP
#png(paste0("./5_Figures/Pancancer_plots/DotHeatmaps/", download.method, "/", Bindea, "_DotHeatmap_ICR_High_vs_Low.png"), res=600,height=10,width=10,unit="in")
png(paste0("./5_Figures/Pancancer_plots/", download.method,"/DotHeatmaps/", Gene.set, "_enrichment_pancancer_Dot_Heatmap/5.5.2.v2.",
           Gene.set, "_DotHeatmap_", Plot_type, "_", exclude_Tgd,".png"), res = 600, height = 12, width = 20, units = "in")

heatmap = Heatmap(Deconvolution_score_for_plot,
                  show_heatmap_legend = TRUE,
                  cluster_rows = FALSE,
                  cluster_columns = FALSE,
                  rect_gp = gpar(type = "none"),
                  column_names_side = "bottom",
                  column_names_gp = gpar(fontsize =20),
                  row_names_gp = gpar(fontsize = 20),
                  col = colorRamp2(c(min, 0, max), c("blue", "white", "red")),
                  heatmap_legend_param = list(at= c(min, 0, max), labels = c(min, "0", max),                             ## Make sure this is the same as col_fun!!
                                              legend_direction = "horizontal", legend_height = unit(5, "in"), title_position = "lefttop"),
                  #top_annotation = ha_column,
                  name = "Mean enrichment score (z-scored by row)",
                  row_names_max_width = unit(7, "in"),
                  #gap = unit(5, "mm"),
                  row_title_gp = gpar(fontsize = 20),
                  cell_fun = function(j, i, x, y, width, height, fill) {
                    grid.rect(x = x, y = y, width = width, height = height, gp = gpar(col = "grey", fill = NA))
                    grid.circle(x = x, y = y, r = 0.015,gp = gpar(fill = col_fun(Deconvolution_score_for_plot[i, j]), col = NA))
                  }
)

print(heatmap)

#title(main = paste0(Bindea, " signature expression in ICR High tumors compared to ICR Low tumors"))
dev.off()

