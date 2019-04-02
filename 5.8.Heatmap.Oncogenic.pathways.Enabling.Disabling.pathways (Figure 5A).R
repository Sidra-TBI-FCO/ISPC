
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path,"R tools/heatmap.3.R"))

required.packages <- c("gtools", "circlize", "dendsort")
ibiopak("ComplexHeatmap")
ipak(required.packages)

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = ""                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized_filtered"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"
Source_surv_data = "Cell_paper"
ICR_k = "HML_classification"
Cutoff_HR = 1
my.palette <- colorRampPalette(c("blue", "white", "red"))(n = 297)
ColsideLabels = c("ICR Cluster", "Cancers", "ICR_Enabled/Disabled")
version = "5.10."            # "5.3." (Pancancer) or "3.17." (Per cancer) or "5.10." (t.test method)

# Load data
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/ssGSEA_Selected.pathways_ES.Rdata"))
load(paste0("./4_Analysis/", download.method,"/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers_Skipped_DLBC.Rdata"))

if(version == "5.3."){
  merged = read.csv(paste0("./4_Analysis/", download.method ,"/Pan_Cancer/Survival_Analysis/", version, "All_Pathways_High_and_Low.csv"))
  pathways = gsub("_cluster_Pancancer", "", as.character(merged$Oncogenic_Pathway[-which(is.na(merged$score_contribution))]))
}
if(version == "3.17."){
  merged = read.csv(paste0("./4_Analysis/", download.method ,"/Pan_Cancer/Survival_Analysis/", version, "All_Pathways_High_and_Low.csv"))
  pathways = gsub("cluster", "", as.character(merged$Oncogenic_Pathway[-which(is.na(merged$score_contribution))]))
}
if(version == "5.10."){
  load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/T.test.ICR.enabled.disabled.cancertypes/",
              "t_tests_ES_Enabled_Disabled_Cancertypes.Rdata"))
  pathways = t.tests$pathway[-which(is.na(t.tests$contribution))]
}

ES.all = ES.all[pathways,]

load(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis/", Source_surv_data, "_Survival_analysis_High_vs_Low_Groups", 
            ICR_k, ".Rdata"))
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = ICR_disabled_cancers[-which(ICR_disabled_cancers == "LAML")]
ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])

annotation = ICR_cluster_assignment_allcancers[,c("HML_cluster", "Cancer")]
annotation$ICR_ED = NA
annotation$ICR_ED[which(annotation$Cancer %in% ICR_enabled_cancers)] = "ICR_enabled"
annotation$ICR_ED[which(annotation$Cancer %in% ICR_disabled_cancers)] = "ICR_disabled"
annotation$ICR_ED[which(annotation$Cancer %in% ICR_neutral_cancers)] = "ICR_neutral"

annotation$HML_cluster.col[annotation$HML_cluster=="ICR High"] = "red"
annotation$HML_cluster.col[annotation$HML_cluster=="ICR Medium"] = "green"
annotation$HML_cluster.col[annotation$HML_cluster=="ICR Low"] = "blue"

annotation$Cancer.col = Cancer_color_table$color[match(annotation$Cancer, Cancer_color_table$Group.1)]

annotation$ICR_ED.col[annotation$ICR_ED == "ICR_enabled"] = "orange"
annotation$ICR_ED.col[annotation$ICR_ED == "ICR_disabled"] = "purple"
annotation$ICR_ED.col[annotation$ICR_ED == "ICR_neutral"] = "grey"

All_survival_analysis_data = All_survival_analysis_data[order(All_survival_analysis_data$HR, decreasing = TRUE),]
Cancer_order = as.character(All_survival_analysis_data$Cancertype[-which(All_survival_analysis_data$Cancertype == "LAML")])
ICR_order = c("ICR Low","ICR Medium","ICR High")
ICR_ED_order = c("ICR_enabled", "ICR_disabled")
##

annotation = annotation[order(match(annotation$ICR_ED, ICR_ED_order), match(annotation$Cancer, Cancer_order)),]
annotation.blot = as.matrix(annotation[,c("HML_cluster.col", "Cancer.col", "ICR_ED.col"), drop = FALSE])
annotation.blot = annotation.blot[which(rownames(annotation.blot) %in% colnames(ES.all)),]
#annotation.blot = annotation.blot[colnames(Expression.data),]                                                                                        # The sample order in annotation.blot needs to be the same as in Expression.data
#Expression.data = Expression.data[colnames(annotation.blot),]

ESz.all = ES.all 
j=1
for(j in 1: nrow(ESz.all))  {
  ESz.all[j,] = (ES.all[j,]-mean(ES.all[j,]))/sd(ES.all[j,]) # z-score the enrichment matrix
}
ESz.all = ESz.all[,rownames(annotation.blot)]

ICR_cluster_assignment_allcancers$ICR_ED = NA
ICR_cluster_assignment_allcancers$ICR_ED[ICR_cluster_assignment_allcancers$Cancer %in% ICR_enabled_cancers] = "ICR_enabled"
ICR_cluster_assignment_allcancers$ICR_ED[ICR_cluster_assignment_allcancers$Cancer %in% ICR_disabled_cancers] = "ICR_disabled"

## Determine order rows/signatures
to_aggregate = t(ESz.all)
to_aggregate = as.data.frame(to_aggregate)
to_aggregate$ICR_ED = ICR_cluster_assignment_allcancers$ICR_ED[match(rownames(to_aggregate), rownames(ICR_cluster_assignment_allcancers))]
mean_ES_all = aggregate(.~ ICR_ED, data = to_aggregate, FUN = mean)
mean_ES_all = t(mean_ES_all)
colnames(mean_ES_all) = mean_ES_all[1,]
mean_ES_all = mean_ES_all[-1,]
mode(mean_ES_all) = "numeric"
mean_ES_all = as.data.frame(mean_ES_all)
mean_ES_all$DeltaED = c(mean_ES_all$ICR_enabled - mean_ES_all$ICR_disabled)
mean_ES_all = mean_ES_all[order(mean_ES_all$DeltaED, decreasing = TRUE),]
cell_order = rownames(mean_ES_all)
ESz.all = ESz.all[c(cell_order),]

save(cell_order, file = paste0("./5_Figures/Pancancer_plots/", download.method, "/Heatmaps_Enabled_Disabled/", version, "_pathway_order.Rdata"))

## Cluster tree
sHc = hclust(ddist <- dist(t(ESz.all)), method = "ward.D2")
plot(sHc,labels=FALSE)

sHc = as.dendrogram(sHc)

#sHc = dendsort(sHc, type = "min", isReverse = FALSE)
#sHc.dd = as.dendrogram(sHc)

x = as.data.frame(t(ESz.all))
wts_c = colSums(x) 
test = reorder(sHc, wts = wts_c)
#test = reorder(sHc.dd, wts = wts_c)

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/5.3.Hallmark_and_ICR_cluster_assignment_allcancers_SkippedDLBC.Rdata"))
#x$cluster = cutree(sHc, k=2)[match(rownames(x), names(cutree(sHc, k=2)))]
#x$cluster[which(x$cluster == 1)] = "beneficial"
#x$cluster[which(x$cluster == 2)] = "non-beneficial"

#Hallmark_and_ICR_cluster_assignment_allcancers$ICR_benefit_cluster = x$cluster[match(rownames(Hallmark_and_ICR_cluster_assignment_allcancers),
#                                                                                     rownames(x))]

#save(Hallmark_and_ICR_cluster_assignment_allcancers, file = paste0("./4_Analysis/", download.method,"/Pan_Cancer/Clustering/", version, "and.5.8.Hallmark_and_ICR_cluster_assignment_allcancers_ICRbenefit_cluster.Rdata"))

## Plot prep
Legend = Cancer_order[-which(Cancer_order %in% c("DLBC", "LAML"))]
rownames(Cancer_color_table) = Cancer_color_table$Group.1
Cancer_color_table = Cancer_color_table[Cancer_order,]
Legend_colors = c(Cancer_color_table$color)

#Legend2 = c("ICR Low","ICR Med","ICR High")
#Legend_colors2 = c("blue","green","red")

Legend3 = c("ICR enabled", "ICR neutral","ICR disabled")
Legend_colors3 = c("orange", "grey", "purple")

rownames(ESz.all) = gsub(".*] ", "",rownames(ESz.all))

dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Heatmaps_Enabled_Disabled"), showWarnings = FALSE)
### Plotting
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Heatmaps_Enabled_Disabled/Final_", version, "and.5.8.Row_clusters_Heatmap_Enabled_Disabled_RNASeq_Pancancer_HML_final_figure.png"), 
    res = 600, height = 15, width = 15, unit = "in")
heatmap.3((as.matrix(ESz.all)),
          main= paste0("Pancancer enrichment scores \nssGSEA/Oncogenic pathway signatures"),
          col=my.palette,
          ColSideColors=annotation.blot,
          font_size_col_Labs = 1.5,
          cex.main = 10,
          ColSideLabs = ColsideLabels,
          #Colv= as.dendrogram(sHc),
          Colv = test,
          dendrogram = "column",
          #Rowv = NULL,
          labCol=NA,
          side.height.fraction = 0.25,
          cexRow = 1.3,
          margins = c(13, 30))

title(sub = list(paste0("Figure: EDAseq normalized, log transformed gene expression data \n was obtained from TCGA, using ", download.method, ". \n",
                        "Oncogenic pathway enrichment z-scores were used to generate this heatmap."), cex = 1), outer = FALSE, line = -1, adj = 0.65)
#legend("bottomleft",legend = Legend,
       #col = Legend_colors,lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
#legend("top", legend = Legend2,
#col = Legend_colors2, lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)
legend("topright", legend = Legend3,
       col = Legend_colors3, lty= 0.5,lwd = 0.5, cex = 0.7, pch= 15, pt.cex = 1.0)

dev.off()

## Color Key

png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Heatmaps_Enabled_Disabled/", version, "and.5.8.Color_Key.png"), 
    res = 600, height = 7, width = 7, unit = "in")
heatmap.3((as.matrix(ESz.all)),
          main= paste0("Pancancer enrichment scores \nssGSEA/Oncogenic pathway signatures"),
          col=my.palette,
          ColSideColors=annotation.blot,
          font_size_col_Labs = 1.5,
          cex.main = 10,
          ColSideLabs = ColsideLabels,
          #Colv= as.dendrogram(sHc),
          #Colv = NULL,
          Rowv = NULL,
          labCol=NA,
          side.height.fraction = 0.25,
          cexRow = 1.3,
          margins = c(13, 30))
dev.off()



