
## Compare ICR score mutated versus wt group

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
Cancer_skip = "DLBC"                                                                                                    # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized_filtered"
Cutoff_HR = 1
ICR_k = "HML_classification"
Source_surv_data = "Cell_paper"
only_significant = "only_significant"
cancer_scaled_ICR = "not_scaled"

# Load data
load("./3_DataProcessing/External/BinaryMatrix.From.Michele.RData")
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers_Skipped_DLBC.Rdata"))
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)
Michele_genes = read.csv(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Mutational_load/Michele_genes_v2_ICRscore.csv"), stringsAsFactors = FALSE,
                         header = FALSE)
Michele_genes = c(Michele_genes[,1])

if(CancerTYPES == "ALL"){
  CancerTYPES = TCGA.cancersets$cancerType
  CancerTYPES = c("Pancancer", CancerTYPES)
}
N.sets = length(CancerTYPES)

# Add cancer scaled ICR score to ICR_cluster_assignment_all
ICR_means = aggregate(ICR_cluster_assignment_allcancers$ICRscore, by = list(ICR_cluster_assignment_allcancers$Cancer), FUN = mean)
colnames(ICR_means) = c("Cancer", "Mean_ICR")
ICR_cluster_assignment_allcancers$ICRscore_scaled = 0
ICR_cluster_assignment_allcancers$ICRscore_scaled = NA
i=2
for(i in 1:N.sets){
  Cancer = CancerTYPES[i]
  if(Cancer == "Pancancer"){next}
  Mean_ICR = ICR_means$Mean_ICR[which(ICR_means$Cancer == Cancer)]
  ICR_cluster_assignment_allcancers$ICRscore_scaled[which(ICR_cluster_assignment_allcancers$Cancer == Cancer)] = ICR_cluster_assignment_allcancers$ICRscore[which(ICR_cluster_assignment_allcancers$Cancer == Cancer)] - Mean_ICR
}

ICR_cluster_assignment_allcancers$ICRscore_scaled_antilog = 2^ICR_cluster_assignment_allcancers$ICRscore_scaled

# Some checks on binary matrix 
dim(binary_matrix) # 8315 501
length(unique(substring(rownames(binary_matrix), 1, 12))) #8315 unique patients
table(substring(rownames(binary_matrix), 14, 15)) #01 (primary tumor)= 7956;  #06 (metastatic =359) ----- Sample types in binary matrix
vial_match1 = rownames(ICR_cluster_assignment_allcancers)[which(substring(rownames(ICR_cluster_assignment_allcancers), 1, 16) %in% rownames(binary_matrix))] #8161 (of 9282) samples in ICR matrix are in binary matrix
vial_match2 = rownames(binary_matrix)[which(rownames(binary_matrix) %in% substring(rownames(ICR_cluster_assignment_allcancers), 1, 16))] #8161
patient_match1 = rownames(ICR_cluster_assignment_allcancers)[which(substring(rownames(ICR_cluster_assignment_allcancers), 1, 12) %in% substring(rownames(binary_matrix), 1, 12))] #8161 (of 9282) patients in ICR matrix are in binary matrix
patient_match2 = rownames(binary_matrix)[which(substring(rownames(binary_matrix), 1, 12) %in% substring(rownames(ICR_cluster_assignment_allcancers), 1, 12))] #8161
## patient match and vial match are both 8161

# subset binary matrix to only include Michele_genes and patients for which ICR score is available
binary_matrix = binary_matrix[which(substring(rownames(binary_matrix), 1, 16) %in% substring(rownames(ICR_cluster_assignment_allcancers), 1, 16)),which(colnames(binary_matrix) %in% Michele_genes)]
binary_matrix = binary_matrix[,order(match(colnames(binary_matrix), Michele_genes))]
dim(binary_matrix)

genes = colnames(binary_matrix)
N.genes = length(genes)

end_result = matrix(ncol = N.sets, nrow = N.genes)
colnames(end_result) = CancerTYPES
rownames(end_result) = genes

i=4
for (i in 1:N.sets){
  Cancer = CancerTYPES[i]
  if (Cancer %in% Cancer_skip) {next}
  if(Cancer == "LAML") {next}
  if(Cancer == "Pancancer"){
    matrix_sub = binary_matrix[which(substring(rownames(binary_matrix), 1, 16) %in% substring(rownames(ICR_cluster_assignment_allcancers), 1, 16)),]
    ICR_cluster_sub = ICR_cluster_assignment_allcancers[which(substring(rownames(ICR_cluster_assignment_allcancers), 1, 16) %in% rownames(matrix_sub)),]
    if(cancer_scaled_ICR == "cancer_scaled_ICR"){
      ICR_cluster_sub$ICRscore = ICR_cluster_sub$ICRscore_scaled_antilog
    }
  }else{
    ICR_cluster_sub = ICR_cluster_assignment_allcancers[which(ICR_cluster_assignment_allcancers$Cancer == Cancer),]
    matrix_sub = binary_matrix[which(substring(rownames(binary_matrix), 1, 16) %in% substring(rownames(ICR_cluster_sub), 1, 16)),]
    ICR_cluster_sub = ICR_cluster_sub[which(substring(rownames(ICR_cluster_sub), 1, 16) %in% rownames(matrix_sub)),]
    if(cancer_scaled_ICR == "cancer_scaled_ICR"){
      ICR_cluster_sub$ICRscore = ICR_cluster_sub$ICRscore_scaled_antilog
    }
  }
  
  cal_df = ICR_cluster_sub[, c("ICRscore", "HML_cluster")]
  k=69
  for (k in 1:N.genes){
    gene = genes[k]
    cal_df[, gene] = NA
    matrix_gene = matrix_sub[,gene, drop = FALSE]
    cal_df[, gene] = matrix_gene[, gene][match(substring(rownames(cal_df), 1, 16), rownames(matrix_gene))]
    if(sum(cal_df[, gene]) == 0){
      #all are unmutated
      next
    }
    #if(sum(cal_df[,gene]) == nrow(cal_df)){
      #all are mutated
      #next
    #}
    if(sum(cal_df[, gene]) < 3 | sum(cal_df[,gene]) == (nrow(cal_df) - 3)){
      next
    }
    p.value = t.test(ICRscore~get(gene), data = cal_df, paired = FALSE, var.equal = FALSE)$p.value
    means = aggregate(cal_df$ICRscore, list(cal_df[,gene]), mean)
    mean.1 = means[means$Group.1== 1, 2]
    mean.0 = means[means$Group.1== 0, 2]
    ratio = mean.1 / mean.0
    if(only_significant == "only_significant"){
      if(p.value < 0.05){
        end_result[gene, Cancer] = ratio
      }
    }else{
      end_result[gene, Cancer] = ratio
    }
  }
  
}

end_result = end_result[,-which(colnames(end_result) %in% c("LAML", "DLBC"))]
ratio_ICRscore_MUT_vs_WT = end_result

dir.create(paste0("./4_Analysis/", download.method,"/Pan_Cancer/Mutational_load/Analysis_Michele_genes/"), showWarnings = FALSE)
save(ratio_ICRscore_MUT_vs_WT, file = paste0("./4_Analysis/", download.method, "/Pan_Cancer/Mutational_load/Analysis_Michele_genes/v2_Mut_vs_wt_ICRScores_", only_significant, ".Rdata"))

# For annotation of df
load(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis/", Source_surv_data, "_Survival_analysis_High_vs_Low_Groups", 
            ICR_k, ".Rdata"))
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = ICR_disabled_cancers[-which(ICR_disabled_cancers == "LAML")]
ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])

All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype %in% c("LAML", "DLBC")),]
All_survival_analysis_data = All_survival_analysis_data[order(All_survival_analysis_data$p_value, decreasing = FALSE),]
ICR_order_part1 = as.character(All_survival_analysis_data$Cancertype[All_survival_analysis_data$HR >1])
All_survival_analysis_data = All_survival_analysis_data[order(All_survival_analysis_data$p_value, decreasing = TRUE),]
ICR_order_part2 = as.character(All_survival_analysis_data$Cancertype[All_survival_analysis_data$HR <1])
#Cancer_order = c("Pancancer", ICR_order_part1, ICR_order_part2)
Cancer_order = c("Pancancer", "SKCM", "BRCA", "LIHC", "BLCA", "HNSC",
                 "CESC", "LUSC", "LUAD", "ESCA", "GBM", "PRAD", "KIRP",
                 "KIRC", "LGG", "SARC", "THCA", "OV", "CHOL", "ACC", "KICH", "PCPG",
                 "TGCT", "THYM", "UCS", "MESO", "PAAD", "UVM", "READ", "COAD",
                 "UCEC", "STAD")

data.info = data.frame(matrix(nrow = ncol(ratio_ICRscore_MUT_vs_WT), ncol = 0), row.names = colnames(ratio_ICRscore_MUT_vs_WT), "Cancer" = NA, "ICR.Enabled.Neutral.Disabled" = NA)
data.info$Cancer = rownames(data.info)
data.info$ICR.Enabled.Neutral.Disabled[which(data.info$Cancer %in% ICR_enabled_cancers)] = "ICR enabled"
data.info$ICR.Enabled.Neutral.Disabled[which(data.info$Cancer %in% ICR_disabled_cancers)] = "ICR disabled"
data.info$ICR.Enabled.Neutral.Disabled[which(data.info$Cancer %in% ICR_neutral_cancers)] = "ICR neutral"
data.info$ICR.Enabled.Neutral.Disabled[which(data.info$Cancer == "Pancancer")] = "Pancancer"
data.info$ICR.Enabled.Neutral.Disabled = factor(data.info$ICR.Enabled.Neutral.Disabled, levels = c("Pancancer", "ICR enabled", "ICR neutral", "ICR disabled"))
data.info = data.info[order(match(data.info$Cancer, Cancer_order)),]
#col_fun = circlize::colorRamp2(c(min(ratio_ICRscore_MUT_vs_WT, na.rm = TRUE), 1, max(ratio_ICRscore_MUT_vs_WT, na.rm = TRUE)), c("blue", "white", "red"))
#col_fun = circlize::colorRamp2(c(0.5, 1, 4), c("blue", "white", "red"))
col_fun = circlize::colorRamp2(c(0.9, 1, max(ratio_ICRscore_MUT_vs_WT, na.rm = TRUE)), c("blue", "white", "red"))
ha_column = HeatmapAnnotation(df = data.frame(`ICR Enabled/Neutral/Disabled` = data.info$ICR.Enabled.Neutral.Disabled),
                              show_annotation_name = TRUE,
                              show_legend = FALSE,
                              col = list(`ICR.Enabled.Neutral.Disabled` = c("Pancancer" = "white", "ICR enabled" = "orange", 
                                                                            "ICR neutral" = "grey", "ICR disabled" = "purple"))
                                         
)

ratio_ICRscore_MUT_vs_WT = ratio_ICRscore_MUT_vs_WT[,order(match(colnames(ratio_ICRscore_MUT_vs_WT), Cancer_order))]

# Heatmap
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Mutation"), showWarnings = FALSE)
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Mutation/v2_Ratio_ICR_score_MUT_vs_WT_", only_significant, "_", cancer_scaled_ICR,".png"), res = 600, width = 14, height = 12, units = "in")
Heatmap(ratio_ICRscore_MUT_vs_WT, 
        name = "Ratio between Mean ICR score \nin Mutated versus WT samples", 
        cluster_rows = FALSE,
        cluster_columns = FALSE,
        row_title_gp = gpar(fontsize = 0.1),
        column_names_gp = gpar(fontsize = 15),
        row_names_gp = gpar(fontsize = 12),
        col = col_fun,
        top_annotation = ha_column,
        column_title = paste0("Ratio between ICR score \nMutated versus WT samples"),
        row_names_max_width = unit(9, "in")
)

dev.off()


