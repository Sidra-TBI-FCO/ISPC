
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
download.method = "Assembler_Panca_Normalized_filtered"

load(paste0(code_path, "Datalists/Selected.pathways.3.4.RData"))

Michele_genes = read.csv(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Mutational_load/Michele_genes_v2_ICRscore.csv"), stringsAsFactors = FALSE,
                         header = FALSE)
Michele_genes = c(Michele_genes[,1])

load('./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Correlation/Correlation_Oncogenic_pathways_and_ICR_irrespective_of_significance.Rdata')
pathways = rownames(pancancer_Geneset_cor_table)[order(pancancer_Geneset_cor_table[, 1], decreasing = FALSE)]

matrix = matrix(nrow = length(Selected.pathways), ncol = length(Michele_genes))
rownames(matrix) = pathways     
colnames(matrix) = Michele_genes

i = 1
for (i in 1:nrow(matrix)){
  pathway_genes = unname(unlist(Selected.pathways[i]))
  for (j in 1:ncol(matrix)){
    gene = colnames(matrix)[j]
    if(gene %in% pathway_genes){
      matrix[i,j] = "YES"
    }else{matrix[i,j] = "NO"}
  }
}
dir.create("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Mutations_in_pathway_genes", showWarnings = FALSE)
write.csv(matrix, file = "./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Mutations_in_pathway_genes/Mutations_GSEA_pathways_overlap.csv")
