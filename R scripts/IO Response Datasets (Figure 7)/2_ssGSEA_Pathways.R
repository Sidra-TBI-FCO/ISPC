## 2_ssGSEA

# GSEA Immune gene lists of expression data, heatmap

# Set-up environment
rm(list = ls())
setwd("~/Dropbox (TBI-Lab)/DB-LAB/Projects - Data/02 Active Projects/SDR400074_PIFR_2018_WH_TCGA_pancancer_immuneprofiling_in-silico/IOResponse_Datasets")
source("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/R scripts/ipak.function.R")

dir.create("./Figures", showWarnings = FALSE)

required.packages = c("GSVA")
ipak(required.packages)
required.bioconductor.packages = c("ComplexHeatmap")
ibiopak(required.bioconductor.packages)

# Define parameters
Gene_set = "Selected.pathways"
dataset = "Prat" # "Chen" "GSE78220" "GSE91061" "vanAllen" "GSE115821"

# Load data
load("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/Selected Pathways/Selected.pathways.3.3.RData")
load(file = paste0("./Data/", dataset, "/", dataset, ".Rdata"))

# Analysis
Gene.set = get(Gene_set)
all_genes_in_set = unname(unlist(Gene.set))

available_genes = intersect(all_genes_in_set, rownames(Expression.matrix))
unavailabe_genes = all_genes_in_set[-which(all_genes_in_set %in% available_genes)]

ES = gsva(Expression.matrix, Gene.set, method="ssgsea")      #Rows are genes, columns are samples

ESz = ES 
for(j in 1: nrow(ESz))  {
  ESz[j,] = (ES[j,]-mean(ES[j,]))/sd(ES[j,]) # z-score the enrichment matrix
}

Clinical.data$Response = Clinical.data$BEST.RESP

#sam.info = data.frame(Response = Clinical.data$Response, Treatment = Clinical.data$Treatment)
#sam.info = data.frame(Sample = Clinical.data$BioSample_s, Response = Clinical.data$Response)
sam.info = data.frame(Response = Clinical.data$Response)

od = hclust(dist(ESz))
#cor_mat = ESz[od, od]
nm = rownames(ESz)

col_fun = circlize::colorRamp2(c(min(ESz), 0, max(ESz)), c("blue", "white", "red"))
#col_fun = circlize::colorRamp2(c(-100, 0, 100), c("blue", "white", "red"))

ha_column = HeatmapAnnotation(df = data.frame(
                                              Response = sam.info$Response),
                              show_annotation_name = TRUE,
                              col = list(Response = c("PD" = "red",
                                                      "SD" = "orange",
                                                      "PR" = "yellow",
                                                      "CR" = "green")
                                         ),
                              show_legend = TRUE
)

# Heatmap
dir.create(paste0("./Figures/", dataset), showWarnings = FALSE)
dir.create(paste0("./Figures/", dataset, "/Heatmaps"), showWarnings = FALSE)
png(paste0("./Figures/", dataset, "/Heatmaps/",  dataset, "_heatmap_Selected_pathways.png"), res = 600, width = 13, height = 12, units = "in")
Heatmap(ESz, 
        name = "z scores", 
        cluster_rows = TRUE,
        cluster_columns = TRUE,
        row_title_gp = gpar(fontsize = 0.1),
        column_names_gp = gpar(fontsize = 10),
        row_names_gp = gpar(fontsize = 10),
        top_annotation = ha_column,
        column_title = paste0(gsub("_", " ", Gene_set), " enrichment scores in ", dataset),
        show_heatmap_legend = TRUE,
        row_names_max_width = unit(4, "in")
)

dev.off()

dir.create("./Analysis", showWarnings = FALSE)
dir.create(paste0("./Analysis/", dataset), showWarnings = FALSE)
dir.create(paste0("./Analysis/", dataset, "/Enrichment_scores"), showWarnings = FALSE)
save(ES, ESz, file = paste0("./Analysis/", dataset, "/Enrichment_scores/", 
                                      dataset ,"_Selected_Pathways_ES.Rdata"))
