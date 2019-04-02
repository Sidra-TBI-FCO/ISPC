
### Prepare Rdata files from RDS files

# Set-up environment
rm(list = ls())
setwd("~/Dropbox (TBI-Lab)/DB-LAB/Projects - Data/02 Active Projects/SDR400074_PIFR_2018_WH_TCGA_pancancer_immuneprofiling_in-silico/IOResponse_Datasets")
source("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/R scripts/ipak.function.R")

required.packages = c("reshape2", "preprocessCore", "ggplot2")
required.bioconductor.packages = c("mygene")
ipak(required.packages = required.packages)
ibiopak(required.bioconductor.packages = required.bioconductor.packages)
  
# Load data (GSE78220)
GSE78220 = readRDS("./Data (Original)/GSE78220/GSE78220/data/GSE78220_normalized_matrix_Design.rds")

# Prepare Rdate file
Expression.matrix = GSE78220[[1]]
Clinical.data = GSE78220[[2]]

result = queryMany(qterms = rownames(Expression.matrix), scopes = "ensembl.gene")
results = data.frame(query = result$query, symbol = result$symbol) 
Expression.matrix = as.data.frame(Expression.matrix)
Expression.matrix$symbol = results$symbol[match(rownames(Expression.matrix), results$query)]
Expression.matrix = aggregate(. ~ symbol, Expression.matrix, FUN = mean)
rownames(Expression.matrix) = Expression.matrix$symbol
Expression.matrix$symbol = NULL
Expression.matrix = as.matrix(Expression.matrix)

dir.create("./Data/GSE78220", showWarnings = FALSE)
save(Expression.matrix, Clinical.data, file = "./Data/GSE78220/GSE78220.Rdata")

rm(list = ls())

## Load data (GSE115821)
load("./Data (Original)/GSE115821/GSE115821_Norm20180830.RData")
Expression.matrix = dataNorm
result = queryMany(qterms = rownames(Expression.matrix), scopes = "ensembl.gene")
results = data.frame(query = result$query, symbol = result$symbol) 
Expression.matrix = as.data.frame(Expression.matrix)
Expression.matrix$symbol = results$symbol[match(rownames(Expression.matrix), results$query)]
Expression.matrix = aggregate(. ~ symbol, Expression.matrix, FUN = mean)
rownames(Expression.matrix) = Expression.matrix$symbol
Expression.matrix$symbol = NULL
Expression.matrix = as.matrix(Expression.matrix)

Clinical.data = read.csv("./Data (Original)/GSE115821/MHH_GSE115821_annotation.csv", stringsAsFactors = FALSE)

Expression.matrix = log(Expression.matrix + 1, 2)

dir.create("./Data/GSE115821", showWarnings = FALSE)
save(Clinical.data, Expression.matrix, file = "./Data/GSE115821/GSE115821.Rdata")


rm(list = ls())
## Load data (GSE91061)
EDAS3 = readRDS("./Data (Original)/GSE91061/GSE91061/data/edas3.rds")
D4 = readRDS("./Data (Original)/GSE91061/GSE91061/data/d4.rds")
Clinical.data = D4

Expression.matrix = EDAS3
result = queryMany(qterms = rownames(Expression.matrix), scopes = "ensembl.gene")
results = data.frame(query = result$query, symbol = result$symbol) 
Expression.matrix = as.data.frame(Expression.matrix)
Expression.matrix$symbol = results$symbol[match(rownames(Expression.matrix), results$query)]
Expression.matrix = aggregate(. ~ symbol, Expression.matrix, FUN = mean)
rownames(Expression.matrix) = Expression.matrix$symbol
Expression.matrix$symbol = NULL
Expression.matrix = as.matrix(Expression.matrix)

dir.create("./Data/GSE91061", showWarnings = FALSE)
save(Expression.matrix, Clinical.data, file = "./Data/GSE91061/GSE91061.Rdata")

rm(list = ls())
## Load data (Van Allen)
matrix = readRDS("./Data (Original)/vanAllen_anti-CTLA4/data/matrix.rds")
design = readRDS("./Data (Original)/vanAllen_anti-CTLA4/data/design.rds")

Clinical.data = design
Expression.matrix = matrix

result = queryMany(qterms = rownames(Expression.matrix), scopes = "entrezgene")
results = data.frame(query = result$query, symbol = result$symbol) 
rownames(Expression.matrix) = results$symbol

dir.create("./Data/vanAllen", showWarnings = FALSE)
save(Clinical.data, Expression.matrix, file = "./Data/vanAllen/vanAllen_anti-CTLA4.Rdata")

rm(list = ls())
## Load data (Chen et al)
table_S6A = read.csv(file = "./Data (Original)/Chen_anti-PD1/supplementary_table_S6A.csv", stringsAsFactors = FALSE, header = TRUE)
table_S6A$Timepoint = gsub("./", "_", table_S6A$Timepoint)

combined = table_S6A
combined$Patient.ID[which(combined$Patient.ID == "6")] = c("6A", "6B", "6C")
rownames(combined) = paste(combined$Patient.ID, combined$Timepoint, sep = "_") 

Expression.matrix.prenorm = as.matrix(t(combined[,-c(1:5)]))
Clinical.data = combined[, c(1:5)]
Clinical.data$Patient.ID = table_S6A$Patient.ID

Expression.matrix = Expression.matrix.prenorm
# Quantile normalisation RNA
#Expression.matrix = normalize.quantiles(Expression.matrix.prenorm)                                                # Quantile normalize
#Expression.matrix = floor(Expression.matrix)                                                    # return non decimal values
#rownames(Expression.matrix) = rownames(Expression.matrix)
#colnames(Expression.matrix) = colnames(Expression.matrix)

dir.create("./Data/Chen", showWarnings = FALSE)
save(Clinical.data, Expression.matrix, file = "./Data/Chen/Chen.Rdata")

rm(list = ls())
# Load data Prat et al
Sheet1 = read.csv("./Data (Original)/Prat/Sheet1_Supplementary_table.csv", stringsAsFactors = FALSE)
Clinical.data = read.csv("./Data (Original)/Prat/Sheet2_Supplementary_table.csv", stringsAsFactors = FALSE)

Expression = Sheet1
rownames(Expression) = Expression$Gene.Name
Expression$Gene.Name = NULL
Expression$Accession = NULL
Expression$Class.Name = NULL
Expression.matrix.prenorm = as.matrix(Expression)

Housekeeping_genes_matrix = Expression.matrix.prenorm[which(Sheet1$Class.Name == "Housekeeping"),]
means_df = data.frame(Sample = colnames(Housekeeping_genes_matrix), Housekeeping_geomean = colMeans(Housekeeping_genes_matrix))
Geom_Mean_all = mean(Housekeeping_genes_matrix)

Expression.matrix_norm = Expression.matrix.prenorm

i=1
for (i in 1:ncol(Expression.matrix.prenorm)){
  sample = colnames(Expression.matrix.prenorm)[i]
  Geomean = means_df$Housekeeping_geomean[which(means_df$Sample == sample)]
  Expression.matrix_norm[,which(colnames(Expression.matrix_norm) == sample)] = (Expression.matrix.prenorm[,which(colnames(Expression.matrix.prenorm) == sample)] / Geomean) * Geom_Mean_all
  
}


Expression.matrix = Expression.matrix_norm

# Quantile normalisation RNA
#Expression.matrix = normalize.quantiles(Expression.matrix_norm)                                                # Quantile normalize
#Expression.matrix = floor(Expression.matrix)                                                    # return non decimal values
#rownames(Expression.matrix) = rownames(Expression.matrix_norm)
#colnames(Expression.matrix) = colnames(Expression.matrix_norm)

Expression.matrix = log(Expression.matrix, 2)

rownames(Clinical.data) = Clinical.data$RCC_FILES

dir.create("./Data/Prat", showWarnings = FALSE)
save(Clinical.data, Expression.matrix, file = "./Data/Prat/Prat.Rdata")

# Density plot
test = melt(Expression.matrix)
plot = ggplot(test, aes(x = value)) +
  geom_density(aes(group = Var2)) +
  theme_bw()
dev.new()
plot(plot)


