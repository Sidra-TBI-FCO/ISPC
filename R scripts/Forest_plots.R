

## Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("RColorBrewer", "forestplot")
ipak(required.packages)
ibiopak(required.bioconductor.packages)

#Set parameters
clustering = "per_cancer"

HR.table = read.csv(file = paste0("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/1.Pancancer_ICR_Clustering/HR_table_", clustering,"_ICR_clusters.csv"),
                    stringsAsFactors = FALSE)

# Cochrane data from the 'rmeta'-package
cochrane_from_rmeta =
  structure(list(
    mean  = as.numeric(c(NA,HR.table$HR[1:4]), NA),
    lower = c(NA,HR.table$CI_lower[c(1:4)], NA),
    upper = c(NA,HR.table$CI_upper[c(1:4)], NA)),
    .Names = c("mean", "lower", "upper"),
    row.names = c(NA, -4),
    class = "data.frame")

tabletext<-cbind(
  c("Cancer", as.character(HR.table$Cancertype)[c(1:4)]),
  c("N", HR.table$N[c(1:4)]),
  c("p-value", HR.table$p_value[c(1:4)]),
  c("HR",      HR.table$HR[c(1:4)]))

#dev.new()

pdf(file = paste0("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/1.Pancancer_ICR_Clustering/", clustering,
"_forestplot.pdf"),
    height = 2, width = 8)

forestplot(tabletext,
           cochrane_from_rmeta,new_page = FALSE,
           is.summary=c(TRUE,rep(FALSE,4),TRUE,rep(FALSE,4),TRUE,FALSE),
           #clip=c(0.1,1),
           xticks = c(0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0),
           xlog=TRUE,
           #xlim = c(0.010, 1.0),
           #alim=c(0,1),
           boxsize = .25,
           vertices = TRUE,
           col=fpColors(box="royalblue",line="darkblue", summary="royalblue"),
           txt_gp = fpTxtGp(label = gpar(fontsize = 15), xlab = gpar(fontsize = 24), cex = 1))
dev.off()
