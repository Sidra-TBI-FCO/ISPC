#################################################################
###
###
### Data input:
###
### Output :
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

required.packages <- c("ggplot2", "graphics")
ipak(required.packages)

# Set Parameters
download.method = "Assembler_Panca_Normalized_filtered"                                                                                       # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"
Source_surv_data = "Cell_paper"
Cutoff_HR = 1

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers_Skipped_DLBC.Rdata"))
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/ssGSEA_Selected.pathways_ES.Rdata"))
load(paste0("./3_DataProcessing/", download.method,"/Pancancer/Pancancer_gene_RNAseq_normalized_Tissue_Type_and_whitelists_filtered.Rdata"))
PanImmune_MS = read.csv("./3_DataProcessing/External/mmc2-PanImmune_MS.csv", stringsAsFactors = FALSE)


ICR_cluster_assignment_allcancers[,"Proliferation"] = ES.all["[LM] Proliferation",][match(rownames(ICR_cluster_assignment_allcancers), colnames(ES.all))]
ICR_cluster_assignment_allcancers[,"TGF beta sig"] = ES.all["[HM] TGF beta signaling",][match(rownames(ICR_cluster_assignment_allcancers), colnames(ES.all))]
ICR_cluster_assignment_allcancers[,"Mutation Rate"] = PanImmune_MS$Nonsilent.Mutation.Rate[match(substring(rownames(ICR_cluster_assignment_allcancers), 1, 12),
                                                                                                 PanImmune_MS$TCGA.Participant.Barcode)]

ICR_cluster_assignment_allcancers = ICR_cluster_assignment_allcancers[-which(ICR_cluster_assignment_allcancers$`Mutation Rate` == 0),]
ICR_cluster_assignment_allcancers = ICR_cluster_assignment_allcancers[-which(is.na(ICR_cluster_assignment_allcancers$`Mutation Rate`)),]

ICR_cluster_assignment_allcancers$`Mutation Rate` = log10(ICR_cluster_assignment_allcancers$`Mutation Rate`)

ICR_cluster_assignment_allcancers$ED = NA
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/", Source_surv_data ,"_Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype %in% c("LAML", "DLBC")),]
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])

ICR_cluster_assignment_allcancers$ED[which(ICR_cluster_assignment_allcancers$Cancer %in% ICR_enabled_cancers)] = "ICR_enabled"
ICR_cluster_assignment_allcancers$ED[which(ICR_cluster_assignment_allcancers$Cancer %in% ICR_neutral_cancers)] = "ICR_neutral"
ICR_cluster_assignment_allcancers$ED[which(ICR_cluster_assignment_allcancers$Cancer %in% ICR_disabled_cancers)] = "ICR_disabled"

ICR_cluster_assignment_allcancers$ED = factor(ICR_cluster_assignment_allcancers$ED, 
                                              levels = c("ICR_enabled", "ICR_neutral",
                                                         "ICR_disabled"))
reg <- function(x, y, col) abline(lm(y~x), col=col)

upper.panel<-function(x, y, col.smooth = "red"){
  #points(x,y, pch=19, cex = 0.1, col=c("orange", "grey", "purple")[ICR_cluster_assignment_allcancers$ED])
  points(x, y, pch=19, cex = 0.1, col = "black")
  r <- round(cor(x, y), digits=2)
  txt <- paste0("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.5, 0.9, txt)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) reg(x[ok], y[ok], col.smooth)
}

#####
reg <- function(x, y, col) abline(lm(y~x), col=col) 

panel.lm =  function (x, y, col = par("col"), bg = NA, pch = par("pch"), 
                      cex = 0.1, col.smooth = "red", span = 2/3, iter = 3, ...)  {
  points(x, y, pch = pch, col = col, bg = bg, cex = cex)
  ok <- is.finite(x) & is.finite(y)
  if (any(ok)) reg(x[ok], y[ok], col.smooth)
}

panel.cor <- function(x, y, digits = 2, prefix = "", cex.cor, ...)
{
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  r <- abs(cor(x, y))
  txt <- format(c(r, 0.123456789), digits = digits)[1]
  txt <- paste0(prefix, txt)
  text(0.5, 0.5, txt, cex = 1.1, font = 4)
}

dev.new()
pairs(ICR_cluster_assignment_allcancers[, c("ICRscore", "Proliferation", "TGF beta sig", "Mutation Rate")], 
      upper.panel = panel.lm,
      lower.panel = panel.cor)


######
i=1
for (i in 1:N.Pathways){
  pathway = Pathways[i]
  ICR_cluster_assignment_allcancers[,target] = ES.all[target,][match(rownames(ICR_cluster_assignment_allcancers), colnames(ES.all))]
  ICR_cluster_assignment_allcancers$ES = ES.all[pathway,][match(rownames(ICR_cluster_assignment_allcancers), colnames(ES.all))]
  dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/", target, "_Pathways_ES_Correlation_Plots/"), showWarnings = FALSE)
  plot =  ggplot(ICR_cluster_assignment_allcancers, aes(x = get(target), y = ES)) +
    geom_point(size = 0.5, shape = 1) +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    geom_smooth(method=lm, color="blue") +
    ylab(paste0(pathway, " ES")) +
    xlab(paste0(target))
  png(paste0("./5_Figures/Pancancer_plots/", download.method,"/", target, "_Pathways_ES_Correlation_Plots/", pathway,"_", target, "_Correlation_plot.png"),
      res = 600, height = 4, width = 4, units = "in")
  plot(plot)
  dev.off()
}


#dev.new()
#plot(ICR_cluster_assignment_allcancers$ICRscore, ICR_cluster_assignment_allcancers$TISscore)
#cor.test(x = ICR_cluster_assignment_allcancers$ICRscore, y = ICR_cluster_assignment_allcancers$`Cytolytic activity score`, method = "pearson")

