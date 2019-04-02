#################################################################
###
###
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

required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg", "stringr")
required.bioconductor.packages = "survival"
#ibiopak("")
ipak(required.packages)
ibiopak(required.bioconductor.packages)

source(paste0(code_path, "R tools/ggkm.R"))

# Set Parameters
CancerTYPES = "ALL"
download.method = "Assembler_Panca_Normalized_filtered"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
Surv_cutoff_years = 10
Source_surv_data = "Cell_paper"
Outcome = "OS"
variable_of_interest = "Nonsilent.Mutation.Rate" #SNV.Neoantigens #Nonsilent.Mutation.Rate
Variable_Category_mut = "Nonsilent.Mutation.Rate Low"
pathway_of_interest = "Proliferation"
Variable_Category_ES = "High"

# Load data
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Mutational_load/6.3.2.", variable_of_interest, "_Pancancer_Categories.Rdata"))
Survival_data = read.csv("./2_Data/TCGA cell 2018 clinical/TCGA_CLINICAL_DATA_CELL_2018_S1.csv",
                         stringsAsFactors = FALSE)
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/5.3.and.Suppl.Enabled_Aggregated_Hallmark_and_ICR_cluster_assignment_allcancers_SkippedDLBC.Rdata"))
colnames(Hallmark_and_ICR_cluster_assignment_allcancers) = gsub(".*] ", "",colnames(Hallmark_and_ICR_cluster_assignment_allcancers))
colnames(Hallmark_and_ICR_cluster_assignment_allcancers) = gsub("\\/", "_",colnames(Hallmark_and_ICR_cluster_assignment_allcancers))
colnames(Hallmark_and_ICR_cluster_assignment_allcancers)[which(colnames(Hallmark_and_ICR_cluster_assignment_allcancers) == "Phopholipase")] = "Phospholipase"
colnames(Hallmark_and_ICR_cluster_assignment_allcancers) = gsub("_cluster_Pancancer", "", colnames(Hallmark_and_ICR_cluster_assignment_allcancers))

Survival_data$X = NULL
df_plot$Pathway_category = Hallmark_and_ICR_cluster_assignment_allcancers[,which(colnames(Hallmark_and_ICR_cluster_assignment_allcancers) == pathway_of_interest)][match(df_plot$sample_barcodes,
                                                                                                                                                                         rownames(Hallmark_and_ICR_cluster_assignment_allcancers))]
df_plot$Pathway_category = word(df_plot$Pathway_category, -1)                                                                                                                                                                  
df_plot$bcr_patient_barcode = substring(df_plot$sample_barcodes, 1, 12)
Survival_data = merge(Survival_data, df_plot, by = "bcr_patient_barcode")
Survival_data$ICR_cluster = factor(Survival_data$ICR_cluster, levels = c("ICR High", "ICR Medium", "ICR Low")) 

Y = Surv_cutoff_years * 365
TS.EventFree = Survival_data[Survival_data[, Outcome] == "0", c(Outcome, paste0(Outcome, ".time"), "ICR_cluster", "Category", "Pathway_category")]
colnames(TS.EventFree) = c("Status","Time", "ICR_cluster", "Variable", "Pathway_category")
TS.EventFree$Time = as.numeric(as.character(TS.EventFree$Time))
TS.EventFree$Time[TS.EventFree$Time > Y] = Y

TS.EventOccured = Survival_data[Survival_data[, Outcome] == "1", c(Outcome, paste0(Outcome, ".time"), "ICR_cluster", "Category", "Pathway_category")]
colnames(TS.EventOccured) = c("Status","Time", "ICR_cluster", "Variable", "Pathway_category")
TS.EventOccured$Time = as.numeric(as.character(TS.EventOccured$Time))
TS.EventOccured$Status[which(TS.EventOccured$Time> Y)] = "0"
TS.EventOccured$Time[TS.EventOccured$Time > Y] = Y

TS.Surv = rbind (TS.EventOccured,TS.EventFree)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "1"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)

# Final filter for variable mutation category
TS.Surv = TS.Surv[which(TS.Surv$Variable == Variable_Category_mut),]
if(length(unique(TS.Surv$ICR_cluster)) == 1){
  print(paste0("For ", Cancer, " all patients of ", variable_of_interest, " ", Variable_Category, " are ", unique(TS.Surv$ICR_cluster)))
  next
}
if(length(unique(TS.Surv$Status)) == 1){
  print(paste0("For ", Cancer, " all patients of ",variable_of_interest, " ", Pathway_Category, " are survival status ", unique(TS.Surv$Status)))
  next
}

# Final filter for variable pathway category
TS.Surv = TS.Surv[which(TS.Surv$Pathway_category == Variable_Category_ES),]

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                                                                                    # calculate the number of months
mfit = survfit(msurv~TS.Surv$ICR_cluster,conf.type = "log-log")

# Calculations
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

#TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR High", "ICR Medium", "ICR Low"))
#TS.Surv = as.factor(TS.Surv[,"Group"])

# Check this!!
##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR High")
mHR = coxph(formula = msurv ~ TS.Surv$ICR_cluster,data = TS.Surv, subset = TS.Surv$ICR_cluster %in% c("ICR High", "ICR Low"))
mHR.extract = extract.coxph(mHR, include.aic = TRUE,
                            include.rsquared = TRUE, include.maxrs=TRUE,
                            include.events = TRUE, include.nobs = TRUE,
                            include.missings = TRUE, include.zph = TRUE)
HRtxt = paste("Hazard-ratio =", signif(exp(mHR.extract@coef),3),"for",names(mHR$coefficients))
beta = coef(mHR)
se   = sqrt(diag(mHR$var))
p    = 1 - pchisq((beta/se)^2, 1)
CI   = confint(mHR)
CI   = round(exp(CI),2)

PLOT_P = p[2]
PLOT_HR = round(signif(exp(mHR.extract@coef),3)[2], 3)
PLOT_CI1 = CI[2,1]
PLOT_CI2 = CI[2,2]      

dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Mutation/Kaplan_Meier_Mutational_load_Pancancer"), showWarnings = FALSE)
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Mutation/Kaplan_Meier_Mutational_load_Pancancer/Kaplan_Meier_", variable_of_interest, "_", Variable_Category_mut, "_",
           pathway_of_interest, "_", Variable_Category_ES,
           "_", Source_surv_data, "_", Outcome , "_ICR_clusters.png"), res=600,height=6,width=8,unit="in")  

ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv$ICR_cluster),
     ystrataname = NULL,
     main= paste0("Survival curve across ICR groups in all cancer samples", "\n of ", variable_of_interest, " ", Variable_Category_mut, "\n and ",
                  pathway_of_interest, " ", Variable_Category_ES),
     xlabs = "Time in months",
     cbPalette = cbPalette,
     PLOT_HR = PLOT_HR,
     PLOT_P = PLOT_P,
     PLOT_CI1 = PLOT_CI1,
     PLOT_CI2 = PLOT_CI2)
dev.off()
