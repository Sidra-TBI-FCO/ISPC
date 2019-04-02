####################################################################
###
### This Script calculates 
### 
### Input data:
### ("./3_DataProcessing/",download.method,"/",Cancer,"/SurvivalData/")
### Output data are saved as Rdata file:
#####################################################################

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg")
required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

source(paste0(code_path, "R tools/ggkm.R"))

# Set Parameters
Pathways = "ALL"
CancerTYPES = "ALL"   # Don't change this variable
Pathway_skip = ""                                                                                                      # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized_filtered"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq" 
Log_file = paste0("./1_Log_Files/", download.method ,"/5.9.Pancancer_Survival_Analysis/5.9.Pancancer_Survival_Analysis_Log_File_",                              # Specify complete name of the logfile that will be saved during this script
                  gsub(":",".",gsub(" ","_",date())),".txt")
ICR_k = "HML_classification"                                                                                            # "HML_classification" or "k3" or "k4" or "k5"
Surv_cutoff_years = 10
Cutoff_HR = 1
ICR_benefit_cluster = "all"
ICR_ID_cancer_type = "all"
Cancer_to_select = "all"
Source_surv_data = "Cell_paper"
Outcome = "OS"
version = "5.10."            # "5.3." (Pancancer) or "3.17." (Per cancer) or "5.10." (t.test method)

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                   # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 
load(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis/", Outcome , "_", Source_surv_data, "_Survival_analysis_High_vs_Low_Groups", ICR_k,".Rdata"))  
load(paste0("./4_Analysis/", download.method,"/Pan_Cancer/Clustering/Try_out_MAPK_", version, "and.5.8.Hallmark_and_ICR_cluster_assignment_allcancers_ICRbenefit_cluster.Rdata"))

if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

# Create folders
dir.create("./4_Analysis/",showWarnings = FALSE)                                                                        # Create folder to save processed data (by Assembler module B)
dir.create(paste0("./4_Analysis/",download.method),showWarnings = FALSE)
dir.create(paste0("./1_Log_Files/"), showWarnings = FALSE)                                                              # Create folder to save logfile
dir.create(paste0("./1_Log_Files/", download.method), showWarnings = FALSE)  
dir.create(paste0("./1_Log_Files/", download.method ,"/5.9.Pancancer_Survival_Analysis/"), showWarnings = FALSE)
cat("This is a log file for Survival Analysis of ",                                                                     # Set-up logfile
    "__________________________________________",
    "",
    "Session Info :",
    capture.output(sessionInfo()),
    "",
    "Script Running Date :",
    capture.output(Sys.time()),
    "",
    "Parameters Used :",
    paste0("CancerTYPES = ", CancerTYPES),                                                          
    paste0("Pathway_skip = ", Pathway_skip),
    paste0("download.method = ", download.method),
    "",
    "Scripts output :",
    "",
    "Clustering",
    file = Log_file,
    append = FALSE, sep= "\n")

#load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = ICR_disabled_cancers[-which(ICR_disabled_cancers == "LAML")]
ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])
ICR_enabled_cancer_samples = rownames(Hallmark_and_ICR_cluster_assignment_allcancers)[which(Hallmark_and_ICR_cluster_assignment_allcancers$Cancer %in% ICR_enabled_cancers)]
ICR_disabled_cancer_samples = rownames(Hallmark_and_ICR_cluster_assignment_allcancers)[which(Hallmark_and_ICR_cluster_assignment_allcancers$Cancer %in% ICR_disabled_cancers)]
ICR_neutral_cancer_samples = rownames(Hallmark_and_ICR_cluster_assignment_allcancers)[which(Hallmark_and_ICR_cluster_assignment_allcancers$Cancer %in% ICR_neutral_cancers)]

if(Source_surv_data == "TCGA_Assembler"){
  ICR_survival_data = read.csv(paste0("./3_DataProcessing/TCGA_Assembler/ACC/SurvivalData/updatedsurvivaldata.csv"))
  for (i in 2:length(CancerTYPES)){
    Cancer = CancerTYPES[i]
    if(Cancer == "LAML"){next}
    survival_data = read.csv(paste0("./3_DataProcessing/TCGA_Assembler/", Cancer, "/SurvivalData/updatedsurvivaldata.csv"))
    ICR_survival_data = rbind(ICR_survival_data, survival_data)
  }
}

if(Source_surv_data == "Cell_paper"){
  ICR_survival_data = read.csv("./2_Data/TCGA cell 2018 clinical/TCGA_CLINICAL_DATA_CELL_2018_S1.csv",
                               stringsAsFactors = FALSE)
  colnames(ICR_survival_data)[which(colnames(ICR_survival_data) == "death_days_to")] = "days_to_death"
  colnames(ICR_survival_data)[which(colnames(ICR_survival_data) == "last_contact_days_to")] = "days_to_last_followup"
}


if(ICR_benefit_cluster == "all"){
  Survival_data = ICR_survival_data
}else{
  selected_patients = rownames(Hallmark_and_ICR_cluster_assignment_allcancers)[which(Hallmark_and_ICR_cluster_assignment_allcancers$ICR_benefit_cluster == ICR_benefit_cluster)]
  Survival_data = ICR_survival_data[which(ICR_survival_data$bcr_patient_barcode %in% substring(selected_patients, 1, 12)),]
}

if(ICR_ID_cancer_type == "all"){
  Survival_data = Survival_data
}else{
  samples_to_select = get(paste0("ICR_", ICR_ID_cancer_type, "_cancer_samples"))
  Survival_data = Survival_data[which(Survival_data$bcr_patient_barcode %in% substring(samples_to_select, 1, 12)),]
}

if(Cancer_to_select == "ALL"){
  Survival_data = Survival_data
}else{
  samples_to_select2 = rownames(Hallmark_and_ICR_cluster_assignment_allcancers)[which(Hallmark_and_ICR_cluster_assignment_allcancers$Cancer == Cancer_to_select)]
  Survival_data = Survival_data[which(Survival_data$bcr_patient_barcode %in% substring(samples_to_select2, 1, 12)),]
}

# Create folders to save the data
dir.create(paste0("./4_Analysis/",download.method,"/Pan_Cancer"),showWarnings = FALSE)
dir.create(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis"),showWarnings = FALSE)
dir.create(paste0("./5_Figures"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/"), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method), showWarnings = FALSE)
dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Try_out_MAPK_Survival_Plots_ICR_benefit_cluster_", version), showWarnings = FALSE)

Survival_data$ICR_cluster = Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster[match(Survival_data$bcr_patient_barcode,substring(rownames(Hallmark_and_ICR_cluster_assignment_allcancers),1,12))]
Survival_data = Survival_data[!is.na(Survival_data$ICR_cluster),]
Survival_data$ICR_cluster = factor(Survival_data$ICR_cluster, levels = c("ICR High", "ICR Medium", "ICR Low")) 
Highest_ICR_group = "ICR High"

Y = Surv_cutoff_years * 365
TS.EventFree = Survival_data[Survival_data[, Outcome] == "0", c(Outcome, paste0(Outcome, ".time"), "ICR_cluster")]
colnames(TS.EventFree) = c("Status","Time", "Group")
TS.EventFree$Time = as.numeric(as.character(TS.EventFree$Time))
TS.EventFree$Time[TS.EventFree$Time > Y] = Y

TS.EventOccured = Survival_data[Survival_data[, Outcome] == "1", c(Outcome, paste0(Outcome, ".time"), "ICR_cluster")]
colnames(TS.EventOccured) = c("Status","Time", "Group")
TS.EventOccured$Time = as.numeric(as.character(TS.EventOccured$Time))
TS.EventOccured$Status[which(TS.EventOccured$Time> Y)] = "EventFree"
TS.EventOccured$Time[TS.EventOccured$Time > Y] = Y

TS.Surv = rbind (TS.EventOccured,TS.EventFree)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "1"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                                                                                         # remove patients with less then 1 day follow up time

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                                                                                    # calculate the number of months
mfit = survfit(msurv~TS.Surv$Group,conf.type = "log-log")

# Calculations
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

#TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR High", "ICR Medium", "ICR Low"))
TS.Surv[,"Group"] = as.factor(TS.Surv[,"Group"])

# Check this!!
##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR High")
mHR = coxph(formula = msurv ~ TS.Surv[,"Group"],data = TS.Surv, subset = TS.Surv$Group %in% c("ICR High", "ICR Low"))
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

PLOT_P = signif(p[2], digits = 3)
PLOT_HR = round(signif(exp(mHR.extract@coef),3)[2], 3)
PLOT_CI1 = CI[2,1]
PLOT_CI2 = CI[2,2]

# plots
png(paste0("./5_Figures/Pancancer_plots/", download.method, "/Try_out_MAPK_Survival_Plots_ICR_benefit_cluster_", version, "/",
           "/Kaplan_Meier_ICR_clusters_", Cancer_to_select ,"_ICR_benefit_cluster_" ,"_", ICR_benefit_cluster, "_", ICR_ID_cancer_type, "_Pancancer.png"),res=600,height=6,width=8,unit="in")                                                                                           # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv[,"Group"]),
     ystrataname = NULL,
     main= paste0("Survival curve across ICR groups of ", Cancer_to_select, " ", ICR_ID_cancer_type, " cancertypes in ", ICR_benefit_cluster, " samples"),
     xlabs = "Time in months",
     cbPalette = cbPalette,
     PLOT_HR = PLOT_HR,
     PLOT_P = PLOT_P,
     PLOT_CI1 = PLOT_CI1,
     PLOT_CI2 = PLOT_CI2)
dev.off()
