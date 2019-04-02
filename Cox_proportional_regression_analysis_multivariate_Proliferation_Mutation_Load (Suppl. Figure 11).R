
####################################################################
###
### Multivariate Cox proportional hazards regression analysis
### 
#####################################################################

# Before running this script, first download TCGA assembler 2.0.3 scripts http://www.compgenome.org/TCGA-Assembler/
# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg", "mi")
required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

source(paste0(code_path, "R tools/ggkm.R"))

# Set Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized_filtered"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq" 
ICR_k = "HML_classification"                                                                                            # "HML_classification" or "k3" or "k4" or "k5"
Surv_cutoff_years = 10
subset = "all"                              #"ICR_enabled", "ICR_disabled", "ICR_neutral", or "all"
exclude_medium = "exclude_medium"                       # For kaplan meijers: include medium, for multivariate regression analysis exclude
Source_surv_data = "Cell_paper"
Cutoff_HR = 1
Stage = "all"     #c("Stage I", "Stage II")   "all"              # Stage for filtering
Stagenames = "All stages"  #"Stage I and Stage II"              # "StageI and StageII" txt for in the file and directory names
Outcome = "OS"
Cancer = NA    # NA or "KIRC"

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/5.3.Hallmark_and_ICR_cluster_assignment_allcancers_SkippedDLBC.Rdata"))
load(paste0("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Signature_Enrichment/ssGSEA_Selected.pathways_ES.Rdata"))
#load("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Signature_Enrichment/Enabled_aggregated_ES_all.Rdata")
PanImmune_MS = read.csv("./3_DataProcessing/External/mmc2-PanImmune_MS.csv", stringsAsFactors = FALSE)
if(CancerTYPES == "ALL"){
  CancerTYPES = TCGA.cancersets$cancerType
}

Survival_df = read.csv("./2_Data/TCGA cell 2018 clinical/TCGA_CLINICAL_DATA_CELL_2018_S1.csv",
                       stringsAsFactors = FALSE)
Survival_df$ICR_cluster = Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster[match(Survival_df$bcr_patient_barcode,substring(rownames(Hallmark_and_ICR_cluster_assignment_allcancers), 1, 12))]
Survival_df = Survival_df[-which(is.na(Survival_df$ICR_cluster)),]
Survival_df$ED = NA

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/", Source_surv_data ,"_Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype %in% c("LAML", "DLBC")),]
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])

Survival_df$ED[which(Survival_df$type %in% ICR_enabled_cancers)] = "ICR_enabled"
Survival_df$ED[which(Survival_df$type %in% ICR_neutral_cancers)] = "ICR_neutral"
Survival_df$ED[which(Survival_df$type %in% ICR_disabled_cancers)] = "ICR_disabled"


if(subset == "all"){
  Survival_df = Survival_df
}else{Survival_df = Survival_df[which(Survival_df$ED == subset),]}

if(exclude_medium == "exclude_medium"){
  Survival_df = Survival_df[-which(Survival_df$ICR_cluster == "ICR Medium"),]
}

if(!is.na(Cancer)){
  Survival_df = Survival_df[which(Survival_df$type == Cancer),]
}

Survival_df$ajcc_pathologic_tumor_stage[which(Survival_df$ajcc_pathologic_tumor_stage %in% c("[Discrepancy]", "[Not Applicable]", "[Not Available]", "[Unknown]"))] = NA
Survival_df$ajcc_pathologic_tumor_stage[which(Survival_df$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA", "Stage IB"))] = "Stage I"
Survival_df$ajcc_pathologic_tumor_stage[which(Survival_df$ajcc_pathologic_tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC"))] = "Stage II"
Survival_df$ajcc_pathologic_tumor_stage[which(Survival_df$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"))] = "Stage III"
Survival_df$ajcc_pathologic_tumor_stage[which(Survival_df$ajcc_pathologic_tumor_stage %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC"))] = "Stage IV"
Survival_df$ajcc_pathologic_tumor_stage = factor(Survival_df$ajcc_pathologic_tumor_stage, levels = c("I/II NOS", "IS", "Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV", "Stage X"))

Survival_df$Proliferation = ES.all[which(rownames(ES.all) == "[LM] Proliferation"),][match(Survival_df$bcr_patient_barcode, substring(colnames(ES.all), 1, 12))]
#Survival_df$Proliferation = ES.all[which(rownames(ES.all) == "enabling_pathways"),][match(Survival_df$bcr_patient_barcode, substring(colnames(ES.all), 1, 12))]

Survival_df$TGF_beta_sig = ES.all[which(rownames(ES.all) == "[HM] TGF beta signaling"),][match(Survival_df$bcr_patient_barcode, substring(colnames(ES.all), 1, 12))]
Survival_df$Mutation_rate = PanImmune_MS$Nonsilent.Mutation.Rate[match(Survival_df$bcr_patient_barcode, PanImmune_MS$TCGA.Participant.Barcode)]
Survival_df$Mutation_rate = log10(Survival_df$Mutation_rate + 1)
Survival_df$Neoantigen_load = PanImmune_MS$SNV.Neoantigens[match(Survival_df$bcr_patient_barcode, PanImmune_MS$TCGA.Participant.Barcode)]
Survival_df$Neoantigen_load = log10(Survival_df$Neoantigen_load + 1)


Y = Surv_cutoff_years * 365
TS.Alive = Survival_df[Survival_df[, Outcome] == "0", c(Outcome,  paste0(Outcome, ".time"), "ICR_cluster", "ajcc_pathologic_tumor_stage", "type", "Proliferation", "TGF_beta_sig","Mutation_rate", "Neoantigen_load")]
colnames(TS.Alive) = c("Status","Time", "ICR_cluster", "pathologic_stage", "Cancer", "Proliferation", "TGF_beta_sig", "Mutation_rate", "Neoantigen_load")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Survival_df[Survival_df[, Outcome] == "1", c(Outcome,  paste0(Outcome, ".time"), "ICR_cluster", "ajcc_pathologic_tumor_stage", "type", "Proliferation", "TGF_beta_sig", "Mutation_rate", "Neoantigen_load")]
colnames(TS.Dead) = c("Status","Time", "ICR_cluster", "pathologic_stage", "Cancer", "Proliferation", "TGF_beta_sig", "Mutation_rate", "Neoantigen_load")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "0"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "1"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                                                                                         # remove patients with less then 1 day follow up time

TS.Surv$Proliferation = as.numeric(TS.Surv$Proliferation)
TS.Surv$Mutation_rate = as.numeric(TS.Surv$Mutation_rate)
TS.Surv$Neoantigen_load = as.numeric(TS.Surv$Neoantigen_load)
TS.Surv$TGF_beta_sig = as.numeric(TS.Surv$TGF_beta_sig)

# Final filter for Stage
if(sum(Stage %in% TS.Surv$pathologic_stage)>=1){
  TS.Surv = TS.Surv[which(TS.Surv$pathologic_stage %in% Stage),]
}else{print(paste0("For ", subset, " no patients with ", Stagenames, " available")) 
  next}

#TS.Surv[,"Group"] = factor(TS.Surv[,"Group"], levels = c("ICR High", "ICR Medium", "ICR Low"))
#TS.Surv[,"ICR_cluster"] = factor(TS.Surv[,"ICR_cluster"], levels = c("ICR High", "ICR Medium", "ICR Low"))

# Multi-variate
#multivariate = coxph(formula = Surv(Time, Status) ~ ICR_cluster + pathologic_stage, data = TS.Surv)
#summary(multivariate)

#Uni-variate
#uni_variate_ICR = coxph(formula = Surv(Time, Status) ~ ICR_cluster, data = TS.Surv)
#summary(uni_variate_ICR)

# Lance miller approach: "Semi-continuous: 1, 2, 3

# stage I or II NOS (T = TX, T2, or T3 / N = N0 / M = M0), for which TNM staging was incomplete 
TS.Surv$pathologic_stage = as.character(TS.Surv$pathologic_stage)
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "I/II NOS")] = NA
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "IS")] = NA
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage 0")] = NA
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage I")] = 1
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage II")] = 2
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage III")] = 3
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage IV")] = 4
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage X")] = NA

TS.Surv$pathologic_stage = as.numeric(TS.Surv$pathologic_stage)
TS.Surv$ICR_cluster = factor(TS.Surv$ICR_cluster, levels = c("ICR High", "ICR Low")) # adjust this by hard coding when "ICR Medium" is included

TS.Surv = TS.Surv[-which(is.na(TS.Surv$Mutation_rate)),]

#TS.Surv = TS.Surv[which(TS.Surv$Cancer %in% ICR_neutral_cancers),]

#Uni-variate
uni_variate_ICR = coxph(formula = Surv(Time, Status) ~ ICR_cluster, data = TS.Surv)
summary(uni_variate_ICR)

uni_variate_proliferation = coxph(formula = Surv(Time, Status) ~ Proliferation, data = TS.Surv)
summary(uni_variate_proliferation)

uni_variate_TGF_beta_sig = coxph(formula = Surv(Time, Status) ~ TGF_beta_sig, data = TS.Surv)
summary(uni_variate_TGF_beta_sig)

uni_variate_mutation = coxph(formula = Surv(Time, Status) ~ Mutation_rate, data = TS.Surv)
summary(uni_variate_mutation)

uni_variate_ps = coxph(formula = Surv(Time, Status) ~ pathologic_stage, data = TS.Surv)
summary(uni_variate_ps)

uni_variate_age = coxph(formula = Surv(Time, Status) ~ age, data = TS.Surv)

# Multi-variate
multivariate = coxph(formula = Surv(Time, Status) ~ ICR_cluster + Proliferation + Mutation_rate + TGF_beta_sig, data = TS.Surv)
summary(multivariate)

multivariate = coxph(formula = Surv(Time, Status) ~ ICR_cluster + Mutation_rate + Proliferation, data = TS.Surv)
summary(multivariate)

class_semi_continuous = missing_variable(TS.Surv$pathologic_stageT, type = "semi-continuous")


###semi-continuous-class

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                                                                                    # calculate the number of months
mfit = survfit(msurv~TS.Surv$ICR_cluster,conf.type = "log-log")

# Calculations
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

# Check this!!
##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR High")
mHR = coxph(formula = msurv ~ TS.Surv[,"ICR_cluster"],data = TS.Surv, subset = TS.Surv$ICR_cluster %in% c("ICR High", "ICR Low"))
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

PLOT_P = signif(p[2],3)
PLOT_HR = round(signif(exp(mHR.extract@coef),3)[2], 3)
PLOT_CI1 = CI[2,1]
PLOT_CI2 = CI[2,2]

dir.create(paste0("./5_Figures/Pancancer_plots/Assembler_Panca_Normalized_filtered/Survival_Plots/Benefit_clusters"), showWarnings = FALSE)

png(paste0("./5_Figures/Pancancer_plots/Assembler_Panca_Normalized/Survival_Plots/Benefit_clusters/",
           "Kaplan_Meier_", Stagenames, "_", subset, "_", Cancer, "_samples.png"),
    res=600,height=6,width=8,unit="in")                                                                                           # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv$ICR_cluster),
     ystrataname = NULL,
     main= paste0("Survival curve across ICR groups (", ICR_k, ") in ", subset, " in ",
                  Stagenames),
     xlabs = "Time in months",
     cbPalette = cbPalette,
     PLOT_HR = PLOT_HR,
     PLOT_P = PLOT_P,
     PLOT_CI1 = PLOT_CI1,
     PLOT_CI2 = PLOT_CI2)
dev.off()

contingency_table = table(patient.table.all$ICR_cluster, patient.table.all$pathologic_stage)
chi2 = chisq.test(contingency_table)
chi2$p.value