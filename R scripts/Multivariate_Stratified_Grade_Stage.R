
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

# Change this parameter!
subset = "ICR_neutral"                              #"ICR_enabled", "ICR_disabled", "ICR_neutral", or "all"
MC = "MC"

# Fixed Parameters
CancerTYPES = "ALL"                                                                                                     # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                        # If CancerTYPES = "ALL", specify here if you want to skip cancertypes
download.method = "Assembler_Panca_Normalized_filtered"                                                                                      # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq" 
ICR_k = "HML_classification"                                                                                            # "HML_classification" or "k3" or "k4" or "k5"
Surv_cutoff_years = 10
exclude_medium = "exclude_medium"                       # For kaplan meijers: include medium, for multivariate regression analysis exclude
Source_surv_data = "Cell_paper"
Cutoff_HR = 1
Stage = "all"     #c("Stage I", "Stage II")   "all"              # Stage for filtering
Stagenames = "All stages"  #"Stage I and Stage II"              # "StageI and StageII" txt for in the file and directory names
Outcome = "OS"
Cancer = NA  #c("LGG", "GBM")     # NA or "KIRC"

# Load data
#load(paste0(code_path, "Datalists/ICR_genes.RData")) 
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/5.3.Hallmark_and_ICR_cluster_assignment_allcancers_SkippedDLBC.Rdata"))
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/ssGSEA_Selected.pathways_ES.Rdata"))
if(CancerTYPES == "ALL"){
  CancerTYPES = TCGA.cancersets$cancerType
}
PanImmune_MS = read.csv("./3_DataProcessing/External/mmc2-PanImmune_MS.csv", stringsAsFactors = FALSE)

Aneuploidy_df = read.csv("./3_DataProcessing/External/mmc2-Aneuploidy.csv", stringsAsFactors = FALSE, header = TRUE)
Survival_df = read.csv("./2_Data/TCGA cell 2018 clinical/TCGA_CLINICAL_DATA_CELL_2018_S1.csv",
                       stringsAsFactors = FALSE)
Survival_df$ICR_cluster = Hallmark_and_ICR_cluster_assignment_allcancers$HML_cluster[match(Survival_df$bcr_patient_barcode,substring(rownames(Hallmark_and_ICR_cluster_assignment_allcancers), 1, 12))]
Survival_df = Survival_df[-which(is.na(Survival_df$ICR_cluster)),]
Survival_df$ED = NA
Survival_df$Mutation_rate = PanImmune_MS$Nonsilent.Mutation.Rate[match(Survival_df$bcr_patient_barcode, PanImmune_MS$TCGA.Participant.Barcode)]
Survival_df$Mutation_rate = log10(Survival_df$Mutation_rate + 0.0001)

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Survival_Analysis/", Source_surv_data ,"_Survival_analysis_High_vs_Low_GroupsHML_classification.Rdata"))
All_survival_analysis_data = All_survival_analysis_data[-which(All_survival_analysis_data$Cancertype %in% c("LAML", "DLBC")),]
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])

Survival_df$ED[which(Survival_df$type %in% ICR_enabled_cancers)] = "ICR_enabled"
Survival_df$ED[which(Survival_df$type %in% ICR_neutral_cancers)] = "ICR_neutral"
Survival_df$ED[which(Survival_df$type %in% ICR_disabled_cancers)] = "ICR_disabled"

Survival_df$TGF_beta_ES = ES.all["[HM] TGF beta signaling",][match(Survival_df$bcr_patient_barcode, 
                                                                   substring(colnames(ES.all), 1, 12))]

Survival_df$Proliferation_ES = ES.all["[LM] Proliferation",][match(Survival_df$bcr_patient_barcode, 
                                                                   substring(colnames(ES.all), 1, 12))]

Survival_df$Aneuploidy_score = Aneuploidy_df$AneuploidyScore.AS.[match(Survival_df$bcr_patient_barcode, substring(Aneuploidy_df$Sample, 1, 12))]

if(MC == "MC"){
  cancer.types = unique(Survival_df$type)
  for (j in cancer.types) {
    cancer.scaled1 = as.numeric(scale(Survival_df$Proliferation_ES[Survival_df$type == j],center = TRUE))
    Survival_df$Proliferation_ES[Survival_df$type==j] = cancer.scaled1
    
    cancer.scaled2 = as.numeric(scale(Survival_df$TGF_beta_ES[Survival_df$type == j],center = TRUE))
    Survival_df$TGF_beta_ES[Survival_df$type==j] = cancer.scaled2
    
    cancer.scaled3 = as.numeric(scale(Survival_df$Aneuploidy_score[Survival_df$type == j],center = TRUE))
    Survival_df$Aneuploidy_score[Survival_df$type==j] = cancer.scaled3
    
    cancer.scaled4 = as.numeric(scale(Survival_df$Mutation_rate[Survival_df$type == j],center = TRUE))
    Survival_df$Mutation_rate[Survival_df$type==j] = cancer.scaled4
  } 
}

Survival_df = Survival_df[-which(is.na(Survival_df$Aneuploidy_score)),]

Survival_df = Survival_df[-which(is.na(Survival_df$Mutation_rate)),]

if(subset == "all"){
  Survival_df = Survival_df
}else{Survival_df = Survival_df[which(Survival_df$ED == subset),]}

if(exclude_medium == "exclude_medium"){
  Survival_df = Survival_df[-which(Survival_df$ICR_cluster == "ICR Medium"),]
}

if(!is.na(Cancer)){
  Survival_df = Survival_df[which(Survival_df$type %in% Cancer),]
}

Survival_df$ajcc_pathologic_tumor_stage[which(Survival_df$ajcc_pathologic_tumor_stage %in% c("[Discrepancy]", "[Not Applicable]", "[Not Available]", "[Unknown]"))] = NA
Survival_df$ajcc_pathologic_tumor_stage[which(Survival_df$ajcc_pathologic_tumor_stage %in% c("Stage I", "Stage IA", "Stage IB"))] = "Stage I"
Survival_df$ajcc_pathologic_tumor_stage[which(Survival_df$ajcc_pathologic_tumor_stage %in% c("Stage II", "Stage IIA", "Stage IIB", "Stage IIC"))] = "Stage II"
Survival_df$ajcc_pathologic_tumor_stage[which(Survival_df$ajcc_pathologic_tumor_stage %in% c("Stage III", "Stage IIIA", "Stage IIIB", "Stage IIIC"))] = "Stage III"
Survival_df$ajcc_pathologic_tumor_stage[which(Survival_df$ajcc_pathologic_tumor_stage %in% c("Stage IV", "Stage IVA", "Stage IVB", "Stage IVC"))] = "Stage IV"
Survival_df$ajcc_pathologic_tumor_stage = factor(Survival_df$ajcc_pathologic_tumor_stage, levels = c("I/II NOS", "IS", "Stage 0", "Stage I", "Stage II", "Stage III", "Stage IV", "Stage X"))

Survival_df$histological_grade[which(Survival_df$histological_grade %in% c("[Discrepancy]", "[Not Available]", "[Unknown]", "GX", "High Grade", "Low Grade", "GB"))] = NA
Survival_df$histological_grade[which(Survival_df$histological_grade %in% c("G1"))] = "Grade 1"
Survival_df$histological_grade[which(Survival_df$histological_grade %in% c("G2"))] = "Grade 2"
Survival_df$histological_grade[which(Survival_df$histological_grade %in% c("G3"))] = "Grade 3"
Survival_df$histological_grade[which(Survival_df$histological_grade %in% c("G4"))] = "Grade 4"
Survival_df$histological_grade[which(Survival_df$type %in% c("GBM"))] = "Grade 4"

Y = Surv_cutoff_years * 365
TS.Alive = Survival_df[Survival_df[, Outcome] == "0", c(Outcome,  paste0(Outcome, ".time"), "ICR_cluster", "ajcc_pathologic_tumor_stage", 
                                                        "histological_grade", "type",
                                                        "TGF_beta_ES", "Proliferation_ES",
                                                        "Aneuploidy_score", "Mutation_rate")]
colnames(TS.Alive) = c("Status","Time", "ICR_cluster", "pathologic_stage", 
                       "Grade", "Cancer",
                       "TGF-beta", "Proliferation",
                       "Aneuploidy_score", "Mutation_rate")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Survival_df[Survival_df[, Outcome] == "1", c(Outcome,  paste0(Outcome, ".time"), "ICR_cluster", "ajcc_pathologic_tumor_stage", 
                                                       "histological_grade", "type",
                                                       "TGF_beta_ES", "Proliferation_ES",
                                                       "Aneuploidy_score", "Mutation_rate")]
colnames(TS.Dead) = c("Status","Time", "ICR_cluster", "pathologic_stage", 
                      "Grade", "Cancer",
                      "TGF-beta", "Proliferation",
                      "Aneuploidy_score", "Mutation_rate")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "0"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "1"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                                                                                         # remove patients with less then 1 day follow up time

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
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage I")] = "Stage I"
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage II")] = "Stage II"
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage III")] = "Stage III"
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage IV")] = "Stage IV"
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage X")] = NA

TS.Surv$Grade_Stage = TS.Surv$pathologic_stage
TS.Surv$Grade_Stage[which(TS.Surv$Cancer %in% c("GBM", "LGG"))] = TS.Surv$Grade[which(TS.Surv$Cancer %in% c("GBM", "LGG"))]

TS.Surv$Grade_Stage[which(TS.Surv$Grade_Stage %in% c("Grade 1", "Stage I"))] = "Stage/Grade 1"
TS.Surv$Grade_Stage[which(TS.Surv$Grade_Stage %in% c("Grade 2", "Stage II"))] = "Stage/Grade 2"
TS.Surv$Grade_Stage[which(TS.Surv$Grade_Stage %in% c("Grade 3", "Stage III"))] = "Stage/Grade 3"
TS.Surv$Grade_Stage[which(TS.Surv$Grade_Stage %in% c("Grade 4", "Stage IV"))] = "Stage/Grade 4"

#TS.Surv$Grade_Stage = factor(TS.Surv$Grade_Stage, levels = c("Stage/Grade 1", "Stage/Grade 2", "Stage/Grade 3", "Stage/Grade 4"))

TS.Surv$Grade_Stage[which(TS.Surv$Grade_Stage %in% c("Stage/Grade 1"))] = 1
TS.Surv$Grade_Stage[which(TS.Surv$Grade_Stage %in% c("Stage/Grade 2"))] = 2
TS.Surv$Grade_Stage[which(TS.Surv$Grade_Stage %in% c("Stage/Grade 3"))] = 3
TS.Surv$Grade_Stage[which(TS.Surv$Grade_Stage %in% c("Stage/Grade 4"))] = 4
TS.Surv$Grade_Stage = as.numeric(TS.Surv$Grade_Stage)

TS.Surv$ICR_cluster = factor(TS.Surv$ICR_cluster, levels = c("ICR High", "ICR Low")) # adjust this by hard coding when "ICR Medium" is included



#TS.Surv = TS.Surv[-which(TS.Surv$Cancer %in% c("LGG")),]
#TS.Surv = TS.Surv[-which(is.na(TS.Surv$Grade_Stage)),]
#TS.Surv = TS.Surv[which(TS.Surv$Cancer %in% ICR_neutral_cancers),]

TS.Surv$`TGF-beta` = (TS.Surv$`TGF-beta` - min(TS.Surv$`TGF-beta`))/(max(TS.Surv$`TGF-beta`)-min(TS.Surv$`TGF-beta`))

TS.Surv$Proliferation <- (TS.Surv$Proliferation-min(TS.Surv$Proliferation))/(max(TS.Surv$Proliferation)
                                                                             -min(TS.Surv$Proliferation))

dir.create("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/9_Figures", showWarnings = FALSE)


#Uni-variate
#1
uni_variate_ICR = coxph(formula = Surv(Time, Status) ~ ICR_cluster, data = TS.Surv)
cox.zph(uni_variate_ICR)
png(paste0("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/9_Figures/",
           subset, "_ICR_univariate.png"), width = 8, height = 4, units = "in", res = 600)
plot(cox.zph(uni_variate_ICR))
dev.off()
summary(uni_variate_ICR)

#2
uni_variate_ps = coxph(formula = Surv(Time, Status) ~ Grade_Stage, data = TS.Surv)
png(paste0("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/9_Figures/",
           subset, "_grade_stage_univariate.png"), width = 8, height = 4, units = "in", res = 600)
plot(cox.zph(uni_variate_ps))
dev.off()
#summary(uni_variate_ps)

#3
uni_variate_mut = coxph(formula = Surv(Time, Status) ~ Mutation_rate, data = TS.Surv)
cox.zph(uni_variate_mut)
png(paste0("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/9_Figures/", MC, "_",
           subset, "_mut_univariate.png"), width = 8, height = 4, units = "in", res = 600)
plot(cox.zph(uni_variate_mut))
dev.off()
summary(uni_variate_mut)

#4
uni_variate_aneup = coxph(formula = Surv(Time, Status) ~ Aneuploidy_score, data = TS.Surv)
cox.zph(uni_variate_aneup)
png(paste0("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/9_Figures/", MC, "_",
           subset, "_aneupl_univariate.png"), width = 8, height = 4, units = "in", res = 600)
plot(cox.zph(uni_variate_aneup))
dev.off()
summary(uni_variate_aneup)

#5 
uni_variate_prolif = coxph(formula = Surv(Time, Status) ~ Proliferation, data = TS.Surv)
cox.zph(uni_variate_prolif)
png(paste0("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/9_Figures/", MC, "_",
           subset, "_proliferation_univariate.png"), width = 8, height = 4, units = "in", res = 600)
plot(cox.zph(uni_variate_prolif))
dev.off()
summary(uni_variate_prolif)

#6
uni_variate_TGF_beta = coxph(formula = Surv(Time, Status) ~ `TGF-beta`, data = TS.Surv)
cox.zph(uni_variate_TGF_beta)
png(paste0("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/9_Figures/", MC, "_",
           subset, "_TGF_beta_univariate.png"), width = 8, height = 4, units = "in", res = 600)
plot(cox.zph(uni_variate_TGF_beta))
dev.off()

summary(uni_variate_TGF_beta)

#123456
multivariate_rev = coxph(formula = Surv(Time, Status) ~ ICR_cluster + Mutation_rate + Aneuploidy_score + Proliferation + `TGF-beta` + strata(Grade_Stage), data = TS.Surv)
cox.zph(multivariate_rev)
#write.csv(as.data.frame(print(cox.zph(multivariate_rev))), file = "~/Desktop/test.csv")
#plot(cox.zph(multivariate_rev))
sink(paste0("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision2/9_Figures/", "v1_", subset, "_", MC,"_multivar_strata_Stage_Grade.txt"))
print(summary(multivariate_rev))
sink()

#multivariate_rev = coxph(formula = Surv(Time, Status) ~ ICR_cluster + Proliferation + Mutation_rate + `TGF-beta` + Aneuploidy_score + Grade_Stage, data = TS.Surv)
#summary(multivariate_rev)
