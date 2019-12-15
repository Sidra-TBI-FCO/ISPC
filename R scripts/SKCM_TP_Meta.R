
## Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("survival","reshape","ggplot2","plyr","Rcpp","colorspace","texreg", "mi", "dplyr", "ggpubr")
required.bioconductor.packages = "survival"
ipak(required.packages)
ibiopak(required.bioconductor.packages)

source("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/R scripts/ggkm_Jessica_Pancancer.R")

# Change this parameter!
subset = "all"                              #"ICR_enabled", "ICR_disabled", "ICR_neutral", or "all"

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
Cancer = "SKCM"  #c("LGG", "GBM")     # NA or "KIRC"
Tissue_type = "Metastatic" #"Primary Solid Tumor" or "Metastatic"

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
Survival_df$tcga_sample_barcode = rownames(Hallmark_and_ICR_cluster_assignment_allcancers)[match(Survival_df$bcr_patient_barcode, substring(rownames(Hallmark_and_ICR_cluster_assignment_allcancers), 1, 12))]
Survival_df$ED = NA
Survival_df$Mutation_rate = PanImmune_MS$Nonsilent.Mutation.Rate[match(Survival_df$bcr_patient_barcode, PanImmune_MS$TCGA.Participant.Barcode)]
Survival_df$Mutation_rate = log10(Survival_df$Mutation_rate + 0.0001)
Survival_df$ICRscore = Hallmark_and_ICR_cluster_assignment_allcancers$ICRscore[match(Survival_df$bcr_patient_barcode, substring(rownames(Hallmark_and_ICR_cluster_assignment_allcancers), 1, 12))]

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

#Survival_df = Survival_df[-which(is.na(Survival_df$Aneuploidy_score)),]

#Survival_df = Survival_df[-which(is.na(Survival_df$Mutation_rate)),]

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
TS.Alive = Survival_df[Survival_df[, Outcome] == "0", c(Outcome,  paste0(Outcome, ".time"), "tcga_sample_barcode", "ICR_cluster", "ajcc_pathologic_tumor_stage", 
                                                        "histological_grade", "type",
                                                        "TGF_beta_ES", "Proliferation_ES",
                                                        "Aneuploidy_score", "Mutation_rate", "ICRscore")]
colnames(TS.Alive) = c("Status","Time", "tcga_sample_barcode", "ICR_cluster", "pathologic_stage", 
                       "Grade", "Cancer",
                       "TGF-beta", "Proliferation",
                       "Aneuploidy_score", "Mutation_rate", "ICRscore")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Survival_df[Survival_df[, Outcome] == "1", c(Outcome,  paste0(Outcome, ".time"), "tcga_sample_barcode", "ICR_cluster", "ajcc_pathologic_tumor_stage", 
                                                       "histological_grade", "type",
                                                       "TGF_beta_ES", "Proliferation_ES",
                                                       "Aneuploidy_score", "Mutation_rate", "ICRscore")]
colnames(TS.Dead) = c("Status","Time", "tcga_sample_barcode", "ICR_cluster", "pathologic_stage", 
                      "Grade", "Cancer",
                      "TGF-beta", "Proliferation",
                      "Aneuploidy_score", "Mutation_rate", "ICRscore")
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
#TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage I")] = "Stage I&II"
#TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage II")] = "Stage I&II"
#TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage III")] = "Stage III&IV"
#TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage IV")] = "Stage III&IV"
TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage == "Stage X")] = NA

#TS.Surv$pathologic_stage = factor(TS.Surv$pathologic_stage, levels = c("Stage I&II", "Stage III&IV"))

TS.Surv$Grade_Stage = TS.Surv$pathologic_stage
TS.Surv$Grade_Stage[which(TS.Surv$Cancer %in% c("GBM", "LGG"))] = TS.Surv$Grade[which(TS.Surv$Cancer %in% c("GBM", "LGG"))]
TS.Surv$Grade_Stage[which(TS.Surv$Grade_Stage %in% c("Grade 1", "Stage I"))] = "Stage/Grade 1"
TS.Surv$Grade_Stage[which(TS.Surv$Grade_Stage %in% c("Grade 2", "Stage II"))] = "Stage/Grade 2"
TS.Surv$Grade_Stage[which(TS.Surv$Grade_Stage %in% c("Grade 3", "Stage III"))] = "Stage/Grade 3"
TS.Surv$Grade_Stage[which(TS.Surv$Grade_Stage %in% c("Grade 4", "Stage IV"))] = "Stage/Grade 4"

TS.Surv$Grade_Stage = factor(TS.Surv$Grade_Stage, levels = c("Stage/Grade 1", "Stage/Grade 2", "Stage/Grade 3", "Stage/Grade 4"))

TS.Surv$ICR_cluster = factor(TS.Surv$ICR_cluster, levels = c("ICR High","ICR Low")) # adjust this by hard coding when "ICR Medium" is included


#TS.Surv = TS.Surv[-which(TS.Surv$Cancer %in% c("LGG")),]
#TS.Surv = TS.Surv[-which(is.na(TS.Surv$Grade_Stage)),]
#TS.Surv = TS.Surv[which(TS.Surv$Cancer %in% ICR_neutral_cancers),]

#TS.Surv$`TGF-beta` = (TS.Surv$`TGF-beta` - min(TS.Surv$`TGF-beta`))/(max(TS.Surv$`TGF-beta`)-min(TS.Surv$`TGF-beta`))
#
#TS.Surv$Proliferation <- (TS.Surv$Proliferation-min(TS.Surv$Proliferation))/(max(TS.Surv$Proliferation)
#                                                                             -min(TS.Surv$Proliferation))

TS.Surv$Tissue_type = substring(TS.Surv$tcga_sample_barcode, 14,15)
TS.Surv$Tissue_type[which(TS.Surv$Tissue_type == "01")] = "Primary Solid Tumor"
TS.Surv$Tissue_type[which(TS.Surv$Tissue_type == "06")] = "Metastatic"
TS.Surv$Tissue_type = factor(TS.Surv$Tissue_type, levels = c("Primary Solid Tumor", "Metastatic"))

if(Tissue_type == "all"){}else{
  TS.Surv = TS.Surv[which(TS.Surv$Tissue_type == Tissue_type),]
}


# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                                                                                    # calculate the number of months
mfit = survfit(msurv~TS.Surv$ICR_cluster,conf.type = "log-log")

# Calculations
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

# Check this!!
##TS.Surv[,"Group"] = relevel(TS.Surv[,"Group"], "ICR High")
mHR = coxph(formula = msurv ~ TS.Surv$ICR_cluster, subset = TS.Surv$ICR_cluster %in% c("ICR High", "ICR Medium", "ICR Low"))
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

PLOT_P = signif(p[1],3)
PLOT_HR = round(signif(exp(mHR.extract@coef),3)[1], 3)
PLOT_CI1 = CI[1,1]
PLOT_CI2 = CI[1,2]

uni_variate_ICRscore = coxph(formula = Surv(Time, Status) ~ ICRscore, data = TS.Surv)
summary(uni_variate_ICRscore)
dir.create(paste0("./5_Figures/Pancancer_plots/Assembler_Panca_Normalized_filtered/Survival_Plots/Benefit_clusters"), showWarnings = FALSE)

#png(paste0("./5_Figures/Pancancer_plots/Assembler_Panca_Normalized/Survival_Plots/Test/",
#     "Kaplan_Meier_", Stagenames, "_", subset, "_", Cancer, "_samples.png"),
#res=600,height=6,width=8,unit="in")                                                                                           # set filename
dev.new()
ggkm(mfit,
     timeby=12,
     ystratalabs = c("ICR High", "ICR Low"),
     ystrataname = NULL,
     xlabs = "Time in months",
     palette = c("red","blue"),
     legend = "none")
     #PLOT_HR = PLOT_HR,
     #PLOT_P = PLOT_P,
     #PLOT_CI1 = PLOT_CI1,
     #PLOT_CI2 = PLOT_CI2)
dev.off()

#TS.Surv$pathologic_stage = as.character(TS.Surv$pathologic_stage)
#TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage %in% c("Stage I"))] = 1
#TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage %in% c("Stage II"))] = 2
#TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage %in% c("Stage III"))] = 3
#TS.Surv$pathologic_stage[which(TS.Surv$pathologic_stage %in% c("Stage IV"))] = 4
#TS.Surv$pathologic_stage = as.numeric(TS.Surv$pathologic_stage)

#uni_variate_stage = coxph(formula = Surv(Time, Status) ~ pathologic_stage, data = TS.Surv)
#cox.zph(uni_variate_stage)

#contingency_table = table(patient.table.all$ICR_cluster, patient.table.all$pathologic_stage)
#chi2 = chisq.test(contingency_table)
#chi2$p.value

table(TS.Surv$Cancer)
table(TS.Surv$ICR_cluster)
table(TS.Surv$Tissue_type)

contingency_table = table(TS.Surv$Tissue_type, TS.Surv$ICR_cluster)
chisq.test(contingency_table)

levels(TS.Surv$Tissue_type) = c("Primary", "Metastatic")

DF1 <- TS.Surv %>%
  group_by(Tissue_type, ICR_cluster) %>%
  summarise(count=n()) %>%
  mutate(perc=count/sum(count))

colors = c("ICR High" = "red", "ICR Medium"= "green", "ICR Low" = "blue")
plot = ggplot(DF1, aes(x = Tissue_type, y =perc*100, fill = ICR_cluster)) + geom_bar(stat="identity") +
  labs(x = "Tissue type", y = "Percentage", fill = "ICR_cluster", face = "bold") +
  theme(panel.grid = element_line(linetype = "solid", colour = "white"), 
        panel.background = element_rect(fill = "white", colour = "grey"),
        axis.text = element_text(size = 19, colour = "black"),
        axis.text.x = element_text(size = 19, colour = "black", angle = 45, vjust = 0.9, hjust = 0.9),
        axis.title = element_text(size = 19, colour = "black"),
        legend.title = element_blank(),
        legend.text = element_text(size = 19, colour = "black"))+
  scale_fill_manual(values= colors)
dev.new()
plot(plot)

plot = ggplot(TS.Surv, aes(x = Tissue_type, y = ICRscore)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2) +
  theme_bw() +
  theme(axis.text.x = element_text(colour = "black", size = 17, angle = 90, vjust = 0.5, hjust = 0.9),
        axis.text.y = element_text(colour = "black", size = 17),
        axis.title.y = element_text(colour = "black", size = 17)) +
  stat_compare_means(method = "t.test") +
  xlab("")

plot(plot)




