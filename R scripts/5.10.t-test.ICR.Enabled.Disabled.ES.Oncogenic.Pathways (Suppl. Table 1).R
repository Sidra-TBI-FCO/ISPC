#################################################################
###
###
##
#################################################################

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))
source(paste0(code_path,"R tools/heatmap.3.R"))

required.packages <- c("gtools", "circlize")
ibiopak("ComplexHeatmap")
ipak(required.packages)

# Set Parameters
Pathways = "ALL"
CancerTYPES = "ALL"
Pathway_skip = ""                                                                                                        
download.method = "Assembler_Panca_Normalized_filtered"                                                                          # Specify download method (this information to be used when saving the file)
assay.platform = "gene_RNAseq"                                                                                          # Specify to which location TCGA-Assembler_v2.0.3 was downloaded
Legend = c("ICR Low","ICR Med","ICR High")
Gene.set = "Selected.pathways"
Source_surv_data = "Cell_paper"
ICR_k = "HML_classification"
Cutoff_HR = 1

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations) 

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/ssGSEA_", Gene.set,"_ES.Rdata"))
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers_Skipped_DLBC.Rdata"))

if (CancerTYPES == "ALL") { 
  CancerTYPES <- TCGA.cancersets$cancerType
}

load(paste0("./4_Analysis/",download.method,"/Pan_Cancer/Survival_Analysis/", Source_surv_data, "_Survival_analysis_High_vs_Low_Groups", 
            ICR_k, ".Rdata"))
ICR_enabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR > Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$HR <= Cutoff_HR & All_survival_analysis_data$p_value < 0.1)])
ICR_disabled_cancers = ICR_disabled_cancers[-which(ICR_disabled_cancers == "LAML")]
ICR_neutral_cancers = as.character(All_survival_analysis_data$Cancertype[which(All_survival_analysis_data$p_value >= 0.1)])

ICR_enabled_samples = rownames(ICR_cluster_assignment_allcancers)[which(ICR_cluster_assignment_allcancers$Cancer %in% ICR_enabled_cancers)]
ICR_disabled_samples = rownames(ICR_cluster_assignment_allcancers)[which(ICR_cluster_assignment_allcancers$Cancer %in% ICR_disabled_cancers)]

annotation = data.frame(sample = c(ICR_enabled_samples, ICR_disabled_samples), ED = NA)
annotation$ED[which(annotation$sample %in% ICR_enabled_samples)] = "ICR_enabled"
annotation$ED[which(annotation$sample %in% ICR_disabled_samples)] = "ICR_disabled"

ES.subset = ES.all[,which(colnames(ES.all) %in% annotation$sample)]    # subset for only ICR enabled and ICR disabled cancer types
ES.subset = t(ES.subset)                             

t.tests <- data.frame (pathway=colnames(ES.subset),p.value=0,p.value.fdr=0,mean.disabled=0,mean.enabled=0, percent_increase = 0,stringsAsFactors = FALSE)
t.tests[t.tests == 0] = NA

i=1
for (i in 1:ncol(ES.subset)){
  subset = annotation
  pathway = colnames(ES.subset)[i]
  subset$expression = ES.subset[, pathway][match(subset$sample, rownames(ES.subset))]
  p.value = t.test(expression~ED, data = subset, paired = FALSE, var.equal = FALSE)$p.value #TRUE or sign.var
  p.value.fdr = p.adjust(p = p.value, method = "fdr", n = nrow(ES.subset))                 # n is number of comparisons
  means = aggregate(subset$expression, list(subset$ED), mean)
  mean.disabled = means[means$Group.1=="ICR_disabled", 2]
  mean.enabled = means[means$Group.1=="ICR_enabled", 2]
  percent_increase = (mean.disabled-mean.enabled)/mean.enabled *100
  t.tests[t.tests$pathway==pathway,c(2:ncol(t.tests))] <- c(p.value,p.value.fdr,mean.disabled,mean.enabled, percent_increase)
}

t.tests = t.tests[order(t.tests$percent_increase),]
t.tests$Upregulated.in = t.tests$mean.enabled<t.tests$mean.disabled
t.tests$Upregulated.in[which(t.tests$Upregulated.in == TRUE)] = "ICR_disabled"
t.tests$Upregulated.in[which(t.tests$Upregulated.in == FALSE)] = "ICR_enabled"

t.tests$contribution = NA
t.tests$contribution[which(t.tests$p.value.fdr < 0.05 & t.tests$Upregulated.in == "ICR_enabled")] = "enabling"
t.tests$contribution[which(t.tests$p.value.fdr < 0.05 & t.tests$Upregulated.in == "ICR_disabled")] = "disabling"

dir.create(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/T.test.ICR.enabled.disabled.cancertypes"), 
           showWarnings = FALSE)
save(t.tests, file = paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/T.test.ICR.enabled.disabled.cancertypes/",
                            "t_tests_ES_Enabled_Disabled_Cancertypes.Rdata"))
write.csv(t.tests,file = paste0("./4_Analysis/", download.method, "/Pan_Cancer/Signature_Enrichment/T.test.ICR.enabled.disabled.cancertypes/",
                                "t_tests_ES_Enabled_Disabled_Cancertypes.csv"))




#subset = annotation
#subset$expression = ES.subset[i,]
#p.value = t.test(expression~ED, data = subset, paired = FALSE, var.equal = FALSE)$p.value #TRUE or sign.var
#p.value.fdr = p.adjust(p = p.value, method = "fdr", n = nrow(ES.subset))                 # n is number of comparisons
#pathway = rownames(ES.subset)[i]
#means = aggregate(subset$expression, list(subset$ED), mean)
#sd = aggregate(subset$expression, list(subset$ED), sd)
#mean.disabled = means[means$Group.1=="ICR_disabled", 2]
#mean.enabled = means[means$Group.1=="ICR_enabled", 2]
#sd.disabled = sd[sd$Group.1=="ICR_disabled", 2]
#sd.enabled = sd[sd$Group.1== "ICR_enabled", 2]
#ratio = mean.disabled/mean.enabled
#t.tests[t.tests$pathway==pathway,c(2:ncol(t.tests))] <- c(p.value,p.value.fdr,ratio,mean.disabled,mean.enabled,sd.disabled,sd.enabled)



