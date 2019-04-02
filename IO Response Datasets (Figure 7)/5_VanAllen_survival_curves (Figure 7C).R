

# Setup environment
rm(list=ls())
setwd("~/Dropbox (TBI-Lab)/DB-LAB/Projects - Data/02 Active Projects/SDR400074_PIFR_2018_WH_TCGA_pancancer_immuneprofiling_in-silico/IOResponse_Datasets")
source("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/R scripts/ipak.function.R")
source("./R tools/ggkm_adapted_for_paper.R")

required.packages = c("survival", "plyr", "raster", "texreg", "stringr")
ipak(required.packages)

# Set Parameters
Surv.cutoff.years = 10                                        # SET cut-off
Pathway = "TGF BETA SIGNALING"  # "TGF BETA SIGNALING" "Proliferation"
Subset = "High" # all High or Low
ICR_readout = "ICRscore"
exclude_medium = "exclude_medium"
#Group.of.interest = "ICR_HML"                                # "tumour_anatomic_site" or "ICR_cluster_k3"

# Create folders and log file
dir.create("./Figures/vanAllen",showWarnings = FALSE)
dir.create("./Figures/vanAllen/Kaplan Meier Plots", showWarnings = FALSE)
dir.create("./Analysis/vanAllen/Survival Analysis", showWarnings = FALSE)

# Load data
load("./Data/vanAllen/vanAllen.Rdata")
load(paste0("./Analysis/vanAllen/Enrichment_scores/vanAllen_Selected_Pathways_ES.Rdata"))
load(paste0("./Analysis/vanAllen/Clustering/vanAllen.ICR.reps5000/vanAllen_ICR_cluster_assignment_k2-6.Rdata"))
load(paste0("./Analysis/vanAllen/Tertiles/group2_", Pathway, "_", ICR_readout, "_tertiles_dataframe.Rdata"))

Survival_data = Clinical.data[,c("overall_survival", "progression_free", "dead", "progression")]
Survival_data$ICR_tertile = tertiles.df$ICR_tertiles[match(rownames(Survival_data), tertiles.df$Sample)]
Survival_data$Pathway = tertiles.df$Pathway_median[match(rownames(Survival_data), tertiles.df$Sample)]

if(exclude_medium == "exclude_medium"){
  Survival_data = Survival_data[-which(Survival_data$ICR_tertile == "Medium"),]
}

Survival_data$ICR_tertile = factor(Survival_data$ICR_tertile, levels = c("High", "Low"))
levels(Survival_data$ICR_tertile) = c("ICR High", "ICR Low")



# Subset to filter for proliferation
if(Subset != "all"){
  Survival_data = Survival_data[which(Survival_data$Pathway == paste0(Pathway, " ", Subset)),]
}

# time / event object creation
Y = Surv.cutoff.years * 365
TS.Alive = Survival_data[Survival_data$dead == "0", c("dead", "overall_survival", "ICR_tertile")]
colnames(TS.Alive) = c("Status","Time", "ICR_tertile")
TS.Alive$Time = as.numeric(as.character(TS.Alive$Time))
TS.Alive$Time[TS.Alive$Time > Y] = Y

TS.Dead = Survival_data[Survival_data$dead == "1", c("dead", "overall_survival", "ICR_tertile")]
colnames(TS.Dead) = c("Status","Time", "ICR_tertile")
TS.Dead$Time = as.numeric(as.character(TS.Dead$Time))
TS.Dead$Status[which(TS.Dead$Time> Y)] = "0"
TS.Dead$Time[TS.Dead$Time > Y] = Y

TS.Surv = rbind (TS.Dead,TS.Alive)
TS.Surv$Time = as.numeric(as.character(TS.Surv$Time))
TS.Surv$Status <- TS.Surv$Status == "1"
TS.Surv = subset(TS.Surv,TS.Surv$Time > 1)                               # remove patients with less then 1 day follow up time

# survival curve
msurv = Surv(TS.Surv$Time/30.4, TS.Surv$Status)                          # calculate the number of months
mfit = survfit(msurv~TS.Surv$ICR_tertile,conf.type = "log-log")

# Calculations (Needs manual adaptation!)
mdiff = survdiff(eval(mfit$call$formula), data = eval(mfit$call$data))
pval = pchisq(mdiff$chisq,length(mdiff$n) - 1,lower.tail = FALSE)
pvaltxt = ifelse(pval < 0.0001,"p < 0.0001",paste("p =", signif(pval, 3)))

Cal.Surv = TS.Surv
#Cal.Surv[,Group.of.interest] = as.factor(Cal.Surv[,Group.of.interest])
#Cal.Surv[,Group.of.interest] = relevel(Cal.Surv[,Group.of.interest], "ICR3")
Cal.Surv$ICR_tertile = as.factor(Cal.Surv$ICR_tertile)
mHR = coxph(formula = msurv ~ Cal.Surv$ICR_tertile,data = Cal.Surv, subset = Cal.Surv$ICR_tertile %in% c("ICR High", "ICR Low"))
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


# plots
png(paste0("./Figures/vanAllen/Kaplan Meier Plots/Paper_KM_ICR_tertiles_", exclude_medium, "_", Pathway, "_", ICR_readout, "_", Subset,".png"),res=600,height=3,width=4,unit="in")  # set filename
ggkm(mfit,
     timeby=12,
     ystratalabs = levels(TS.Surv$ICR_tertile),
     ystrataname = NULL,
     main= NULL,
     xlabs = "Time in months",
     cbPalette = cbPalette)
     #PLOT_P = signif(p[2],3),
     #PLOT_HR = round(signif(exp(mHR.extract@coef),3)[2], 3),
     #PLOT_CI1 = CI[2,1],
     #PLOT_CI2 = CI[2,2])
dev.off()
