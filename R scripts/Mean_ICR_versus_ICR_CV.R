

# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located

source(paste0(code_path, "R tools/ipak.function.R"))

required.packages = c("raster", "ggplot2", "ggrepel")
ipak(required.packages)

#Set parameters
ICR_disabled_cancers = c("UVM", "LGG", "PAAD", "KIRC")
ICR_enabled_cancers = c("BRCA", "SKCM", "UCEC", "SARC", "LIHC", "HNSC", "STAD", "BLCA")

# Load data
results_cv = read.csv("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision/New_results/Variation/cv_results.csv",
                      stringsAsFactors = FALSE)
results_cv$ID = "ICR_neutral"
results_cv$ID[which(results_cv$Cancer %in% ICR_disabled_cancers)] = "ICR_disabled"
results_cv$ID[which(results_cv$Cancer %in% ICR_enabled_cancers)] = "ICR_enabled"

results_cv$ID = factor(results_cv$ID, levels = c("ICR_enabled", "ICR_neutral", "ICR_disabled"))
results_cv$Cancer[which(results_cv$ID == "ICR_neutral")] = ""

plot = ggplot(results_cv, aes(x= Mean_ICR, y = Coefficient_of_variation, label = Cancer)) +
  geom_point(aes(color = results_cv$ID)) +
  theme_bw() +
  theme(legend.position = "none",
        axis.title.x = element_text(color = "black", size = 16),
        axis.title.y = element_text(color = "black", size = 16),
        axis.text.x = element_text(color = "black", size = 16),
        axis.text.y = element_text(color = "black", size = 16)) +
  geom_text_repel() +
  scale_color_manual(values = c("orange", "grey", "purple")) +
  xlab("Mean ICR") +
  ylab("Coefficient of Variation") +
  ylim(10, 30) +
  xlim(3.5, 9)

png("~/Dropbox (TBI-Lab)/DB-LAB/Manuscripts/Pancancer Jessica/For Submission JITC/Revision/New_results/Variation/Mean_ICR_vs_CV.png", res = 600,
    width = 4, height = 3.5, units = "in")
plot(plot)
dev.off()
