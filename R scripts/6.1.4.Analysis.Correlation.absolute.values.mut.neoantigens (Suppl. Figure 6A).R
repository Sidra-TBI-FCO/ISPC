# Setup environment
rm(list=ls())

setwd("~/Dropbox (TBI-Lab)/TCGA Analysis pipeline/")                                                                    # Setwd to location were output files have to be saved.
#setwd("~/Dropbox (TBI-Lab)/External Collaborations/TCGA Analysis pipeline/")    

code_path = "~/Dropbox (Personal)/Jessica PhD Project/QCRI-SIDRA-ICR-Jessica/"                                          # Set code path to the location were the R code is located
#code_path = "~/Dropbox (Personal)/R-projects/QCRI-SIDRA-ICR/" 
#code_path = "C:/Users/whendrickx/R/GITHUB/TCGA_Pipeline/"

source(paste0(code_path, "R tools/ipak.function.R"))
required.packages = c("ggrepel")
ipak(required.packages)

#Set parameters
download.method = "Assembler_Panca_Normalized_filtered"
log_transform = "not_transformed" #"not_transformed" or "log_transform"
#variable_of_interest = "Nonsilent.Mutation.Rate" #SNV.Neoantigens #Nonsilent.Mutation.Rate

#Load data
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Mutational_load/t.test.table_Nonsilent.Mutation.Rate_", log_transform,".Rdata"))
t.tests.mut = t.tests
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Mutational_load/t.test.table_SNV.Neoantigens_", log_transform,".Rdata"))
t.tests.neo = t.tests
rm(list = c("t.tests"))
load("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Cancer_color_table.Rdata")

plot_df = t.tests.mut[, c("Cancer", "mean.all")]
colnames(plot_df) = c("Cancer", "Mean Mutation Rate")
plot_df$`Mean Neoantigen count` = t.tests.neo$mean.all[match(plot_df$Cancer, t.tests.neo$Cancer)]

dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/Correlation_Neo_Mut_FC_Plots/"), showWarnings = FALSE)
png(paste0("./5_Figures/Pancancer_plots/", download.method,"/Correlation_Neo_Mut_FC_Plots/mean_mut_vs_mean_neo_each_cancer_log_scale.png"),
    res = 600, height = 8, width = 8, units = "in")
ggplot(plot_df, aes(x = `Mean Mutation Rate`, y = `Mean Neoantigen count`)) +
  geom_point(size = 0.5, shape = 1) +
  geom_text_repel(aes(label = Cancer), size = 4) +
  theme_bw() +
  #scale_color_manual(breaks = c(Cancer_color_table$Group.1),
  #values = c(Cancer_color_table$color)) + (doesn't work yet)
  theme(panel.grid = element_blank(),
        axis.text=element_text(size=15),
        axis.title=element_text(size=16)) +
  scale_x_log10()+
  scale_y_log10() +
  labs(x = paste0("Mean Mutation Rate"), y = paste0("Mean Neoantigen count"))
#geom_smooth(method=lm, se = FALSE, color="blue")

dev.off()
