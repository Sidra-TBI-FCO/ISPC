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

required.packages <- c("ggplot2", "Hmisc")
ipak(required.packages)

## Set Parameters
CancerTYPES = c("ALL")                                                                                                   # Specify the cancertypes that you want to download or process, c("...","...") or "ALL"
Cancer_skip = c("")                                                                                                      # If CancerTYPES = "ALL", specify cancertypes that you do not want to download
download.method = "Assembler_Panca_Normalized_filtered"                                                                                       # Specify download method (this information to be used when saving the file)
variable_of_interest = "SNV.Neoantigens" #SNV.Neoantigens #Nonsilent.Mutation.Rate
log_transform = "log_transform" #"not_transformed" #"log_transform"

# Load data
TCGA.cancersets = read.csv(paste0(code_path, "Datalists/TCGA.datasets.csv"),stringsAsFactors = FALSE)                    # TCGA.datasets.csv is created from Table 1. (Cancer Types Abbreviations)
PanImmune_MS = read.csv("./3_DataProcessing/External/mmc2-PanImmune_MS.csv", stringsAsFactors = FALSE)
load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Clustering/ICR_cluster_assignment_allcancers_Skipped_DLBC.Rdata"))

df_plot = data.frame(sample_barcodes = rownames(ICR_cluster_assignment_allcancers), ICR_cluster = ICR_cluster_assignment_allcancers$HML_cluster,
                     Cancer = ICR_cluster_assignment_allcancers$Cancer,
                     variable = PanImmune_MS[,variable_of_interest][match(substring(rownames(ICR_cluster_assignment_allcancers), 1, 12),
                                                                   PanImmune_MS$TCGA.Participant.Barcode)],
                     stringsAsFactors = FALSE)

df_plot = df_plot[-which(is.na(df_plot$variable)),]

df_plot$variable = as.numeric(df_plot$variable)

dir.create("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Mutational_load", showWarnings = FALSE)
save(df_plot, file = paste0("./4_Analysis/", download.method, "/Pan_Cancer/Mutational_load/", variable_of_interest, ".Rdata"))

if(log_transform == "log_transform"){
  df_plot$variable = log10(df_plot$variable + 1)
}
## For log transformation of y-axis if not transformed change 0 to 0.0001
df_plot$variable[which(df_plot$variable == 0)] = 0.0001

df_plot_agg = aggregate(df_plot$variable, by = list(df_plot$Cancer), FUN = mean)
df_plot_agg$Group.1 = as.character(df_plot_agg$Group.1)
cancer_order = df_plot_agg$Group.1[order(df_plot_agg$x)]
save(cancer_order, file = paste0("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Mutational_load/cancer_order_", variable_of_interest))

if(variable_of_interest == "SNV.Neoantigens"){
  load(paste0("./4_Analysis/Assembler_Panca_Normalized_filtered/Pan_Cancer/Mutational_load/cancer_order_Nonsilent.Mutation.Rate"))
}

df_plot$Cancer = factor(df_plot$Cancer, levels = cancer_order)

plot = ggplot(df_plot, aes(x=reorder(Cancer, variable, FUN = mean), y = variable)) +
  geom_jitter(width = 0.2, size = 0.3) +
  stat_summary(fun.y =  mean, fun.ymin = mean,
               fun.ymax = mean,geom = "crossbar",
               width = 0.5, color = "red")

load(paste0("./4_Analysis/", download.method, "/Pan_Cancer/Mutational_load/t.test.table_", variable_of_interest, "_", log_transform,".Rdata"))
t.tests$sig = NA
t.tests$sig[which(t.tests$p.value < 0.05)] = "*"
t.tests$Cancer = factor(t.tests$Cancer, levels = cancer_order)

dir.create(paste0("./5_Figures/Pancancer_plots/", download.method, "/", variable_of_interest, "_scatterplots"), showWarnings = FALSE)

if(variable_of_interest == "Nonsilent.Mutation.Rate"){
  y = 430
}
if(variable_of_interest == "SNV.Neoantigens"){
  y = 3.7
}

if(log_transform == "log_transform"){
  axis_title = paste0(gsub("\\.", " ", variable_of_interest), "\n log10(Value + 1)")
}
if(log_transform == "not_transformed" & variable_of_interest == "Nonsilent.Mutation.Rate"){
  axis_title = paste0("Nonsilent Mutations \nper Megabase")
  }

png(paste0("./5_Figures/Pancancer_plots/", download.method, "/", variable_of_interest, "_scatterplots/final_", variable_of_interest, 
           "_", log_transform, "_scatterplot.png"), res = 600, width = 17, height = 5, units = "in")

plot = ggplot(df_plot, aes(x=ICR_cluster, y = variable)) +
  geom_jitter(aes(color = ICR_cluster), width = 0.02, size = 0.3) +
  theme_bw() +
  theme(axis.text.x=element_blank()) +
  labs(title = variable_of_interest,x = "Cancer types", y = axis_title) +
  scale_color_manual(values = c("blue", "green", "red")) +
  stat_summary(fun.y =  mean, fun.ymin = mean,
               fun.ymax = mean,geom = "crossbar",
               width = 0.5, color = "red")

plot = plot + facet_grid(. ~ Cancer, switch = "x")

if(log_transform == "not_transformed"){
plot = plot + scale_y_continuous(breaks = c(0,1,10,100,1000)) +
    coord_trans(y = "log10", limy=c(0.01,1000)) 
}

# The difference between transforming the scales and
# transforming the coordinate system is that scale
# transformation occurs BEFORE statistics, and coordinate
# transformation afterwards. https://ggplot2.tidyverse.org/reference/coord_trans.html
# In this case we need to use coordinate transformation and no transformation of the scales as below:
  #scale_y_log10(breaks = scales::trans_breaks("log10", function(x) 10^x),
  #labels = scales::trans_format("log10", scales::math_format(10^.x)),
  #limits = c(0.01,500))
  

plot = plot + 
  geom_text(data = t.tests,
            mapping = aes(x= 1.5, y = y, label = sig),
            hjust = -0.1,
            vjust = -0.1,
            size = 7) +
  theme(
    strip.background = element_blank(),
    text = element_text(size = 20),
    strip.text = element_text(angle = 90),
    axis.ticks.x = element_blank()
  )

plot(plot)
dev.off()


