
## Predictive impact of ICR alone

# VanAllen dataset

# Setup environment
rm(list = ls())
setwd("~/Dropbox (TBI-Lab)/DB-LAB/Projects - Data/02 Active Projects/SDR400074_PIFR_2018_WH_TCGA_pancancer_immuneprofiling_in-silico/IOResponse_Datasets")
source("~/Dropbox (TBI-Lab)/DB-LAB/bioinformatics tools/R scripts/ipak.function.R")

required.packages = c("ggplot2", "ggpubr")
ipak(required.packages)

# Set parameters
dataset = "vanAllen"
clinical_outcome = "group2" #"Response_PD_PRCR" #"Response_PD_others" #"RECIST" or "group"
extremes = "extremes2" #extremes2 "no"
ICR_readout = "ICRscore" # "ICRscore"

# Load data
load(paste0("./Analysis/", dataset, "/Enrichment_scores/", dataset, "_Selected_Pathways_ES.Rdata"))
load(paste0("./Analysis/", dataset, "/Clustering/", dataset, ".ICR.reps5000/", dataset, "_ICR_cluster_assignment_k2-6.Rdata"))
load(paste0("./Data/", dataset, "/", dataset, ".Rdata"))

# Create directories
dir.create("./Figures/vanAllen/ICR_Boxplots", showWarnings = FALSE)
dir.create(paste0("./Figures/vanAllen/ICR_Boxplots/", clinical_outcome), showWarnings = FALSE)

Clinical.data$group2 = Clinical.data$group
Clinical.data$group2 = factor(Clinical.data$group2, levels = c("nonresponse", "long-survival", "response"))
levels(Clinical.data$group2) = c("nonresponse", "long-survival or response", "long-survival or response")
Clinical.data$group3 = Clinical.data$group
Clinical.data$group3 = factor(Clinical.data$group3, levels = c("nonresponse", "long-survival", "response"))
levels(Clinical.data$group3) = c("nonresponse", "nonresponse", "response")
Clinical.data$RECIST = factor(Clinical.data$RECIST, levels = c("PD", "SD", "PR", "CR", "X"))
Clinical.data$group = factor (Clinical.data$group, levels = c("nonresponse",  "long-survival", "response"))

if(clinical_outcome == "Response_PD_others"){
  Clinical.data$Response_PD_others = Clinical.data$RECIST
  Clinical.data = Clinical.data[which(Clinical.data$Response_PD_others %in% c("PD", "SD", "PR", "CR")),]
  levels(Clinical.data$Response_PD_others) = c("PD", "SDPRCR", "SDPRCR", "SDPRCR", "X")
  ES = ES[,which(colnames(ES) %in% rownames(Clinical.data))]
  table_cluster_assignment = table_cluster_assignment[which(rownames(table_cluster_assignment) %in% rownames(Clinical.data)),]
}

if(clinical_outcome == "Response_PD_PRCR"){
  Clinical.data$Response_PD_PRCR = Clinical.data$RECIST
  Clinical.data = Clinical.data[which(Clinical.data$Response_PD_PRCR %in% c("PD", "PR", "CR")),]
  levels(Clinical.data$Response_PD_PRCR) = c("PD", "SD", "PRCR", "PRCR", "X")
  ES = ES[,which(colnames(ES) %in% rownames(Clinical.data))]
  table_cluster_assignment = table_cluster_assignment[which(rownames(table_cluster_assignment) %in% rownames(Clinical.data)),]
}
if(clinical_outcome == "group" & extremes == "extremes"){
  Clinical.data = Clinical.data[which(Clinical.data$group %in% c("nonresponse", "response")),]
  ES = ES[,which(colnames(ES) %in% rownames(Clinical.data))]
  table_cluster_assignment = table_cluster_assignment[which(rownames(table_cluster_assignment) %in% rownames(Clinical.data)),]
}

rownames(ES) = gsub(".*] ", "",rownames(ES))
rownames(ES) = gsub("./", "_", rownames(ES))

Pathways = rownames(ES)
N.pathways = length(Pathways)

i=1
for (i in 1:N.pathways){
  Pathway = Pathways[i]
  tertiles.df = data.frame(Sample = rownames(table_cluster_assignment), ICR_score = table_cluster_assignment[,ICR_readout],
                           Pathway_ES = NA, ICR_tertiles = NA, Pathway_tertiles = NA)
  tertiles.df$Pathway_ES = ES[which(rownames(ES) == Pathway),]
  
  # ICR tertiles
  first = unname(quantile(tertiles.df$ICR_score, probs = seq(0, 1, length.out = 4))[2])
  second = unname(quantile(tertiles.df$ICR_score, probs = seq(0, 1, length.out = 4))[3])
  tertiles.df$ICR_tertiles = NA
  tertiles.df$ICR_tertiles[which(tertiles.df$ICR_score < first)] = "Low"
  tertiles.df$ICR_tertiles[which(tertiles.df$ICR_score >= first & tertiles.df$ICR_score < second)] = "Medium"
  tertiles.df$ICR_tertiles[which(tertiles.df$ICR_score >= second)] = "High"
  
  # Pathway tertiles
  first = unname(quantile(tertiles.df$Pathway_ES, probs = seq(0, 1, length.out = 4))[2])
  second = unname(quantile(tertiles.df$Pathway_ES, probs = seq(0, 1, length.out = 4))[3])
  tertiles.df$Pathway_tertiles = NA
  tertiles.df$Pathway_tertiles[which(tertiles.df$Pathway_ES < first)] = "Low"
  tertiles.df$Pathway_tertiles[which(tertiles.df$Pathway_ES >= first & tertiles.df$Pathway_ES < second)] = "Medium"
  tertiles.df$Pathway_tertiles[which(tertiles.df$Pathway_ES >= second)] = "High"
  
  tertiles.df$Pathway_median = NA
  median = median(tertiles.df$Pathway_ES)
  tertiles.df$Pathway_median[which(tertiles.df$Pathway_ES < median)] = paste0(Pathway, " Low")
  tertiles.df$Pathway_median[which(tertiles.df$Pathway_ES >= median)] = paste0(Pathway, " High")
  
  # Plot df1
  plot_df = data.frame(Sample = rownames(table_cluster_assignment), ICR_score = table_cluster_assignment[,ICR_readout],
                       Response = NA)
  plot_df$Response = Clinical.data[, clinical_outcome][match(plot_df$Sample, rownames(Clinical.data))]
  
  plot = ggplot(plot_df, aes(x = Response, y = ICR_score, fill = Response), outlier.shape = NA) +
    ylim(-2, max(plot_df$ICR_score)*1.2) +
    geom_boxplot(outlier.shape = NA) + geom_jitter(position = position_jitter(width = 0.15)) + theme_bw() + 
    stat_compare_means(method = "t.test", comparisons = list(c("nonresponse", "long-survival or response")),
                       size = 6.5, label.y = max(plot_df$ICR_score * 1.1)) + # remove comparisons when clinical outcome is different from group
    #ylab(gsub("_", " ", ICR_readout)) + scale_fill_manual(values=c("#E8A025", "#099E74")) +
    ylab("ICR score") + scale_fill_manual(values=c("#E8A025", "#099E74")) +
    theme(text = element_text(size = 18, colour = "black"), axis.text.x = element_blank())
  
  
  png(filename = paste0("./Figures/vanAllen/ICR_Boxplots/", clinical_outcome, "/1_Boxplot_Overall_", ICR_readout, "_by_Response.png"), res = 600,
      width = 5.5, height = 4, units = "in")
      #width = 5, height = 5, units = "in")
  plot(plot)
  dev.off()
  
  # Plot facet by Pathway
  plot_df = data.frame(Sample = tertiles.df$Sample, Pathway = tertiles.df$Pathway_median, ICR_score = tertiles.df$ICR_score,
                       Response = NA)
  plot_df$Response = Clinical.data[, clinical_outcome][match(plot_df$Sample, rownames(Clinical.data))]
  
  plot = ggplot(plot_df, aes(x = Response, y = ICR_score, fill = Response), outlier.shape = NA) +
    ylim(-2, max(plot_df$ICR_score)*1.2) +
    facet_grid(rows = plot_df$Pathway) +
    geom_boxplot(outlier.shape = NA) +
    geom_jitter(position = position_jitter(width = 0.15)) +
    theme_bw() + #scale_fill_brewer(palette = "Set1") +
    stat_compare_means(method = "t.test", comparisons = list(c("nonresponse", "long-survival or response")),
                       size = 7, label.y = max(plot_df$ICR_score * 1.1)) + # remove comparisons when clinical outcome is different from group
    ylab("ICR score") + 
    scale_fill_manual(values=c("#E8A025", "#099E74")) +
    theme(
      strip.background = element_blank(),
      text = element_text(size = 20, colour = "black"),
      strip.text = element_text(angle = 0),
      axis.ticks.x = element_blank(),
      axis.text.x = element_blank()
    )
    
    
  png(filename = paste0("./Figures/vanAllen/ICR_Boxplots/", clinical_outcome, "/Facet_by_", Pathway, "_Boxplot_Overall_", ICR_readout,"_by_Response.png"), res = 600,
      width = 7, height = 5, units = "in")
      #width = 6.5, height = 6, units = "in")
  plot(plot)
  dev.off()
  
  dir.create("./Analysis/vanAllen/Tertiles", showWarnings = FALSE)
  save(tertiles.df, file = paste0("./Analysis/vanAllen/Tertiles/", clinical_outcome, "_", Pathway, "_", ICR_readout, "_tertiles_dataframe.Rdata"))
}

