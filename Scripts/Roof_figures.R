# Author: M, de Celis Rodriguez
# Date: 30/01/2024
# Project: MadreenRoof - Microbial Communities

library(vegan)
library(ggplot2)
library(cowplot)
library(hillR)
library(reshape2)
library(agricolae)
library(igraph)
library(piecewiseSEM)
library(Hmisc)
library(randomForest)
library(nlme)


rm(list=ls()) #Clear R environment

# Set the project location as working directory
setwd("~/../OneDrive/Proyecto - MadreenRoof/Roof/")

#
#### CUSTOM FUNCTIONS       ####
#Veech probability function for network calculation
veech_prob <- function(asv_df, conf.level = 0.05, type = "both") {
  
  ## asv_df:  ASV/OTU table with rows as samples and columns as ASV/OTU
  
  # Will compile C++ functions
  require(Rcpp)
  
  # Species co-occurrence table (presence-absence matrix multiplication)
  sourceCpp("Scripts/eigenMapMatMult.cpp") #C++
  cooc_df <- eigenMapMatMult(t(asv_df), asv_df, n_cores = 4)
  row.names(cooc_df) <- colnames(cooc_df) <- colnames(asv_df)
  
  # Allocate resources for results table
  calc_list <- vector("list", nrow(cooc_df)**2)
  
  # Total number of sites
  N_tot <- nrow(asv_df)
  j_tot <- 0:N_tot
  # Number of ways that sp1 and sp2 could be arranged among N_tot sites
  choose_N_tot <- choose(N_tot, j_tot)
  
  i <- 0
  
  if (type == "both") {
    
    print("Calculating Co-occurrences/Co-exclusions")
    
    # Calculate Pgl and Plt according to Veech (2013) for each pair of species
    for (sp1 in 1:nrow(cooc_df)) {
      
      # Number of sites occupied for species 1
      N_sp1 <- cooc_df[sp1,sp1]
      
      for (sp2 in 1:ncol(cooc_df)) {
        if (sp1 >= sp2) {
          next
        }
        
        # Number of sites occupied for species 2
        N_sp2 <- cooc_df[sp2,sp2]
        
        # Number of ways that sites having sp2 and not sp1 could be arranged among N_tot sites
        choose_sp2_tot <- choose(N_tot-j_tot, N_sp2-j_tot)
        # Number of ways that sites having sp1 could be arranged among sites not having sp2
        choose_sp1_tot <- choose(N_tot - N_sp2, N_sp1 - j_tot)
        
        # Number of samples where sp1 and sp2 co-occur
        N_sp1.2 <- cooc_df[sp1,sp2]
        
        # Number of ways that N_sp1.2 co-occurrence sites could be arranged among the N_tot sites
        # choose_N_tot[which(j_tot == N_sp1.2)]
        # Number of ways that sites having sp2 and not sp1 could be arranged among sites not having sp1 and sp2
        # choose_sp2_tot[which(j_tot == N_sp1.2)]
        # Number of ways that sites having sp1 could be arranged among sites not having sp2
        # choose_sp1_tot[which(j_tot == N_sp1.2)]
        
        ## NUMERATOR: Total number of ways that sp1 and sp2 could be distributed among N_tot sites for a given N_sp1, N_sp2 and N_sp1.2
        ## DENOMINATOR: The total number of ways that sp1 and sp2 can be arranged among N_tot sites without regard for N_sp1.2
        # Probability sp1 and sp2 co-occur at exactly N_sp1.2 sites
        P_sp1.sp2 <- (choose_N_tot[which(j_tot == N_sp1.2)] * choose_sp2_tot[which(j_tot == N_sp1.2)] * choose_sp1_tot[which(j_tot == N_sp1.2)]) / 
          (choose_N_tot[which(j_tot == N_sp1)] * choose_N_tot[which(j_tot == N_sp2)])
        
        
        # Probability of sp1 and sp2 co-occurring significantly more often than expected
        p.val_Pgt <- P_sp1.sp2 + sum((choose_N_tot[which(j_tot > N_sp1.2)] * choose_sp2_tot[which(j_tot > N_sp1.2)] * choose_sp1_tot[which(j_tot > N_sp1.2)]) / 
                                       (choose_N_tot[which(j_tot == N_sp1)] * choose_N_tot[which(j_tot == N_sp2)]))
        
        if (p.val_Pgt <= conf.level) {
          
          i <- i+1
          calc_list[[i]] <- data.frame(sp1 = row.names(cooc_df)[sp1], sp2 = row.names(cooc_df)[sp2],
                                       p.value = p.val_Pgt, type = "Co-occurrence")
          
        } else {
          
          # Probability of sp1 and sp2 co-occurring significantly less often than expected
          p.val_Plt <- P_sp1.sp2 + sum((choose_N_tot[which(j_tot < N_sp1.2)] * choose_sp2_tot[which(j_tot < N_sp1.2)] * choose_sp1_tot[which(j_tot < N_sp1.2)]) / 
                                         (choose_N_tot[which(j_tot == N_sp1)] * choose_N_tot[which(j_tot == N_sp2)]))
          
          if (p.val_Plt <= conf.level) {
            
            i <- i+1
            calc_list[[i]] <- data.frame(sp1 = row.names(cooc_df)[sp1], sp2 = colnames(cooc_df)[sp2],
                                         p.value = p.val_Plt, type = "Co-exclusion")
          }
          
        }
        
      }
      print(sp1*100/nrow(cooc_df))
    }
    
  } else if (type == "Co-occurrences") {
    
    print("Calculating Co-occurrences")
    
    # Calculate Pgl and Plt according to Veech (2013) for each pair of species
    for (sp1 in 1:nrow(cooc_df)) {
      
      # Number of sites occupied for species 1
      N_sp1 <- cooc_df[sp1,sp1]
      
      for (sp2 in 1:ncol(cooc_df)) {
        if (sp1 >= sp2) {
          next
        }
        
        # Number of sites occupied for species 2
        N_sp2 <- cooc_df[sp2,sp2]
        
        # Number of ways that sites having sp2 and not sp1 could be arranged among N_tot sites
        choose_sp2_tot <- choose(N_tot-j_tot, N_sp2-j_tot)
        # Number of ways that sites having sp1 could be arranged among sites not having sp2
        choose_sp1_tot <- choose(N_tot - N_sp2, N_sp1 - j_tot)
        
        # Number of samples where sp1 and sp2 co-occur
        N_sp1.2 <- cooc_df[sp1,sp2]
        
        # Number of ways that N_sp1.2 co-occurrence sites could be arranged among the N_tot sites
        # choose_N_tot[which(j_tot == N_sp1.2)]
        # Number of ways that sites having sp2 and not sp1 could be arranged among sites not having sp1 and sp2
        # choose_sp2_tot[which(j_tot == N_sp1.2)]
        # Number of ways that sites having sp1 could be arranged among sites not having sp2
        # choose_sp1_tot[which(j_tot == N_sp1.2)]
        
        ## NUMERATOR: Total number of ways that sp1 and sp2 could be distributed among N_tot sites for a given N_sp1, N_sp2 and N_sp1.2
        ## DENOMINATOR: The total number of ways that sp1 and sp2 can be arranged among N_tot sites without regard for N_sp1.2
        # Probability sp1 and sp2 co-occur at exactly N_sp1.2 sites
        P_sp1.sp2 <- (choose_N_tot[which(j_tot == N_sp1.2)] * choose_sp2_tot[which(j_tot == N_sp1.2)] * choose_sp1_tot[which(j_tot == N_sp1.2)]) / 
          (choose_N_tot[which(j_tot == N_sp1)] * choose_N_tot[which(j_tot == N_sp2)])
        
        
        # Probability of sp1 and sp2 co-occurring significantly more often than expected
        p.val_Pgt <- P_sp1.sp2 + sum((choose_N_tot[which(j_tot > N_sp1.2)] * choose_sp2_tot[which(j_tot > N_sp1.2)] * choose_sp1_tot[which(j_tot > N_sp1.2)]) / 
                                       (choose_N_tot[which(j_tot == N_sp1)] * choose_N_tot[which(j_tot == N_sp2)]))
        
        if (p.val_Pgt <= conf.level) {
          
          i <- i+1
          calc_list[[i]] <- data.frame(sp1 = row.names(cooc_df)[sp1], sp2 = row.names(cooc_df)[sp2],
                                       p.value = p.val_Pgt, type = "Co-occurrence")
          
        } 
      }
      print(sp1*100/nrow(cooc_df))
    }
    
  } else if (type == "Co-exclusions") {
    
    print("Calculating Co-exclusions")
    
    # Calculate Pgl and Plt according to Veech (2013) for each pair of species
    for (sp1 in 1:nrow(cooc_df)) {
      
      # Number of sites occupied for species 1
      N_sp1 <- cooc_df[sp1,sp1]
      
      for (sp2 in 1:ncol(cooc_df)) {
        if (sp1 >= sp2) {
          next
        }
        
        # Number of sites occupied for species 2
        N_sp2 <- cooc_df[sp2,sp2]
        
        # Number of ways that sites having sp2 and not sp1 could be arranged among N_tot sites
        choose_sp2_tot <- choose(N_tot-j_tot, N_sp2-j_tot)
        # Number of ways that sites having sp1 could be arranged among sites not having sp2
        choose_sp1_tot <- choose(N_tot - N_sp2, N_sp1 - j_tot)
        
        # Number of samples where sp1 and sp2 co-occur
        N_sp1.2 <- cooc_df[sp1,sp2]
        
        # Number of ways that N_sp1.2 co-occurrence sites could be arranged among the N_tot sites
        # choose_N_tot[which(j_tot == N_sp1.2)]
        # Number of ways that sites having sp2 and not sp1 could be arranged among sites not having sp1 and sp2
        # choose_sp2_tot[which(j_tot == N_sp1.2)]
        # Number of ways that sites having sp1 could be arranged among sites not having sp2
        # choose_sp1_tot[which(j_tot == N_sp1.2)]
        
        ## NUMERATOR: Total number of ways that sp1 and sp2 could be distributed among N_tot sites for a given N_sp1, N_sp2 and N_sp1.2
        ## DENOMINATOR: The total number of ways that sp1 and sp2 can be arranged among N_tot sites without regard for N_sp1.2
        # Probability sp1 and sp2 co-occur at exactly N_sp1.2 sites
        P_sp1.sp2 <- (choose_N_tot[which(j_tot == N_sp1.2)] * choose_sp2_tot[which(j_tot == N_sp1.2)] * choose_sp1_tot[which(j_tot == N_sp1.2)]) / 
          (choose_N_tot[which(j_tot == N_sp1)] * choose_N_tot[which(j_tot == N_sp2)])
        
        
        # Probability of sp1 and sp2 co-occurring significantly less often than expected
        p.val_Plt <- P_sp1.sp2 + sum((choose_N_tot[which(j_tot < N_sp1.2)] * choose_sp2_tot[which(j_tot < N_sp1.2)] * choose_sp1_tot[which(j_tot < N_sp1.2)]) / 
                                       (choose_N_tot[which(j_tot == N_sp1)] * choose_N_tot[which(j_tot == N_sp2)]))
        
        if (p.val_Plt <= conf.level) {
          
          i <- i+1
          calc_list[[i]] <- data.frame(sp1 = row.names(cooc_df)[sp1], sp2 = colnames(cooc_df)[sp2],
                                       p.value = p.val_Plt, type = "Co-exclusion")
        }
      }
      print(sp1*100/nrow(cooc_df))
    }
    
  } else if (!type %in% c("both", "Co-occurrence", "Co-exclusions")) {
    
    print("Only \"both\", \"Co-occurrences\" and \"Co-exclusions\" are currently supported")
    
  }
  
  print("Merging List to a Matrix")
  calc_df <- data.table::rbindlist(calc_list)
  return(calc_df)
  
}

#Set Unique Unidentified
set_unid <- function(tax_df) {
  
  for (i in 1:ncol(tax_df)) {
    if (i != 7) {
      tax_df[,i] <- ifelse(is.na(tax_df[,i]),
                           ifelse(grepl("Unidentified", tax_df[,i-1]) == TRUE,
                                  tax_df[,i-1],
                                  paste("Unidentified", tax_df[,i-1], sep = " ")),
                           tax_df[,i])
    } else{
      tax_df[,i] <- ifelse(is.na(tax_df[,i]),
                           ifelse(grepl("Unidentified", tax_df[,i-1]) == TRUE,
                                  tax_df[,i-1],
                                  "sp."),
                           tax_df[,i])
      
    }
    
  }
  
  return(tax_df)
  
}

#Rarefaction
rare_list <- function(asv_df, depth, bootstrap = 100) {
  
  asv.r_df <- vector(mode = "list")
  for (boot in 1:bootstrap) {
    set.seed(boot)
    asv.r_df[[boot]] <- rrarefy(asv_df[rowSums(asv_df) >= depth, ], depth)
    
    print(boot*100/bootstrap)
    
  }
  
  return(asv.r_df)
  
}


#
#### DATA LOADING           ####

## Sample data
sample_df <- read.table("Data/sample_df.txt", sep = "\t", header = TRUE, encoding = "latin1", na.strings = "")
sample_df$sample_time <- factor(sample_df$sample_time, levels = c("Initial", "Final"))
row.names(sample_df) <- sample_df$ID_well
sample_df$Material <- factor(sample_df$Material, levels = c("SW", "CS", "SCG", "Peat"))
sample_df$Treatment <- factor(sample_df$Treatment, levels = c("SW", "SW-BC", "CS", "CS-BC", "SCG", "SCG-BC", "Peat"))

colnames(sample_df)[13] <- "EC"

## Community data
#Bacteria
asv_prok <- as.data.frame(readRDS("Data/Sequencing/Outputs/ASV_prok.rds"))

tax_prok <- readRDS("Data/Sequencing/Outputs/tax_prok.rds")
tax.ud_prok <- tax_prok
tax.ud_prok[tax.ud_prok == "Unidentified"] <- NA
tax.ud_prok <- set_unid(tax.ud_prok)


#Fungi
asv_fun <- as.data.frame(readRDS("Data/Sequencing/Outputs/ASV_fun.rds"))

tax_fun <- readRDS("Data/Sequencing/Outputs/tax_fun.rds")
tax.ud_fun <- tax_fun
tax.ud_fun[tax.ud_fun == "Unidentified"] <- NA
tax.ud_fun <- set_unid(tax.ud_fun)
tax.ud_fun <- apply(tax.ud_fun, 2, function(x) gsub("[a-z]__", "", x))
tax.ud_fun <- as.data.frame(tax.ud_fun)

## Colors
col_mat <- c("#613b25", "#facb7a", "#2ac219", "#ae0e36")
names(col_mat) <- c("SCG", "CS", "SW", "Peat")

col_trt <- c("#613b25", "#453227", "#facb7a", "#cfb17e", "#2ac219", "#43803c", "#ae0e36")
names(col_trt) <- c("SCG", "SCG-BC", "CS", "CS-BC", "SW", "SW-BC", "Peat")

col_t <- c("#000000", "#999999")
names(col_t) <- c("Initial", "Final")

#
################################################################################ FIGURE 1                     ####
################################### SUBSTRATE DIFFERENCES        ####
#### OVERALL NUTRIENT POOL  ####

## PRINCIPAL COMPONENT ANALYSIS - PCA

soil_pca <- prcomp(sample_df[complete.cases(sample_df[,8:13]),8:13], scale. = TRUE)

soil.pca_df <- scores(soil_pca, display = "sites")
soil.pca_df <- merge(sample_df[,1:7], soil.pca_df, by = "row.names")
soil.pca_df <- soil.pca_df[,-1]

soil.pca_vec <- as.data.frame(scores(soil_pca, display = "species"))
soil.pca_vec$property <- row.names(soil.pca_vec)

soil.pca_vec2 <- soil.pca_vec
soil.pca_vec2["P", "PC2"] <- soil.pca_vec2["P", "PC2"] + 0.02
soil.pca_vec2["EC", "PC1"] <- soil.pca_vec2["EC", "PC1"] - 0.02

gg.pca_soil <- ggplot() +
  geom_segment(data = soil.pca_vec, aes(x = 0, y = 0, xend = PC1*5, yend = PC2*5), 
               arrow = arrow(length = unit(0.2, "cm")), size = 1, alpha = 0.8) + 
  geom_line(data = soil.pca_df, aes(x = PC1, y = PC2, group = ID_block), color = "gray", show.legend = FALSE) +
  geom_point(data = soil.pca_df, aes(x = PC1, y = PC2, fill = Material, color = sample_time, shape = Biochar, size = Biochar), stroke = 1.5) +
  geom_text(data = soil.pca_vec2, aes(x = PC1*5.25, y = PC2*5.25, label = property), size = 6) +
  scale_shape_manual(values = c(21, 23)) +
  scale_size_manual(values = c(4, 3), guide = "none") +
  scale_fill_manual(values = col_mat) +
  scale_color_manual(values = c("#000000", "#999999")) +
  xlab(paste("PC1 (", round(summary(soil_pca)$importance[2,]*100,2)[1], "%)", sep = "")) + 
  ylab(paste("PC2 (", round(summary(soil_pca)$importance[2,]*100,2)[2], "%)", sep = "")) + 
  guides(shape = guide_legend(override.aes = list(size = c(3, 2), stroke = 2)),
         color = guide_legend(title = "Time", override.aes = list(shape = 21, size = 3, stroke = 2)),
         fill = guide_legend(override.aes = list(shape = 21, size = 3, color = NA))) +
  theme_bw() +
  theme(legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"))

gg.pca_soil

set.seed(1)
adonis2(scale(sample_df[complete.cases(sample_df[,8:13]),8:13]) ~ Material + Biochar + sample_time, 
        data = sample_df[complete.cases(sample_df[,8:13]),], method = "eu")

#
#### VARIANCE PARTITIONING  ####

substrate_initial <- subset(sample_df, sample_time == "Initial")
substrate_initial.dist <- as.matrix(vegdist(scale(substrate_initial[,8:13]), method = "euclidean"))

## Variance partitioning analysis
varp_initial <- varpart(substrate_initial.dist, 
                        ~ Material,
                        ~ Biochar,
                        data = substrate_initial[row.names(substrate_initial.dist),])

## Set results in a graphicable format
varp.plot_initial <- varp_initial$part$indfract

varp.plot_initial$Adj.R.squared <- (varp.plot_initial$Adj.R.squared - min(varp.plot_initial$Adj.R.squared)) /
  sum(varp.plot_initial$Adj.R.squared - min(varp.plot_initial$Adj.R.squared))

varp.plot_initial$variable <- c("Compost", "Biochar", "Shared", "Residuals")
varp.plot_initial$variable <- factor(varp.plot_initial$variable, levels = c("Compost", "Biochar", "Shared", "Residuals"))



substrate_final <- subset(sample_df, sample_time == "Final")
substrate_final.dist <- as.matrix(vegdist(scale(substrate_final[-15,8:13]), method = "euclidean"))

## Variance partitioning analysis
varp_final <- varpart(substrate_final.dist, 
                      ~ Material,
                      ~ Biochar,
                      data = substrate_final[row.names(substrate_final.dist),])

## Set results in a graphicable format
varp.plot_final <- varp_final$part$indfract

varp.plot_final$Adj.R.squared <- (varp.plot_final$Adj.R.squared - min(varp.plot_final$Adj.R.squared)) /
  sum(varp.plot_final$Adj.R.squared - min(varp.plot_final$Adj.R.squared))

varp.plot_final$variable <- c("Compost", "Biochar", "Shared", "Residuals")
varp.plot_final$variable <- factor(varp.plot_final$variable, levels = c("Compost", "Biochar", "Shared", "Residuals"))


varp.plot_total <- rbind(cbind.data.frame(Time = "Initial", varp.plot_initial),
                         cbind.data.frame(Time = "Final", varp.plot_final))

varp.plot_total$Time <- factor(varp.plot_total$Time, levels = c("Initial", "Final"))

gg.varp_functioning <- ggplot(data = varp.plot_total) +
  geom_bar(aes(x = Time, y = Adj.R.squared, fill = variable), stat = "identity", position = position_fill(reverse = TRUE)) +
  scale_fill_manual(values = c("#d5954f", "#00954f", "#726200", "#4d4d4d")) +
  guides(fill = guide_legend(nrow = 2)) +
  ylab("Variance explained") +
  theme_minimal() +
  theme(legend.position = "bottom",
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        legend.title = element_blank(), 
        legend.text = element_text(size = 15, color = "black"),
        panel.grid.major.y = element_line(linewidth = 1, color = "gray80"),
        panel.grid.minor.y = element_line(linewidth = 0.75, color = "gray90"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank())

gg.varp_functioning

#
#### TOMATO PROPERTIES      ####

tomato_df <- melt(sample_df[complete.cases(sample_df),c(1:7, 14,15)], 
                  id.vars = c("ID_well", "ID_block", "ID_sample", "sample_time", "Biochar", "Material", "Treatment"))

tomato_df$Treatment <- factor(tomato_df$Treatment, 
                              levels = c("SW", "SW-BC", "CS", "CS-BC", "SCG", "SCG-BC", "Peat"))


stats_df <- NULL
for (var in unique(tomato_df$variable)) {
  
  stats_sub <- merge(aggregate(value ~ Treatment, subset(tomato_df, variable == var), function(x) quantile(x)[4]),
                     LSD.test(aov(value ~ Treatment, subset(tomato_df, variable == var)), trt = "Treatment")$groups[,2, drop = FALSE],
                     by.x = "Treatment", by.y = "row.names")
  
  stats_df <- rbind(stats_df, cbind.data.frame(variable = var, stats_sub)) 
  
}
stats_df$pos <- ifelse(stats_df$variable == "Height", stats_df$value + 10, stats_df$value + 0.5)

label <- c("Height (cm)", "Yield (Kg)")
names(label) <- c("Height", "Yield")

gg.box_tomato <- ggplot(data = tomato_df) +
  geom_point(aes(x = Treatment, y = value, color = Treatment), position = position_jitterdodge(), size = 3, show.legend = FALSE) +
  geom_boxplot(aes(x = Treatment, y = value, color = Treatment), alpha = 0.25, size = 2, show.legend = FALSE) +
  geom_label(data = stats_df, aes(x = Treatment, y = pos, label = groups), 
             color = "black", fill = NA, label.size = NA, size = 7) +
  facet_wrap(~ variable, nrow = 1, scales = "free_y", labeller = as_labeller(label)) +
  scale_color_manual(values = col_trt) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16, color = "black"),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 16, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"))

gg.box_tomato

#
################################### EXPORT FIGURE                ####

gg.figure1 <- plot_grid(plot_grid(gg.pca_soil, gg.varp_functioning, rel_widths = c(1, 0.45), labels = c("A", "B"), label_size = 18), 
                        gg.box_tomato, ncol = 1, rel_heights = c(1, 0.5), labels = c("", "C"), label_size = 18)
gg.figure1

ggsave("../Manuscript/Figures/Figure_1.png", gg.figure1, bg = "white", width = 12, height = 10)

#
################################################################################ SUPPLEMENTARY FIGURE S2      ####
#### OVERALL NUTRIENT POOL  ####

## BOXPLOT

sample_df$C_N <- sample_df$C / sample_df$N

sample_df_try <- sample_df
sample_df_try$empty <- NA


soil_df <- melt(sample_df_try[,c(1:13,16:17)], id.vars = c("ID_well", "ID_block", "ID_sample", "sample_time", "Biochar", "Material", "Treatment"))
soil_df$Treatment <- factor(soil_df$Treatment, 
                            levels = c("SW", "SW-BC", "CS", "CS-BC", "SCG", "SCG-BC", "Peat"))

soil_df$variable <- factor(soil_df$variable, levels = c("C", "N", "C_N", "empty", "P", "K", "pH", "EC"))


stats_df <- NULL
for (var in c("C", "N", "C_N", "P", "K", "pH", "EC")) {
  for (trt in unique(soil_df$Treatment)) {
    
    stats_sub0 <- subset(soil_df, variable == var & Treatment == trt & sample_time == "Initial")
    stats_sub6 <- subset(soil_df, variable == var & Treatment == trt & sample_time == "Final")
    
    stats_sub <- merge(stats_sub0[,c(2,6:9)], stats_sub6[,c(2,9)], by = "ID_block")
    stats_sub.df <- cbind.data.frame(max = max(max(stats_sub$value.x, na.rm = TRUE), max(stats_sub$value.y, na.rm = TRUE)),
                                     position = max(quantile(stats_sub$value.x, na.rm = TRUE)[4], quantile(stats_sub$value.y, na.rm = TRUE)[4]),
                                     p_value = wilcox.test(stats_sub$value.x, stats_sub$value.y)$p.value)
    
    stats_df <- rbind(stats_df, cbind.data.frame(variable = var, Treatment = trt, stats_sub.df)) 
    
  }
}

stats_df$variable <- factor(stats_df$variable, levels = c("C", "N", "C_N", "empty", "P", "K", "pH", "EC"))

label <- c("C (%)", "N (%)", "C/N", "", "P (ppm)", "K (ppm)", "pH", "EC (mS/m)")
names(label) <- c("C", "N", "C_N", "empty", "P", "K", "pH", "EC")

gg.box_soil <- ggplot(data = soil_df) +
  geom_point(aes(x = Treatment, y = value, color = sample_time), position = position_jitterdodge(), size = 2) +
  geom_boxplot(aes(x = Treatment, y = value, color = sample_time), size = 1, alpha = 0.75) +
  geom_label(data = stats_df, aes(x = Treatment, y = position+0.025*max, label = ifelse(p_value < 0.05, "*", "")), 
             color = "black", fill = NA, label.size = NA, size = 7) +
  guides(color = guide_legend(title = "Time")) +
  facet_wrap(~ variable, ncol = 2, scales = "free_y", labeller = as_labeller(label)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 16, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"))

gg.box_soil

#
################################### EXPORT FIGURE                ####

gg.figureS2 <- gg.box_soil
gg.figureS2

ggsave("../Manuscript/Figures/Supp_Figure_S2.pdf", gg.figureS2, bg = "white", width = 8, height = 10)

#
################################################################################ SUPPLEMENTARY FIGURE S3      ####
#### TOMATO YIELD vs HEIGH  ####

gg.height_yield <- ggplot(data = substrate_final) +
  geom_point(aes(x = Height, y = Yield, color = Material, shape = Biochar, size = Biochar)) +
  scale_shape_manual(values = c(16, 18)) +
  scale_size_manual(values = c(4, 5), guide = "none") +
  scale_color_manual(values = col_mat) +
  guides(color = guide_legend(override.aes = list(size = 4), order = 1),
         shape = guide_legend(override.aes = list(size = c(3, 4)), order = 2)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"))

gg.height_yield

#
#### INITIAL vs FINAL  ####

substrate_df <- merge(subset(sample_df, sample_time == "Initial")[,c(2,5:13,16)], 
                      subset(sample_df, sample_time == "Final")[,c(2,5:13,16)], by = c("ID_block", "Biochar", "Material", "Treatment"))

substrate_df <- cbind.data.frame(substrate_df[,1:4], C = substrate_df$C.y/substrate_df$C.x, N = substrate_df$N.y/substrate_df$N.x,
                                 P = substrate_df$P.y/substrate_df$P.x, K = substrate_df$P.y/substrate_df$P.x, 
                                 pH = substrate_df$pH.y/substrate_df$pH.x, EC = substrate_df$EC.y/substrate_df$EC.x, 
                                 C_N = substrate_df$C_N.y/substrate_df$C_N.x)

substrate_df$empty <- NA

substrate_df.plot <- melt(substrate_df, id.vars = c("ID_block", "Biochar", "Material", "Treatment"))
substrate_df.plot$variable <- factor(substrate_df.plot$variable, levels = c("C", "N", "C_N", "empty", "P", "K", "pH", "EC"))

stats_df <- NULL
for (var in unique(substrate_df.plot$variable)[-8]) {
  
  stats_sub <- merge(aggregate(value ~ Treatment, subset(substrate_df.plot, variable == var), function(x) quantile(x)[4]),
                     LSD.test(aov(value ~ Treatment, subset(substrate_df.plot, variable == var)), trt = "Treatment")$groups[,2, drop = FALSE],
                     by.x = "Treatment", by.y = "row.names")
  
  stats_df <- rbind(stats_df, cbind.data.frame(variable = var, stats_sub)) 
  
}

stats_df <- rbind.data.frame(stats_df, cbind.data.frame(variable = "empty", Treatment = "SW", value = NA, groups = NA))
stats_df$variable <- factor(stats_df$variable, levels = c("C", "N", "C_N", "empty", "P", "K", "pH", "EC"))

gg.substrate_cultivation <- ggplot(data = substrate_df.plot) +
  geom_point(aes(x = Treatment, y = value, color = Treatment), position = position_jitterdodge(), alpha = 0.25, size = 2, show.legend = FALSE) +
  geom_boxplot(aes(x = Treatment, y = value, color = Treatment), size = 1, show.legend = FALSE) +
  geom_label(data = stats_df, aes(x = Treatment, y = value*1.1, label = groups), 
             color = "black", fill = NA, label.size = NA, size = 5) +
  facet_wrap(~ variable, ncol = 2, scales = "free_y", labeller = as_labeller(label)) +
  scale_color_manual(values = col_trt) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 16, color = "black"),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 16, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"))

gg.substrate_cultivation

#
################################### EXPORT FIGURE                ####

gg.figureS3 <- plot_grid(gg.height_yield, gg.substrate_cultivation, ncol = 1, rel_heights = c(0.5, 1), labels = c("A", "B"), label_size = 18)
gg.figureS3

ggsave("../Manuscript/Figures/Supp_Figure_S3.pdf", gg.figureS3, bg = "white", width = 8, height = 10)

#
################################################################################ SUPPLEMENTARY FIGURE S1      ####
################################### RAREFACTION                  ####
#### BACTERIA               ####

rare_prok <- rarecurve(asv_prok, step = 500, 40000, xlab = "Sample Size", ylab = "Species", label = TRUE)

rare_prok.plot <- NULL
for (n in 1:length(rare_prok)) {
  rare_prok.plot <- rbind(rare_prok.plot, cbind.data.frame(Sample = n, N.reads = attributes(rare_prok[[n]])$Subsample, N.ASVs = rare_prok[[n]]))
}

rare_prok.plot <- merge(rare_prok.plot, aggregate(N.reads ~ Sample, rare_prok.plot, max), by = "Sample")
colnames(rare_prok.plot) <- c("Sample", "N.reads", "N.ASVs", "max.N.reads")

gg.rare_prok <- ggplot(rare_prok.plot, aes(x = N.reads, y = N.ASVs, group = Sample, color = ifelse(max.N.reads < 36997, "red", "black"))) + 
  geom_vline(xintercept = 36997, linewidth = 1, linetype = "dashed", color = "blue") +
  geom_line(linewidth = 1, show.legend = FALSE) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Sample Size", y = "Species", title = "Bacteria") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", color = "black"),
        axis.title = element_text(size = 17, color = "black"),
        axis.text = element_text(size = 15, color = "black"))

gg.rare_prok

#asv.r_prok <- rare_list(asv_prok, 36997, bootstrap = 100) 
#saveRDS(asv.r_prok, "Outputs/asv.r_prok.rds")
asv.r_prok <- readRDS("Outputs/asv.r_prok.rds")

#
#### FUNGI                  ####

rare_fun <- rarecurve(asv_fun, step = 500, 40000, xlab = "Sample Size", ylab = "Species", label = TRUE)

rare_fun.plot <- NULL
for (n in 1:length(rare_fun)) {
  rare_fun.plot <- rbind(rare_fun.plot, cbind.data.frame(Sample = n, N.reads = attributes(rare_fun[[n]])$Subsample, N.ASVs = rare_fun[[n]]))
}

rare_fun.plot <- merge(rare_fun.plot, aggregate(N.reads ~ Sample, rare_fun.plot, max), by = "Sample")
colnames(rare_fun.plot) <- c("Sample", "N.reads", "N.ASVs", "max.N.reads")

gg.rare_fun <- ggplot(rare_fun.plot, aes(x = N.reads, y = N.ASVs, group = Sample, color = ifelse(max.N.reads < 47909, "red", "black"))) + 
  geom_vline(xintercept = 47909, linewidth = 1, linetype = "dashed", color = "blue") +
  geom_line(linewidth = 1, show.legend = FALSE) +
  scale_color_manual(values = c("black", "red")) +
  labs(x = "Sample Size", y = "Species", title = "Fungi") +
  theme_bw() +
  theme(plot.title = element_text(size = 18, face = "bold", color = "black"),
        axis.title = element_text(size = 17, color = "black"),
        axis.text = element_text(size = 15, color = "black"))

gg.rare_fun

#asv.r_fun2 <- rare_list(asv_fun, 47909, bootstrap = 100) 
#saveRDS(asv.r_fun, "Outputs/asv.r_fun.rds")
asv.r_fun <- readRDS("Outputs/asv.r_fun.rds")

#
################################### EXPORT FIGURE                ####

rare_tot.plot <- rbind(cbind.data.frame(Group = "Bacteria", Threshold = 36997, rare_prok.plot),
                       cbind.data.frame(Group = "Fungi", Threshold = 47909, rare_fun.plot))
rare_tot.plot$filter <- ifelse(rare_tot.plot$max.N.reads >= rare_tot.plot$Threshold, "Kept", "Wiped")

gg.rare_tot <- ggplot(rare_tot.plot, aes(x = N.reads, y = N.ASVs, group = Sample, color = filter)) + 
  geom_vline(aes(xintercept = Threshold), linewidth = 1, linetype = "dashed", color = "blue") +
  geom_line(linewidth = 1, show.legend = FALSE) +
  scale_color_manual(values = c("black", "red")) +
  facet_wrap(~ Group, nrow = 1, scales = "free") +
  labs(x = "Sample Size", y = "Species") +
  theme_bw() +
  theme(axis.title = element_text(size = 17, color = "black"),
        axis.text = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 17, color = "black", hjust = 0, face = "bold"),
        strip.background = element_rect(fill = "white", linewidth = NA))

gg.rare_tot

gg.figureS1 <- gg.rare_tot

ggsave("../Manuscript/Figures/Supp_Figure_S1.png", gg.figureS1, bg = "white", width = 12, height = 4)

#
################################################################################ FIGURE 2                     ####
################################### BETA DIVERSITY               ####
#### BACTERIA               ####

bray_prok <- as.matrix(Reduce("+",lapply(asv.r_prok, function(x) vegdist(x, method = "bray", binary = FALSE)))/100)

set.seed(1)
nMDS_prok <- metaMDS(bray_prok)
nMDS_prok$stress

nMDS_prok.plot <- as.data.frame(nMDS_prok[["points"]])
nMDS_prok.plot$ID_well <- rownames(nMDS_prok.plot)
nMDS_prok.plot <- merge(nMDS_prok.plot, sample_df, by = "ID_well")

nMDS_prok.plot$sample_time <- factor(nMDS_prok.plot$sample_time, levels = c("Initial", "Final"))

gg.nmds.bray_prok <- ggplot(nMDS_prok.plot) + 
  geom_line(aes(x = MDS1, y = MDS2, fill = Material, group = ID_block), show.legend = FALSE) +
  geom_point(aes(x = MDS1, y = MDS2, fill = Material, shape = Biochar, size = Biochar, color = sample_time), 
             stroke = 1.5, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 23)) +
  scale_size_manual(values = c(4, 3), guide = "none") +
  scale_fill_manual(values = col_mat) +
  scale_color_manual(values = c("#000000", "#999999")) +
  annotate(geom = "text", x = Inf, y = Inf, label = paste("stress =", round(nMDS_prok$stress, 3)), size = 6, vjust = 1.2, hjust = 1.05) +
  theme_bw() +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  theme(plot.title = element_text(size = 18, color = "black", face = "bold"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"))

gg.nmds.bray_prok

#
#### FUNGI                  ####

bray_fun <- as.matrix(Reduce("+",lapply(asv.r_fun, function(x) vegdist(x, method = "bray", binary = FALSE)))/100)

set.seed(1)
nMDS_fun <- metaMDS(bray_fun)
nMDS_fun$stress

nMDS_fun.plot <- as.data.frame(nMDS_fun[["points"]])
nMDS_fun.plot$ID_well <- rownames(nMDS_fun.plot)
nMDS_fun.plot <- merge(nMDS_fun.plot, sample_df, by = "ID_well")

nMDS_fun.plot$sample_time <- factor(nMDS_fun.plot$sample_time, levels = c("Initial", "Final"))

gg.nmds.bray_fun <- ggplot(nMDS_fun.plot) + 
  geom_line(aes(x = MDS1, y = MDS2, fill = Material, group = ID_block), show.legend = FALSE) +
  geom_point(aes(x = MDS1, y = MDS2, fill = Material, shape = Biochar, size = Biochar, color = sample_time), 
             stroke = 1.5, show.legend = FALSE) +
  scale_shape_manual(values = c(21, 23)) +
  scale_size_manual(values = c(4, 3), guide = "none") +
  scale_fill_manual(values = col_mat) +
  scale_color_manual(values = c("#000000", "#999999")) +
  annotate(geom = "text", x = Inf, y = Inf, label = paste("stress =", round(nMDS_fun$stress, 3)), size = 6, vjust = 1.2, hjust = 1.05) +
  theme_bw() +
  guides(shape = guide_legend(override.aes = list(size = c(3, 1))),
         color = guide_legend(override.aes = list(shape = 16, size = 3))) +
  theme(plot.title = element_text(size = 18, color = "black", face = "bold"),
        axis.text.y = element_text(size = 16, color = "black"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 16, color = "black"))

gg.nmds.bray_fun

#
#### STATISTICS             ####

## Bacteria
set.seed(1)
adonis2(bray_prok ~ Material * Biochar + sample_time, data = sample_df[row.names(bray_prok),])
set.seed(1)
adonis2(bray_prok ~ Biochar, data = sample_df[row.names(bray_prok),])


## Fungi
set.seed(1)
adonis2(bray_fun ~ Material * Biochar + sample_time, data = sample_df[row.names(bray_fun),])
set.seed(1)
adonis2(bray_fun ~ Biochar, data = sample_df[row.names(bray_fun),])

#
################################### ALPHA DIVERSITY              ####
#### BACTERIA               ####

alpha_prok <- cbind.data.frame(Richness = Reduce("+",lapply(asv.r_prok, function(x) hill_taxa(x, q = 0)))/100,
                               `Hill–Shannon` = Reduce("+",lapply(asv.r_prok, function(x) hill_taxa(x, q = 1)))/100,
                               `Hill–Simpson` = Reduce("+",lapply(asv.r_prok, function(x) hill_taxa(x, q = 2)))/100)

alpha_prok$ID_well <- row.names(alpha_prok)
alpha_prok <- merge(alpha_prok, sample_df, by = "ID_well")

alpha_prok.plot <- melt(alpha_prok[,1:10], id.vars = c("ID_well", "ID_block", "ID_sample", "sample_time",
                                                       "Biochar", "Material", "Treatment"))

alpha_prok.plot$sample_time <- factor(alpha_prok.plot$sample_time, levels = c("Initial", "Final"))

## BOXPLOT
gg.alpha_prok <- ggplot(data = alpha_prok.plot) +
  geom_point(aes(x = Treatment, y = value, fill = Treatment, color = sample_time), 
             size = 3, shape = 21, position = position_jitterdodge(), show.legend = FALSE) +  
  geom_boxplot(aes(x = Treatment, y = value, fill = Treatment, color = sample_time), 
             position = position_dodge(), alpha = 0.75, outlier.shape = NA) +
  scale_fill_manual(values = col_trt) +
  scale_color_manual(values = col_t) +
  facet_wrap( ~ variable, scales = "free_y", ncol = 1) +
  ylab("") + xlab("") + 
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"))

gg.alpha_prok

#
#### FUNGI                  ####

alpha_fun <- cbind.data.frame(Richness = Reduce("+",lapply(asv.r_fun, function(x) hill_taxa(x, q = 0)))/100,
                               `Hill–Shannon` = Reduce("+",lapply(asv.r_fun, function(x) hill_taxa(x, q = 1)))/100,
                               `Hill–Simpson` = Reduce("+",lapply(asv.r_fun, function(x) hill_taxa(x, q = 2)))/100)

alpha_fun$ID_well <- row.names(alpha_fun)
alpha_fun <- merge(alpha_fun, sample_df, by = "ID_well")

alpha_fun.plot <- melt(alpha_fun[,1:10], id.vars = c("ID_well", "ID_block", "ID_sample", "sample_time",
                                                       "Biochar", "Material", "Treatment"))

alpha_fun.plot$sample_time <- factor(alpha_fun.plot$sample_time, levels = c("Initial", "Final"))

## BOXPLOT
gg.alpha_fun <- ggplot(data = alpha_fun.plot) +
  geom_point(aes(x = Treatment, y = value, fill = Treatment, color = sample_time), 
             size = 3, shape = 21, position = position_jitterdodge(), show.legend = FALSE) +  
  geom_boxplot(aes(x = Treatment, y = value, fill = Treatment, color = sample_time), 
               position = position_dodge(), alpha = 0.75, outlier.shape = NA) +
  scale_fill_manual(values = col_trt) +
  scale_color_manual(values = col_t) +
  facet_wrap( ~ variable, scales = "free_y", ncol = 1) +
  ylab("") + xlab("") + 
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"))

gg.alpha_fun

#
#### STATISTICS             ####

alpha_tot.plot <- rbind(cbind.data.frame(Group = "Bacteria", alpha_prok.plot),
                        cbind.data.frame(Group = "Fungi", alpha_fun.plot))

alpha_tot.stat <- merge(subset(alpha_tot.plot, sample_time == "Initial")[,c(1,3,6:10)], subset(alpha_tot.plot, sample_time == "Final")[,c(1,3,6:10)], 
                        by = c("Group", "ID_block", "Biochar", "Material", "Treatment", "variable"))

colnames(alpha_tot.stat)[7:8] <- c("value_0", "value_6") 
alpha_tot.stat$Increase <- alpha_tot.stat$value_6 / alpha_tot.stat$value_0


## Bacteria
wilcox.test(subset(alpha_tot.stat, variable == "Richness" & Group == "Bacteria")$value_0,
            subset(alpha_tot.stat, variable == "Richness" & Group == "Bacteria")$value_6, alternative = "less", paired = TRUE)
wilcox.test(subset(alpha_tot.stat, variable == "Hill–Shannon" & Group == "Bacteria")$value_0,
            subset(alpha_tot.stat, variable == "Hill–Shannon" & Group == "Bacteria")$value_6, alternative = "less", paired = TRUE)
wilcox.test(subset(alpha_tot.stat, variable == "Hill–Simpson" & Group == "Bacteria")$value_0,
            subset(alpha_tot.stat, variable == "Hill–Simpson" & Group == "Bacteria")$value_6, alternative = "less", paired = TRUE)

## Fungi
wilcox.test(subset(alpha_tot.stat, variable == "Richness" & Group == "Fungi")$value_0,
            subset(alpha_tot.stat, variable == "Richness" & Group == "Fungi")$value_6, alternative = "less", paired = TRUE)
wilcox.test(subset(alpha_tot.stat, variable == "Hill–Shannon" & Group == "Fungi")$value_0,
            subset(alpha_tot.stat, variable == "Hill–Shannon" & Group == "Fungi")$value_6, alternative = "less", paired = TRUE)
wilcox.test(subset(alpha_tot.stat, variable == "Hill–Simpson" & Group == "Fungi")$value_0,
            subset(alpha_tot.stat, variable == "Hill–Simpson" & Group == "Fungi")$value_6, alternative = "less", paired = TRUE)

#
################################### EXPORT FIGURE                ####

## BETA DIVERSITY
gg.div_beta <- ggplot(nMDS_prok.plot) + 
  geom_line(aes(x = MDS1, y = MDS2, fill = Material, group = ID_block), show.legend = FALSE) +
  geom_point(aes(x = MDS1, y = MDS2, fill = Material, shape = Biochar, size = Biochar, color = sample_time), 
             stroke = 1.5) +
  scale_shape_manual(values = c(21, 23)) +
  scale_size_manual(values = c(4, 3)) +
  scale_fill_manual(values = col_mat) +
  scale_color_manual(values = c("#000000", "#999999")) +
  annotate(geom = "text", x = Inf, y = Inf, label = paste("stress =", round(nMDS_prok$stress, 3)), size = 6, vjust = 1.2, hjust = 1.05) +
  theme_bw() +
  guides(shape = guide_legend(override.aes = list(size = c(4, 3))),
         fill = guide_legend(override.aes = list(shape = 21, size = 4, 
                                                 fill = c("#2ac219", "#facb7a", "#613b25", "#ae0e36"))),
         color = guide_legend(title = "Time", override.aes = list(shape = 21, size = 4))) +
  theme(legend.position = "right",
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"))

gg.div_beta.leg <- get_legend(gg.div_beta)

gg.beta <- plot_grid(gg.nmds.bray_prok, gg.div_beta.leg, gg.nmds.bray_fun, nrow = 1, rel_widths = c(1, 0.25, 1),
                     labels = c("A", "", "B"), label_size = 18)
gg.beta


## ALPHA DIVERSITY
gg.div_alpha.leg <- ggplot(data = alpha_prok.plot) +
  geom_point(aes(x = Treatment, y = value, fill = Treatment, color = sample_time), 
             size = 3, shape = 21, position = position_jitterdodge(), show.legend = FALSE) +  
  geom_boxplot(aes(x = Treatment, y = value, color = sample_time), 
               position = position_dodge(), alpha = 0.75, outlier.shape = NA) +
  scale_fill_manual(values = col_trt) +
  scale_color_manual(values = col_t) +
  facet_wrap( ~ variable, scales = "free_y") +
  ylab("") + xlab("") + 
  guides(color = guide_legend(title = "Time", override.aes = list(linewidth = 1.5)),
         fill = guide_legend(nrow = 1)) +
  theme_bw() +
  theme(legend.position = "bottom",
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"))

gg.div_alpha.leg <- get_plot_component(gg.div_alpha.leg, "guide-box", return_all = TRUE) 

gg.alpha <- plot_grid(plot_grid(gg.alpha_prok, gg.alpha_fun, labels = c("C", "D"), label_size = 18), 
                      gg.div_alpha.leg[[3]], ncol = 1, rel_heights = c(1, 0.1))
gg.alpha


## FIGURE 2

gg.figure2 <- plot_grid(gg.beta, gg.alpha, ncol = 1, rel_heights = c(1,1.25))
gg.figure2

ggsave("../Manuscript/Figures/Figure_2.png", gg.figure2, bg = "white", width = 10, height = 10)

#
################################################################################ SUPPLEMENTARY FIGURE S4                     ####
################################### TAXONOMIC EXPLORATION        ####
#### BACTERIA               ####

asv.t_prok <- Reduce("+", asv.r_prok) / length(asv.r_prok)
asv.t_prok <-  apply(asv.t_prok, 1, function(x) x/sum(x))

asv.tax_prok <- merge(tax.ud_prok, asv.t_prok, by.x = "Seq_ID", by.y = "row.names")

asv.tax_prok.plot <- melt(asv.tax_prok) 

## PHYLUM

asv.tax_prok.plot <- aggregate(value ~ variable + Phylum, asv.tax_prok.plot, sum)

## Assess minimum abundance for further thresholds
max.ab_phy <- aggregate(value ~ Phylum, subset(asv.tax_prok.plot, value > 0), max)
asv.tax_prok.plot[!asv.tax_prok.plot$Phylum %in% max.ab_phy[max.ab_phy$value >= 2.5e-2,]$Phylum, "Phylum"] <- "Other"

asv.tax_prok.plot <- aggregate(value ~ variable + Phylum, asv.tax_prok.plot, sum)
asv.tax_prok.plot <- merge(sample_df, asv.tax_prok.plot, by.x = "ID_well", by.y = "variable")

asv.tax_prok.plot$Block <- gsub("[0-9]", "", asv.tax_prok.plot$ID_block)

orderP <- levels(factor(asv.tax_prok.plot$Phylum))
orderP <- orderP[! orderP %in% c("Other", "Unidentified")]
orderP <- append(orderP, c("Other"))

## FULL TAXONOMY
gg.tax_prok <- ggplot(asv.tax_prok.plot, aes(x = Block, y = value, fill = factor(Phylum, levels = orderP))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Phylum", values = c("#f05b5b", "#ccebc5", "#ffff33", "#9b38c2", "#e6ab02", "#1ee85b", "#80b1d3", "#ed244c", 
                                                "#ffaf0f", "#e6f5c9", "#9c9cdb", "#fff2ae", "#bf5b17")) +
  guides(fill = guide_legend(nrow = 3)) + 
  facet_grid(sample_time ~ Treatment, scales = "free") +
  ylab("Abundance") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 16, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 15, color = "black"))

gg.tax_prok

#
#### FUNGI                  ####

asv.t_fun <- Reduce("+", asv.r_fun) / length(asv.r_fun)
asv.t_fun <-  apply(asv.t_fun, 1, function(x) x/sum(x))

asv.tax_fun <- merge(tax_fun, asv.t_fun, by.x = "Seq_ID", by.y = "row.names")

asv.tax_fun.plot <- melt(asv.tax_fun) 

## CLASS

asv.tax_fun.plot <- aggregate(value ~ variable + Class, asv.tax_fun.plot, sum)
asv.tax_fun.plot$Class <- gsub("[a-z]__", "", asv.tax_fun.plot$Class)

## Assess minimum abundance for further thresholds
max.cl_phy <- aggregate(value ~ Class, subset(asv.tax_fun.plot, value > 0), max)
asv.tax_fun.plot[!asv.tax_fun.plot$Class %in% max.cl_phy[max.cl_phy$value >= 5e-2,]$Class, "Class"] <- "Other"

asv.tax_fun.plot <- aggregate(value ~ variable + Class, asv.tax_fun.plot, sum)

asv.tax_fun.plot <- merge(sample_df, asv.tax_fun.plot, by.x = "ID_well", by.y = "variable")

asv.tax_fun.plot$Block <- gsub("[0-9]", "", asv.tax_fun.plot$ID_block)

orderC <- levels(factor(asv.tax_fun.plot$Class))
orderC <- orderC[! orderC %in% c("Other", "Unidentified")]
orderC <- append(orderC, c("Other", "Unidentified"))

## FULL TAXONOMY
gg.tax_fun <- gg.tax_fun <- ggplot(asv.tax_fun.plot, 
                                   aes(x = Block, y = value, fill = factor(Class, levels = orderC))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Class", values = c("#a6d483", "#e84f61", "#ffff33", "#e6ab02", "#ecf299", "#65c4c9", "#f5d0bc", 
                                               "#bebada", "#feffba", "#caede8", "#bf5b17", "#6a3d9a")) +
  guides(fill = guide_legend(nrow = 3)) + 
  facet_grid(sample_time ~ Treatment, scales = "free") +
  ylab("Abundance") +
  theme_bw() + 
  theme(title = element_text(size = 16, color = "black"),
        legend.position = "bottom", 
        legend.text.align = 0,
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 14, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 15, color = "black"),
        legend.title = element_text(size = 16, color = "black"))

gg.tax_fun

#
################################### EXPORT FIGURE                ####

## TAXONOMIC EXPLORATION
gg.figureS4 <- plot_grid(gg.tax_prok, gg.tax_fun, nrow = 2,
                         labels = c("A", "B"), label_size = 18)

gg.figureS4

ggsave("../Manuscript/Figures/Supp_Figure_S4.png", gg.figureS4, bg = "white", width = 10, height = 8)

#
################################################################################ SUPPLEMENTARY FIGURE S5      ####
#### ALPHA INITIAL vs FINAL ####

## Bacteria
wilcox.test(subset(alpha_tot.stat, variable == "Richness" & Group == "Bacteria")$value_0,
            subset(alpha_tot.stat, variable == "Richness" & Group == "Bacteria")$value_6, alternative = "less", paired = TRUE)
wilcox.test(subset(alpha_tot.stat, variable == "Hill–Shannon" & Group == "Bacteria")$value_0,
            subset(alpha_tot.stat, variable == "Hill–Shannon" & Group == "Bacteria")$value_6, alternative = "less", paired = TRUE)
wilcox.test(subset(alpha_tot.stat, variable == "Hill–Simpson" & Group == "Bacteria")$value_0,
            subset(alpha_tot.stat, variable == "Hill–Simpson" & Group == "Bacteria")$value_6, alternative = "less", paired = TRUE)

## Fungi
wilcox.test(subset(alpha_tot.stat, variable == "Richness" & Group == "Fungi")$value_0,
            subset(alpha_tot.stat, variable == "Richness" & Group == "Fungi")$value_6, alternative = "less", paired = TRUE)
wilcox.test(subset(alpha_tot.stat, variable == "Hill–Shannon" & Group == "Fungi")$value_0,
            subset(alpha_tot.stat, variable == "Hill–Shannon" & Group == "Fungi")$value_6, alternative = "less", paired = TRUE)
wilcox.test(subset(alpha_tot.stat, variable == "Hill–Simpson" & Group == "Fungi")$value_0,
            subset(alpha_tot.stat, variable == "Hill–Simpson" & Group == "Fungi")$value_6, alternative = "less", paired = TRUE)


#
################################### EXPORT FIGURE                ####

gg.alpha_line <- ggplot(data = alpha_tot.plot) +
  geom_point(aes(x = sample_time, y = value, color = Material, shape = Biochar), size = 3) +
  geom_line(aes(x = sample_time, y = value, color = Material, group = ID_block), linewidth = 1) +
  scale_color_manual(values = col_mat) +
  facet_wrap(Group ~ variable, scales = "free_y") +
  ylab("") + xlab("") + 
  theme_bw() +
  theme(legend.position = "right",
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black"),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"))

gg.figureS5 <- gg.alpha_line

ggsave("../Manuscript/Figures/Supp_Figure_S5.png", gg.figureS5, bg = "white", width = 12, height = 10)

################################################################################ SUPPLEMENTARY FIGURE S6      ####
################################### NETWORK ANALYSIS - GLOBAL    ####
#### BACTERIA               ####

genus_prok <- asv.t_prok
genus_prok <- aggregate(genus_prok ~ Genus, tax.ud_prok, sum)
row.names(genus_prok) <- genus_prok$Genus
genus_prok <- genus_prok[,-1]

## Occurrence filter (ASVs present in more than 10 and less than 90% of samples)
core.g_prok <- genus_prok[rowSums(genus_prok > 0) > 0.33*ncol(genus_prok) & rowSums(genus_prok > 0), ]

mean(colSums(core.g_prok))
sd(colSums(core.g_prok))
nrow(core.g_prok)/nrow(genus_prok)

# 40.06% of ASVs (498/1243 total) meaning 96.67 ± 1.73% of reads

## CO-OCCURRENCE NETWORKS
genus.net_prok <- core.g_prok
genus.net_prok[genus.net_prok > 0] <- 1

veech_prok <- veech_prob(t(genus.net_prok), type = "both")

net_prok <- graph_from_data_frame(veech_prok[veech_prok$type == "Co-occurrence", 1:2], directed = FALSE)
cw_prok <- cluster_walktrap(net_prok, steps = 20)

mod_prok <- cbind.data.frame(names = cw_prok$names, module = cw_prok$membership)
max(mod_prok$module)

# Layout
set.seed(1)
l <- layout_with_fr(net_prok, grid = "nogrid")
l <- norm_coords(l, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
l_df <- as.data.frame(l)  ## convert the layout to a data.frame
l_df$names <- names(V(net_prok))  ## add in the species codes

## NODE TABLE
node_prok <- merge(mod_prok, l_df, by = "names")
node_prok <- node_prok[order(node_prok$module, decreasing = FALSE),]

## EDGE TABLE
edge_prok <- veech_prok
edge_prok <- merge(veech_prok, node_prok[,-2], by.x = "sp2", by.y = "names")
edge_prok <- merge(edge_prok, node_prok[,-2], by.x = "sp1", by.y = "names")
colnames(edge_prok) <- c("sp1", "sp2", "p.value", "type", "sp2.x", "sp2.y", "sp1.x", "sp1.y")

## GET MODULES WITH MORE THAN 10 NODES
aggregate(names ~ module, mod_prok, length)

node_prok$module <- ifelse(node_prok$module == 1, "Module Bac-1", 
                           ifelse(node_prok$module == 4, "Module Bac-2", ifelse(node_prok$module == 5, "Module Bac-3", NA)))

## DRAW NETWORK
gg.net_prok <- ggplot() +
  geom_curve(data = edge_prok, aes(x = sp1.x, xend = sp2.x, y = sp1.y, yend = sp2.y), 
             color = "gray70", linewidth = 0.75, alpha = 0.8) +
  geom_point(data = node_prok, aes(x = V1, y = V2, fill = factor(module)), 
             size = 5, shape = 21, stroke = 1.25, color = "black", show.legend = FALSE) +
  scale_shape_manual(values = c(22, 21)) +
  scale_fill_manual(values = c("#910a25", "#00bfff", "#ff8247"), na.value = "white") +
  theme_void() 

gg.net_prok

## NETWORK PROPERTIES
length(V(net_prok))                                     # Number of nodes = 429
length(E(net_prok))                                     # Number of nodes = 17196
edge_density(net_prok)                                  # Edge density    = 0.187308  Edges/((Nodes^2-Nodes)/2)
transitivity(net_prok, type = "average")                # Clustering      = 0.5684415
max(membership(cw_prok))                                # Module number   = 12
modularity(cw_prok)                                     # Modularity      = 0.2710322
mean_distance(net_prok, directed = FALSE)               # Average path    = 2.115592
assortativity_degree(net_prok, directed = "FALSE")      # Assortativity   = 0.4410126

#
#### FUNGI                  ####

genus_fun <- asv.t_fun
genus_fun <- aggregate(genus_fun ~ Genus, tax.ud_fun, sum)
row.names(genus_fun) <- genus_fun$Genus
genus_fun <- genus_fun[,-1]

## Occurrence filter (ASVs present in more than 10 and less than 90% of samples)
core.g_fun <- genus_fun[rowSums(genus_fun > 0) > 0.33*ncol(genus_fun) & rowSums(genus_fun > 0), ]

mean(colSums(core.g_fun))
sd(colSums(core.g_fun))
nrow(core.g_fun)/nrow(genus_fun)

# 22.83% of ASVs (174/762 total) meaning 98.01 ± 2.09% of reads

## CO-OCCURRENCE NETWORKS
genus.net_fun <- core.g_fun
genus.net_fun[genus.net_fun > 0] <- 1

veech_fun <- veech_prob(t(genus.net_fun), type = "both")

net_fun <- graph_from_data_frame(veech_fun[veech_fun$type == "Co-occurrence", 1:2], directed = FALSE)
cw_fun <- cluster_walktrap(net_fun, steps = 20)

mod_fun <- cbind.data.frame(names = cw_fun$names, module = cw_fun$membership)
max(mod_fun$module)

# Layout
set.seed(1)
l <- layout_with_fr(net_fun, grid = "nogrid")
l <- norm_coords(l, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
l_df <- as.data.frame(l)  ## convert the layout to a data.frame
l_df$names <- names(V(net_fun))  ## add in the species codes

## NODE TABLE
node_fun <- merge(mod_fun, l_df, by = "names")
node_fun <- node_fun[order(node_fun$module, decreasing = FALSE),]

## EDGE TABLE
edge_fun <- veech_fun
edge_fun <- merge(veech_fun, node_fun[,-2], by.x = "sp2", by.y = "names")
edge_fun <- merge(edge_fun, node_fun[,-2], by.x = "sp1", by.y = "names")
colnames(edge_fun) <- c("sp1", "sp2", "p.value", "type", "sp2.x", "sp2.y", "sp1.x", "sp1.y")

## GET MODULES WITH MORE THAN 10 NODES
aggregate(names ~ module, mod_fun, length)

node_fun$module <- ifelse(node_fun$module == 1, "Module Fun-3", 
                          ifelse(node_fun$module == 2, "Module Fun-2", ifelse(node_fun$module == 3, "Module Fun-1", NA)))

## DRAW NETWORK
gg.net_fun <- ggplot() +
  geom_curve(data = edge_fun, aes(x = sp1.x, xend = sp2.x, y = sp1.y, yend = sp2.y), 
             color = "gray70", linewidth = 0.75, alpha = 0.8) +
  geom_point(data = node_fun, aes(x = V1, y = V2, fill = factor(module)), 
             size = 5, shape = 22, stroke = 1.25, color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("#910a25", "#00bfff", "#ff8247")) +
  theme_void() 

gg.net_fun

## NETWORK PROPERTIES
length(V(net_fun))                                     # Number of nodes = 158
length(E(net_fun))                                     # Number of nodes = 2824
edge_density(net_fun)                                  # Edge density    = 0.2276868  Edges/((Nodes^2-Nodes)/2)
transitivity(net_fun, type = "average")                # Clustering      = 0.6157927
max(membership(cw_fun))                                # Module number   = 3
modularity(cw_fun)                                     # Modularity      = 0.3283893
mean_distance(net_fun, directed = FALSE)               # Average path    = 2.079739
assortativity_degree(net_fun, directed = "FALSE")      # Assortativity   = 0.358618

#
#### BIPARTITE              ####

core_bi <- rbind(core.g_prok, core.g_fun[,colnames(core.g_prok)])

## CO-OCCURRENCE NETWORKS
genus.net_bi <- core_bi
genus.net_bi[genus.net_bi > 0] <- 1

veech_bi <- veech_prob(t(genus.net_bi), type = "both")

veech_bi$sp1.k <- ifelse(veech_bi$sp1 %in% tax.ud_prok$Genus, "Bacteria", ifelse(veech_bi$sp1 %in% tax.ud_fun$Genus, "Fungi", ""))
veech_bi$sp2.k <- ifelse(veech_bi$sp2 %in% tax.ud_prok$Genus, "Bacteria", ifelse(veech_bi$sp2 %in% tax.ud_fun$Genus, "Fungi", ""))

## EDGE TABLE
edge_bi <- subset(veech_bi, type == "Co-occurrence")
edge_bi <- edge_bi[edge_bi$sp1.k != edge_bi$sp2.k,]

mod_prok$module <- ifelse(mod_prok$module %in% c(1,4,5), paste(mod_prok$module, "_bac", sep = ""), NA)
mod_fun$module <- ifelse(mod_fun$module %in% c(1,2,3), paste(mod_fun$module, "_fun", sep = ""), NA)

mod_bi <- rbind.data.frame(mod_prok, mod_fun)
mod_bi <- mod_bi[mod_bi$names %in% edge_bi$sp1 | mod_bi$names %in% edge_bi$sp2,]

net_bi <- graph_from_data_frame(edge_bi[, 1:2], directed = FALSE)

# Layout
set.seed(1)
l <- layout_with_fr(net_bi, grid = "nogrid")
l <- norm_coords(l, ymin = -1, ymax = 1, xmin = -1, xmax = 1)
l_df <- as.data.frame(l)  ## convert the layout to a data.frame
l_df$names <- names(V(net_bi))  ## add in the species codes

## NODE TABLE
node_bi <- merge(mod_bi, l_df, by = "names", all = TRUE)
node_bi <- node_bi[order(node_bi$module, decreasing = FALSE),]
node_bi$Kingdom <- ifelse(node_bi$names %in% tax.ud_prok$Genus, "Bacteria", "Fungi")

node_bi <- node_bi[!node_bi$names %in% c("Lacunisphaera", "Phialemonium", "Unidentified 0319-6G20", "Unidentified IMCC26256"),]

## EDGE TABLE
edge_bi <- merge(edge_bi, node_bi[,c(1,3,4)], by.x = "sp2", by.y = "names")
edge_bi <- merge(edge_bi, node_bi[,c(1,3,4)], by.x = "sp1", by.y = "names")
colnames(edge_bi) <- c("sp1", "sp2", "p.value", "type", "sp1.k", "sp2.k", "sp2.x", "sp2.y", "sp1.x", "sp1.y")

## DRAW NETWORK
gg.net_bi <- ggplot() +
  geom_curve(data = edge_bi, aes(x = sp1.x, xend = sp2.x, y = sp1.y, yend = sp2.y), 
             color = "gray70", linewidth = 0.75, alpha = 0.8) +
  geom_point(data = node_bi, aes(x = V1, y = V2, color = factor(module), fill = factor(module)), 
             size = 5, shape = ifelse(node_bi$Kingdom == "Bacteria", 21, 22), color = "black", stroke = 1.25, show.legend = FALSE) +
  scale_fill_manual(values = c("#910a25", "#ff8247", "#00bfff", "#910a25", "#00bfff", "#ff8247"), na.value = "white") +
  guides(fill = guide_legend(override.aes = list(shape = 21, size = 3, color = NA))) +
  theme_void() 

gg.net_bi

#
################################### EXPORT FIGURE                ####

gg.figureS6 <- plot_grid(plot_grid(gg.net_prok, gg.net_fun, nrow = 1, labels = c("A", "B"), label_size = 18), 
                         ncol = 1, gg.net_bi, labels = c("", "C"), label_size = 18)

ggsave("../Manuscript/Figures/Supp_Figure_S6.png", gg.figureS6, bg = "white", width = 12, height = 12)

#
################################################################################ FIGURE 3                     ####
################################### NETWORK ANALYSIS - LOCAL     ####
#### BACTERIA               ####

## MODULE COMPLETENESS

mod.comp_prok <- cbind.data.frame(Sample = colnames(genus.net_prok),
                                  `Module Bac-1` = colSums(genus.net_prok[subset(node_prok, module == "Module Bac-1")$names,]) / 
                                    colSums(genus.net_prok),
                                  `Module Bac-2` = colSums(genus.net_prok[subset(node_prok, module == "Module Bac-2")$names,]) / 
                                    colSums(genus.net_prok),
                                  `Module Bac-3` = colSums(genus.net_prok[subset(node_prok, module == "Module Bac-3")$names,]) / 
                                    colSums(genus.net_prok))

mod.comp_prok <- merge(mod.comp_prok, sample_df[,1:7], by.x = "Sample", by.y = "ID_well")

mod.comp_prok.plot <- melt(mod.comp_prok)

## BOXPLOT
gg.mod_prok <- ggplot(data = mod.comp_prok.plot) +
  geom_point(aes(x = Treatment, y = value, fill = Treatment, color = sample_time), 
             size = 2, shape = 21, position = position_jitterdodge(), show.legend = FALSE) +  
  geom_boxplot(aes(x = Treatment, y = value, fill = Treatment, color = sample_time), 
               position = position_dodge(), alpha = 0.75, outlier.shape = NA) +
  scale_fill_manual(values = col_trt) +
  scale_color_manual(values = col_t) +
  facet_wrap( ~ variable, scales = "free_y") +
  ylab("") + xlab("") + 
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"))

gg.mod_prok

#
#### FUNGI                  ####

## MODULE COMPLETENESS

mod.comp_fun <- cbind.data.frame(Sample = colnames(genus.net_fun),
                                 `Module Fun-3` = colSums(genus.net_fun[subset(node_fun, module == "Module Fun-3")$names,]) / 
                                   colSums(genus.net_fun),
                                 `Module Fun-2` = colSums(genus.net_fun[subset(node_fun, module == "Module Fun-2")$names,]) / 
                                   colSums(genus.net_fun),
                                 `Module Fun-1` = colSums(genus.net_fun[subset(node_fun, module == "Module Fun-1")$names,]) / 
                                   colSums(genus.net_fun))

mod.comp_fun <- merge(mod.comp_fun, sample_df[,1:7], by.x = "Sample", by.y = "ID_well")

mod.comp_fun.plot <- melt(mod.comp_fun)
mod.comp_fun.plot$variable <- factor(mod.comp_fun.plot$variable, levels = c("Module Fun-1", "Module Fun-2", "Module Fun-3"))

## BOXPLOT
gg.mod_fun <- ggplot(data = mod.comp_fun.plot) +
  geom_point(aes(x = Treatment, y = value, fill = Treatment, color = sample_time), 
             size = 2, shape = 21, position = position_jitterdodge(), show.legend = FALSE) +  
  geom_boxplot(aes(x = Treatment, y = value, fill = Treatment, color = sample_time), 
               position = position_dodge(), alpha = 0.75, outlier.shape = NA) +
  scale_fill_manual(values = col_trt) +
  scale_color_manual(values = col_t) +
  facet_wrap( ~ variable, scales = "free_y") +
  ylab("") + xlab("") + 
  theme_bw() +
  theme(legend.position = "none",
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 90, vjust = 0.5, hjust = 1),
        strip.text = element_text(size = 15, color = "black"),
        strip.background = element_rect(fill = "white"))

gg.mod_fun

#
################################### CORRELATIONS                 ####
#### HEATMAP                ####

mod.comp_stat <- merge(mod.comp_prok[,1:4], mod.comp_fun, by = "Sample")

mod.comp_stat.6 <- merge(mod.comp_stat, sample_df[,c(1,8:15)], by.x = "Sample", by.y = "ID_well")
mod.comp_stat.6 <- mod.comp_stat.6[complete.cases(mod.comp_stat.6),]

cor_df <- rcorr(as.matrix(mod.comp_stat.6[,c(2:7,14:21)]), type = "spearman")

cor.p_df <- as.matrix(cor_df$P)
cor.p_df[is.na(cor.p_df)] <- 1
cor.p_df <- apply(cor.p_df, 1, function(x) p.adjust(x, method = "fdr"))

cor.r_df <- as.matrix(cor_df$r)
diag(cor.r_df) <- NA

cor.p_df <- cor.p_df[c(1,2,3,6,5,4),7:14]
cor.r_df <- cor.r_df[c(1,2,3,6,5,4),7:14]

cor.r_df.plot <- merge(melt(cor.r_df), melt(cor.p_df), by = c("Var1", "Var2"))
cor.r_df.plot$signif <- ifelse(cor.r_df.plot$value.y < 0.001, "***", 
                               ifelse(cor.r_df.plot$value.y < 0.01, "**", ifelse(cor.r_df.plot$value.y < 0.05, "*", "")))

gg.heat <- ggplot(cor.r_df.plot, aes(x = Var1, y = Var2, fill = value.x)) + 
  geom_tile(color = "gray50") + 
  scale_fill_gradient2(name = "Spearman\nCorrelation", high = "#c0272f", mid = "#ffffff", low = "#3f4eaa", midpoint = 0) + 
  geom_text(aes(x = Var1, y = Var2, label = signif), color = "black", size = 4) +
  coord_fixed() +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20)) +
  scale_y_discrete(limits=rev)+
  theme_void() + 
  theme(legend.position = "right",
        legend.title = element_text(size = 17, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"))

gg.heat

#
################################### EXPORT FIGURE                ####

gg.figure3 <- plot_grid(plot_grid(gg.mod_prok, gg.mod_fun, ncol = 1), gg.heat, labels = c("A", "B"), label_size = 18)
gg.figure3

ggsave("../Manuscript/Figures/Figure_3.png", gg.figure3, bg = "white", width = 12, height = 9)

#
################################################################################ SUPPLEMENTARY FIGURE S7      ####
################################### MODULE TAXONOMY              ####
#### BACTERIA               ####

tax.mod_prok <- merge(node_prok[,1:2], genus_prok, by.x = "names", by.y = "row.names")

## Module Bac-2
tax.mod_prok <- subset(tax.mod_prok, module == "Module Bac-2")
tax.mod_prok[,3:71] <- apply(tax.mod_prok[,3:71], 2, function(x) x/sum(x))

tax.mod_prok.plot <- merge(unique(tax.ud_prok[,1:6]), melt(tax.mod_prok, id.vars = c("names", "module")), by.x = "Genus", by.y = "names")

## PHYLUM

mod.phy_prok.plot <- aggregate(value ~ variable + Phylum + module, tax.mod_prok.plot, sum)

## Assess minimum abundance for further thresholds
mod.phy_prok.plot[!mod.phy_prok.plot$Phylum %in% unique(asv.tax_prok.plot$Phylum), "Phylum"] <- "Other"
mod.phy_prok.plot <- aggregate(value ~ variable + Phylum + module, mod.phy_prok.plot, sum)

mod.phy_prok.plot <- merge(sample_df, mod.phy_prok.plot, by.x = "ID_well", by.y = "variable")

mod.phy_prok.plot$Block <- gsub("[0-9]", "", mod.phy_prok.plot$ID_block)

mod.sub.phy_prok.plot <- subset(mod.phy_prok.plot, sample_time == "Final")

comp.phy_prok.plot <- rbind(cbind.data.frame(Block = mod.sub.phy_prok.plot$Block, sample_time = "Module Bac-2", mod.sub.phy_prok.plot[,c(7,19,17)]), 
                            cbind.data.frame(subset(asv.tax_prok.plot, sample_time == "Final")[,c(19,4,7,18,17)]))

orderP <- levels(factor(comp.phy_prok.plot$Phylum))
orderP <- orderP[! orderP %in% c("Other", "Unidentified")]
orderP <- append(orderP, c("Other"))

gg.tax.mod_prok <- ggplot(comp.phy_prok.plot, aes(x = Block, y = value, fill = factor(Phylum, levels = orderP))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Phylum", values = c("#f05b5b", "#ccebc5", "#ffff33", "#9b38c2", "#e6ab02", "#1ee85b", "#80b1d3", "#ed244c", 
                                                "#ffaf0f", "#e6f5c9", "#9c9cdb", "#fff2ae", "#bf5b17")) +
  guides(fill = guide_legend(nrow = 3)) + 
  facet_grid(sample_time ~ Treatment, scales = "free") +
  ylab("Abundance") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 16, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 15, color = "black"))

gg.tax.mod_prok

#
#### FUNGI                  ####

tax.mod_fun <- merge(node_fun[,1:2], genus_fun, by.x = "names", by.y = "row.names")

## Module Fun-2
tax.mod_fun <- subset(tax.mod_fun, module == "Module Fun-2")
tax.mod_fun[,3:72] <- apply(tax.mod_fun[,3:72], 2, function(x) x/sum(x))

tax.mod_fun.plot <- merge(unique(tax.ud_fun[,1:6]), melt(tax.mod_fun, id.vars = c("names", "module")), by.x = "Genus", by.y = "names")

## CLASS

mod.phy_fun.plot <- aggregate(value ~ variable + Class + module, tax.mod_fun.plot, sum)

## Assess minimum abundance for further thresholds
mod.phy_fun.plot[!mod.phy_fun.plot$Class %in% unique(asv.tax_fun.plot$Class), "Class"] <- "Other"
mod.phy_fun.plot <- aggregate(value ~ variable + Class + module, mod.phy_fun.plot, sum)

mod.phy_fun.plot <- merge(sample_df, mod.phy_fun.plot, by.x = "ID_well", by.y = "variable")

mod.phy_fun.plot$Block <- gsub("[0-9]", "", mod.phy_fun.plot$ID_block)

mod.sub.phy_fun.plot <- subset(mod.phy_fun.plot, sample_time == "Final")

comp.phy_fun.plot <- rbind(cbind.data.frame(Block = mod.sub.phy_fun.plot$Block, sample_time = "Module Fun-2", mod.sub.phy_fun.plot[,c(7,19,17)]), 
                            cbind.data.frame(subset(asv.tax_fun.plot, sample_time == "Final")[,c(19,4,7,18,17)]))

orderC <- levels(factor(comp.phy_fun.plot$Class))
orderC <- orderC[! orderC %in% c("Other", "Unidentified")]
orderC <- append(orderC, c("Other", "Unidentified"))

gg.tax.mod_fun <- ggplot(comp.phy_fun.plot, aes(x = Block, y = value, fill = factor(Class, levels = orderC))) + 
  geom_bar(stat = "identity", position = "stack") + 
  scale_fill_manual(name = "Class", values = c("#a6d483", "#e84f61", "#ffff33", "#e6ab02", "#ecf299", "#65c4c9", "#f5d0bc", 
                                               "#bebada", "#feffba", "#caede8", "#bf5b17", "#6a3d9a")) +
  guides(fill = guide_legend(nrow = 3)) + 
  facet_grid(sample_time ~ Treatment, scales = "free") +
  ylab("Abundance") +
  theme_bw() + 
  theme(legend.position = "bottom", 
        legend.text.align = 0,
        axis.title.y = element_text(size = 16, color = "black"),
        axis.text.y = element_text(size = 15, color = "black"),
        axis.title.x = element_blank(),
        axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        strip.text = element_text(size = 16, color = "black"),
        strip.background = element_rect(fill = "white"),
        legend.text = element_text(size = 14, color = "black"),
        legend.title = element_text(size = 15, color = "black"))

gg.tax.mod_fun

#
################################### EXPORT FIGURE                ####

gg.figureS7 <- plot_grid(gg.tax.mod_prok, gg.tax.mod_fun, ncol = 1, labels = c("A", "B"), label_size = 18)
gg.figureS7

ggsave("../Manuscript/Figures/Supp_Figure_S7.png", gg.figureS7, bg = "white", width = 12, height = 11)

#
################################################################################ FIGURE 4                     ####
################################### SEM                          ####
#### SEM                    ####
## ALPHA DIVERSITY

sem_df <- merge(alpha_prok[,1:4], alpha_fun, by = "ID_well", all = TRUE)


## MODULES
sem_df <- merge(mod.comp_prok[,1:7], sem_df, by.x = "Sample", by.y = "ID_well", all = TRUE)
sem_df <- merge(sem_df, mod.comp_fun[,1:4], by = "Sample", all = TRUE)

colnames(sem_df) <- gsub("\\.x", "_bac", colnames(sem_df))
colnames(sem_df) <- gsub("\\.y", "_fun", colnames(sem_df))


sem_df <- sem_df[complete.cases(sem_df),]

colnames(sem_df) <- gsub("Module ", "", colnames(sem_df))
colnames(sem_df) <- gsub("-", "_", colnames(sem_df))

sem_sem <- psem(
  
  lm(Yield ~ `Bac_2` + `Fun_2` + pH + EC + C + N + P + K, data = sem_df),
  
  lm(`Bac_2` ~ Material + Biochar + pH + EC + C + N + P + K, data = sem_df),
  lm(`Fun_2` ~ Material + Biochar + pH + EC + C + N + P + K, data = sem_df),
  Bac_2 %~~% Fun_2,
  
  lm(pH ~ Material + Biochar + C + N + P + K, data = sem_df),
  lm(EC ~ Material + Biochar + C + N + P + K, data = sem_df),
  
  lm(C ~ Material + Biochar, data = sem_df),
  
  lm(N ~ Material + Biochar + C, data = sem_df),
  lm(P ~ Material + Biochar + C + N, data = sem_df),
  lm(K ~ Material + Biochar + C, data = sem_df),
  
  data = sem_df
)

summary(sem_sem)

coefs_sem <- coefs(sem_sem)
rsq_sem <- rsquared(sem_sem)

std.est <- NULL
for (res in 1:nrow(coefs_sem)) {
  if (!grepl("Material|Biochar", coefs_sem[res,"Predictor"])) {
    
    beta <- as.numeric(coefs_sem[res, "Estimate"]) * (sd(as.numeric(sem_df[,coefs_sem[res,"Predictor"]])) /
                                                        as.numeric(sd(sem_df[,coefs_sem[res,"Response"]])))
    std.est <- rbind(std.est, cbind.data.frame(Response = coefs_sem[res,"Response"], Predictor = coefs_sem[res,"Predictor"], std.est = beta))
    
    
  }
}

coefs_sem$order <- as.numeric(row.names(coefs_sem))
coefs_sem <- merge(coefs_sem, std.est, by = c("Response", "Predictor"), all = TRUE)
coefs_sem <- coefs_sem[order(coefs_sem$order), -10]
colnames(coefs_sem)[9] <- ""

#
################################### EXPORT TABLES                ####

write.table(coefs_sem, "Figures/Figure_4_coefs.txt", sep = "\t", quote = FALSE, row.names = FALSE)
write.table(rsq_sem, "Figures/Figure_4_rsq.txt", sep = "\t", quote = FALSE, row.names = FALSE)

#
################################################################################ SUPPLEMENTARY FIGURE S8      ####
#### HEATMAP                ####

cor.sem_df <- rcorr(as.matrix(sem_df[,c(27,3,29,20:25)]), type = "spearman")

cor.sem.p_df <- as.matrix(cor.sem_df$P)
cor.sem.p_df[is.na(cor.sem.p_df)] <- 1
cor.sem.p_df <- apply(cor.sem.p_df, 1, function(x) p.adjust(x, method = "fdr"))

cor.sem.r_df <- as.matrix(cor.sem_df$r)
diag(cor.sem.r_df) <- NA

cor.sem.r_df.plot <- merge(melt(cor.sem.r_df), melt(cor.sem.p_df), by = c("Var1", "Var2"))
cor.sem.r_df.plot$signif <- ifelse(cor.sem.r_df.plot$value.y < 0.001, "***", 
                                   ifelse(cor.sem.r_df.plot$value.y < 0.01, "**", ifelse(cor.sem.r_df.plot$value.y < 0.05, "*", "")))

gg.heat_sem <- ggplot(cor.sem.r_df.plot, aes(x = Var1, y = Var2, fill = value.x)) + 
  geom_tile(color = "gray50") + 
  scale_fill_gradient2(name = "Spearman\nCorrelation", high = "#c0272f", mid = "#ffffff", low = "#3f4eaa", midpoint = 0) + 
  geom_text(aes(x = Var1, y = Var2, label = signif), color = "black", size = 4) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20)) +
  scale_y_discrete(limits=rev)+
  theme_void() + 
  theme(plot.margin = ggplot2::margin(l = 30, b = 15, r = 15, t = 15, unit = "pt"),
        legend.position = "right",
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15, color = "black"),
        axis.text.y = element_text(size = 15, color = "black", hjust = 1))

gg.heat_sem

#
#### RANDOM FOREST          ####

set.seed(1)
rf.yield <- randomForest(formula = Yield ~ `Bac_2` + `Fun_2` + pH + EC + C + N + P + K, data = sem_df,  
                         ntree = 1000, mtry = 1, nodesize = 5, importance = TRUE, keep.forest = TRUE)

rf.yield

# Tuning the model

hyper_grid <- expand.grid(mtry = seq(1, 5, by = 1), 
                          node_size = seq(3, 15, by = 2),
                          OOB_MSE = 0)

for(i in 1:nrow(hyper_grid)) {
  
  set.seed(1)
  
  # train model
  model <- randomForest(formula = Yield ~ `Bac_2` + `Fun_2` + pH + EC + C + N + P + K, data = sem_df, ntree = 1000,
                        keep.forest = FALSE, importance = TRUE, mtry = hyper_grid$mtry[i], nodesize = hyper_grid$node_size[i])
  
  # add OOB error to grid
  hyper_grid$OOB_MSE[i] <- model$mse[length(model$mse)]
  
}

head(hyper_grid[order(hyper_grid$OOB_MSE), ])

set.seed(1)
rf.yield_fit <- randomForest(formula = Yield ~ `Bac_2` + `Fun_2` + pH + EC + C + N + P + K, data = sem_df,  
                             ntree = 1000, mtry = 4, nodesize = 11, importance = TRUE, keep.forest = TRUE)

rf.yield_fit


# Get variable importance from the model fit
ImpData_yield <- as.data.frame(randomForest::importance(rf.yield_fit, scale = TRUE, type = 1))
ImpData_yield$Var.Names <- row.names(ImpData_yield)

ImpData_yield <- ImpData_yield[order(ImpData_yield$`%IncMSE`, decreasing = FALSE),]

ImpData_yield$Var.Names <- factor(ImpData_yield$Var.Names, levels = ImpData_yield$Var.Names)

gg.imp_yield <- ggplot(ImpData_yield) +
  geom_bar(aes(x = Var.Names, y = `%IncMSE`), fill = "#328798", stat = "identity", position = "dodge") +
  xlab("Variable") + ylab("Importance (%IncMSE)") +
  theme_light() +
  coord_flip() +
  theme(plot.margin = ggplot2::margin(l = 30, b = 15, r = 15, t = 15, unit = "pt"),
        axis.title.x = element_text(size = 17, color = "black"),
        axis.title.y = element_blank(),
        axis.text = element_text(size = 15, color = "black"))

gg.imp_yield

#
################################### EXPORT FIGURE                ####

gg.figureS8 <- plot_grid(gg.heat_sem, gg.imp_yield, ncol = 1, rel_heights = c(1,0.75), labels = c("A", "B"), label_size = 18)
gg.figureS8

ggsave("../Manuscript/Figures/Supp_Figure_S8.png", gg.figureS8, bg = "white", width = 12, height = 10)

#
################################################################################ SUPPLEMENTARY FIGURE S9     ####
################################### DIFF ABUND - WILCOXON TEST   ####
#### BACTERIA               ####

## Initial time
mod.t0_prok <- melt(tax.mod_prok[,c("names", colnames(tax.mod_prok)[colnames(tax.mod_prok) %in% 
                                                                      subset(sample_df, sample_time == "Initial")$ID_well])])
mod.t0_prok <- merge(sample_df[,1:7], mod.t0_prok, by.x = "ID_well", by.y = "variable")

## Final time
mod.t6_prok <- melt(tax.mod_prok[,c("names", colnames(tax.mod_prok)[colnames(tax.mod_prok) %in% subset(sample_df, sample_time == "Final")$ID_well])])
mod.t6_prok <- merge(sample_df[,1:7], mod.t6_prok, by.x = "ID_well", by.y = "variable")

mod.t_prok <- merge(mod.t0_prok[,c(2,5:9)], mod.t6_prok[,c(2,5:9)], by = c("ID_block", "Biochar", "Material", "Treatment", "names"))
colnames(mod.t_prok)[6:7] <- c("value_t0", "value_t6")

wilcox_prok <- NULL
for (gen in unique(mod.t_prok$names)) {
  
  gen_prok <- subset(mod.t_prok, names == gen)
  gen_sub <- gen_prok[gen_prok$value_t0 > 0 & gen_prok$value_t6 > 0,]
  
  wilcox_sub <- cbind.data.frame(Genus = gen, mean_t0 = mean(gen_prok$value_t0), mean_t6 = mean(gen_prok$value_t6),
                                 z_value = mean((gen_prok$value_t6 - mean(gen_prok$value_t0)) / sd(gen_prok$value_t0)),
                                 z_sd = sd((gen_prok$value_t6 - mean(gen_prok$value_t0)) / sd(gen_prok$value_t0)),
                                 p.value = wilcox.test(gen_prok$value_t0, gen_prok$value_t6, paired = TRUE)$p.value)
  wilcox_prok <- rbind(wilcox_prok, wilcox_sub)
  
}

wilcox_prok$fdr <- p.adjust(wilcox_prok$p.value, method = "fdr")

wilcox_prok <- wilcox_prok[order(wilcox_prok$z_value, decreasing = TRUE), ]
wilcox_prok <- subset(wilcox_prok, fdr <= 0.05)

gg.wilcox_prok <- ggplot(subset(wilcox_prok, z_value > 0)) +
  geom_bar(aes(x = factor(Genus, wilcox_prok$Genus), y = z_value), stat = "identity",
           fill = "#de3e49") +
  geom_errorbar(aes(x = factor(Genus, wilcox_prok$Genus), ymin = z_value - z_sd, ymax = z_value + z_sd), width = 0.2,
                position = position_dodge(0.05), color = "black") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 15, color = "black"))

gg.wilcox_prok

#
#### FUNGI                  ####

## Initial time
mod.t0_fun <- melt(tax.mod_fun[,c("names", colnames(tax.mod_fun)[colnames(tax.mod_fun) %in% subset(sample_df, sample_time == "Initial")$ID_well])])
mod.t0_fun <- merge(sample_df[,1:7], mod.t0_fun, by.x = "ID_well", by.y = "variable")

## Final time
mod.t6_fun <- melt(tax.mod_fun[,c("names", colnames(tax.mod_fun)[colnames(tax.mod_fun) %in% subset(sample_df, sample_time == "Final")$ID_well])])
mod.t6_fun <- merge(sample_df[,1:7], mod.t6_fun, by.x = "ID_well", by.y = "variable")

mod.t_fun <- merge(mod.t0_fun[,c(2,5:9)], mod.t6_fun[,c(2,5:9)], by = c("ID_block", "Biochar", "Material", "Treatment", "names"))
colnames(mod.t_fun)[6:7] <- c("value_t0", "value_t6")

wilcox_fun <- NULL
for (gen in unique(mod.t_fun$names)) {
  
  gen_fun <- subset(mod.t_fun, names == gen)
  gen_sub <- gen_fun[gen_fun$value_t0 > 0 & gen_fun$value_t6 > 0,]
  
  wilcox_sub <- cbind.data.frame(Genus = gen, mean_t0 = mean(gen_fun$value_t0), mean_t6 = mean(gen_fun$value_t6),
                                 z_value = mean((gen_fun$value_t6 - mean(gen_fun$value_t0)) / sd(gen_fun$value_t0)),
                                 z_sd = sd((gen_fun$value_t6 - mean(gen_fun$value_t0)) / sd(gen_fun$value_t0)),
                                 p.value = wilcox.test(gen_fun$value_t0, gen_fun$value_t6, paired = TRUE)$p.value)
  wilcox_fun <- rbind(wilcox_fun, wilcox_sub)
  
}

wilcox_fun$fdr <- p.adjust(wilcox_fun$p.value, method = "fdr")

wilcox_fun <- wilcox_fun[order(wilcox_fun$z_value, decreasing = TRUE), ]
wilcox_fun <- subset(wilcox_fun, fdr <= 0.05)

gg.wilcox_fun <- ggplot(subset(wilcox_fun, z_value > 0)) +
  geom_bar(aes(x = factor(Genus, wilcox_fun$Genus), y = z_value), stat = "identity",
           fill = "#de3e49") +
  geom_errorbar(aes(x = factor(Genus, wilcox_fun$Genus), ymin = z_value - z_sd, ymax = z_value + z_sd), width = 0.2,
                position = position_dodge(0.05), color = "black") +
  geom_hline(yintercept = 0) +
  theme_bw() +
  theme(axis.title.x = element_blank(),
        axis.title.y = element_text(size = 17, color = "black"),
        axis.text.x = element_text(size = 15, color = "black", angle = 90, hjust = 1, vjust = 0.5),
        axis.text.y = element_text(size = 15, color = "black"))

gg.wilcox_fun


#
#### HEATMAP                ####

z.genus_prok <- t(genus_prok)[as.character(unique(mod.t6_prok$ID_well)),subset(wilcox_prok, z_value > 0)$Genus]
z.genus_fun <- t(genus_fun)[as.character(unique(mod.t6_fun$ID_well)),subset(wilcox_fun, z_value > 0)$Genus]

z.genus_tot <- merge(z.genus_prok, z.genus_fun, by = "row.names")

z.genus_tot <- merge(sample_df[,c(1,8:15)], z.genus_tot, by.x = "ID_well", by.y = "Row.names")

cor.z_df <- rcorr(as.matrix(z.genus_tot[,2:53]), type = "spearman")

cor.z.p_df <- as.matrix(cor.z_df$P)
cor.z.p_df[is.na(cor.z.p_df)] <- 1
cor.z.p_df <- apply(cor.z.p_df, 1, function(x) p.adjust(x, method = "fdr"))

cor.z.r_df <- as.matrix(cor.z_df$r)
diag(cor.z.r_df) <- NA

cor.z.p_df <- cor.z.p_df[9:52,1:8]
cor.z.r_df <- cor.z.r_df[9:52,1:8]

cor.z.r_df.plot <- merge(melt(cor.z.r_df), melt(cor.z.p_df), by = c("Var1", "Var2"))
cor.z.r_df.plot$signif <- ifelse(cor.z.r_df.plot$value.y < 0.001, "***", 
                               ifelse(cor.z.r_df.plot$value.y < 0.01, "**", ifelse(cor.z.r_df.plot$value.y < 0.05, "*", "")))

gg.heat_z <- ggplot(cor.z.r_df.plot, aes(x = Var1, y = Var2, fill = value.x)) + 
  geom_tile(color = "gray50") + 
  scale_fill_gradient2(name = "Spearman\nCorrelation", high = "#c0272f", mid = "#ffffff", low = "#3f4eaa", midpoint = 0) + 
  geom_text(aes(x = Var1, y = Var2, label = signif), color = "black", size = 4) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20)) +
  scale_y_discrete(limits=rev)+
  theme_void() + 
  theme(legend.position = "none",
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 15, color = "black", hjust = 1))

gg.heat_z

#
################################### EXPORT FIGURE                ####

gg.heat_z.legend <- ggplot(cor.z.r_df.plot, aes(x = Var1, y = Var2, fill = value.x)) + 
  geom_tile(color = "gray50") + 
  scale_fill_gradient2(name = "Spearman\nCorrelation", high = "#c0272f", mid = "#ffffff", low = "#3f4eaa", midpoint = 0) + 
  geom_text(aes(x = Var1, y = Var2, label = signif), color = "black", size = 4) +
  guides(fill = guide_colourbar(barwidth = 0.5, barheight = 20)) +
  scale_y_discrete(limits=rev)+
  theme_void() + 
  theme(legend.position = "right",
        legend.title = element_text(size = 15, color = "black"),
        legend.text = element_text(size = 15, color = "black"),
        axis.title = element_blank(),
        axis.text.x = element_text(size = 15, angle = 90, vjust = 0.5, hjust = 1, color = "black"),
        axis.text.y = element_text(size = 15, color = "black", hjust = 1))

gg.heat_z.legend <- get_legend(gg.heat_z.legend)

gg.figureS9 <- plot_grid(plot_grid(gg.wilcox_prok, gg.wilcox_fun, nrow = 1, labels = c("A", "B"), label_size = 18, rel_widths = c(1, 0.65)),
                         plot_grid(gg.heat_z, gg.heat_z.legend, nrow = 1, rel_widths = c(1, 0.15)), ncol = 1, labels = c("", "C"), label_size = 18)
gg.figureS9

ggsave("../Manuscript/Figures/Supp_Figure_S9.png", gg.figureS9, bg = "white", width = 12, height = 12)

#
################################################################################  SAVE WORKSPACE              ####

load("Outputs/Roof_data.RData")
#save.image("Outputs/Roof_data.RData")
