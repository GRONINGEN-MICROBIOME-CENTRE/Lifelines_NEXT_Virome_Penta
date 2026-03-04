################################################################################
##### LL-NEXT: vOTU persistance: Phage mechanisms
### Author(s): Asier Fernández-Pato
### Last updated: 12th December, 2025
################################################################################

#****************
# Load libraries
#****************
library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(broom.mixed)
library(data.table)
library(ggplot2)
library(glmmTMB)
library(stringr)
library(scales)

#****************
# Define functions
#****************

# Function to perform GLM associations with binary persistence (adjusting for abundance)
run_glm_binary <- function(metadata, abundance_col, features) {
  
  # Initialize results table
  glm_results <- data.frame(
    system = features,
    estimate = NA,
    std_error = NA,
    z_value = NA,
    p_value = NA,
    stringsAsFactors = FALSE
  )
  
  # Loop over each feature
  for (i in seq_along(features)) {
    sys <- features[i]
    formula <- as.formula(paste0("Persistent ~ ", sys, " + ", abundance_col))
    
    # Fit GLM (logistic regression)
    model <- glm(formula, data = metadata, family = binomial)
    coef_summary <- summary(model)$coefficients
    
    # Handle possible row name differences
    rowname <- paste0(sys, "Yes")
    if (!rowname %in% rownames(coef_summary)) rowname <- sys
    
    # Store results
    glm_results[i, c("estimate","std_error","z_value","p_value")] <-
      coef_summary[rowname, c("Estimate","Std. Error","z value","Pr(>|z|)")]
  }
  
  # FDR correction
  glm_results$FDR <- p.adjust(glm_results$p_value, method = "BH")
  
  return(glm_results)
}

# Function to perform a GLMM for associations of viral features with binary persistence (adjusting for abundance and host as random effect)
run_glmmTMB_binary <- function(metadata, abundance_col, DGR_ADF_cols) {
  
  # Ensure predictors are factors
  metadata[DGR_ADF_cols] <- lapply(metadata[DGR_ADF_cols], factor)
  
  # Initialize results table
  glmm_results <- data.frame(
    system = character(),
    estimate = numeric(),
    std_error = numeric(),
    z_value = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop over each system
  for (sys in DGR_ADF_cols) {
    formula <- as.formula(paste0("Persistent ~ ", sys, " + ", abundance_col, " + (1 | Bacterial_genus_host)"))
    
    # Fit GLMM
    model <- glmmTMB(formula, data = metadata, family = binomial)
    
    # Extract coefficients
    coef_table <- summary(model)$coefficients$cond
    rowname <- intersect(c(paste0(sys, "Yes"), sys), rownames(coef_table))
    
    if (length(rowname) > 0) {
      glmm_results <- rbind(glmm_results, data.frame(
        system = sys,
        estimate = coef_table[rowname, "Estimate"],
        std_error = coef_table[rowname, "Std. Error"],
        z_value = coef_table[rowname, "z value"],
        p_value = coef_table[rowname, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # FDR correction
  glmm_results$FDR <- p.adjust(glmm_results$p_value, method = "BH")
  glmm_results <- glmm_results[order(glmm_results$FDR), ]
  
  return(glmm_results)
}


# Function to perform a GLMM for associations of viral features with the number of persistent individuals (adjusting for abundance and host as random effect)
run_glmmTMB_nb <- function(metadata, response_col, abundance_col, DGR_ADF_cols) {
  
  # Ensure predictors are factors
  metadata[DGR_ADF_cols] <- lapply(metadata[DGR_ADF_cols], factor)
  
  # Initialize results table
  nb_results <- data.frame(
    system = character(),
    estimate = numeric(),
    std_error = numeric(),
    z_value = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )
  
  # Loop over each system
  for (sys in DGR_ADF_cols) {
    formula <- as.formula(paste0(response_col, " ~ ", sys, " + ", abundance_col, " + (1 | Bacterial_genus_host)"))
    
    # Fit zero-inflated negative binomial GLMM
    model <- glmmTMB(formula, data = metadata, family = nbinom2, ziformula = ~1)
    
    coef_table <- summary(model)$coefficients$cond
    row_match <- grep(paste0("^", sys, "Yes$"), rownames(coef_table))
    
    if (length(row_match) == 1) {
      nb_results <- rbind(nb_results, data.frame(
        system = sys,
        estimate = coef_table[row_match, "Estimate"],
        std_error = coef_table[row_match, "Std. Error"],
        z_value = coef_table[row_match, "z value"],
        p_value = coef_table[row_match, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      ))
    }
  }
  
  # FDR correction
  nb_results$FDR <- p.adjust(nb_results$p_value, method = "BH")
  nb_results <- nb_results[order(nb_results$FDR), ]
  
  return(nb_results)
}


#****************
# Load data
#****************

# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_Analysis/VIRUSES/")

#Load metadata and abundance table for the LL-NEXT samples 

# Read metadata tables 
Sample_metadata <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt")
Sample_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_07052025.txt")
Sample_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_07052025.txt")

# Read abundance table
Abundance_table <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_07052025.txt")
Abundance_table_infants <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_INFANTS_07052025.txt")
Abundance_table_mothers <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_MOTHERS_07052025.txt")

# Read viral metadata tables
Viral_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")
Viral_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_INFANTS_07052025.txt")
Viral_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_MOTHERS_07052025.txt")


#**************************************************************
# 1. Associations of phage traits with vOTU-level persistence
#**************************************************************

# We explore the associations between vOTU persistence and:
# -- A) Viral lifestyle
# -- B) Anti-defense systems
# -- C) DGRs

###############
# A. INFANTS
###############

# Ensure Persistent is binary
Viral_metadata_infants$Persistent <- factor(Viral_metadata_infants$Persistent, levels = c("No", "Yes"))

# Encode Lifestyle as binary
Viral_metadata_infants$Lifestyle_binary <- ifelse(
  is.na(Viral_metadata_infants$Lifestyle),NA,
  ifelse(Viral_metadata_infants$Lifestyle == "Temperate", "Yes", "No")
)

# Convert to factor, keeping NAs
Viral_metadata_infants$Lifestyle_binary <- factor(Viral_metadata_infants$Lifestyle_binary,
                                                  levels = c("No", "Yes"))

#-------------
# A.1. Association with persistence (Binary) adjusting for infant mean abundance (GLM)
#-------------

# Features to test
features <- c("Lifestyle_binary", "Antidefense_system", "DGR")

# Run the association using a GLM
glm_results_infants <- run_glm_binary(
  metadata = Viral_metadata_infants,
  abundance_col = "Mean_abundance_infants",
  features = features
)

#-------------
# A.2. Association with persistence (Binary) adjusting for infant mean abundance + viral host at genus level (GLMM)
#-------------

# Here we test ADS (including common subtypes) and DGRs

# Define columns for DGR and Anti-Defense systems
DGR_ADF_cols <- c("DGR", grep("Anti", colnames(Viral_metadata_infants), value = TRUE))

# Run the association using a binomial GLMM
glmm_results_infants <- run_glmmTMB_binary(Viral_metadata_infants, "Mean_abundance_infants", DGR_ADF_cols)


#-------------
# A.3. Association with N infants persistent correcting for viral host at genus level
#-------------

# Run the association using a Negative binomial regression with random effects
zinb_results_infants <- run_glmmTMB_nb(Viral_metadata_infants,
                                       "N_Infants_Persistent", "Mean_abundance_infants", DGR_ADF_cols)


###############
# B. MOTHERS
###############

# Ensure Persistent is binary
Viral_metadata_mothers$Persistent <- factor(Viral_metadata_mothers$Persistent, levels = c("No", "Yes"))

# Encode Lifestyle as binary
Viral_metadata_mothers$Lifestyle_binary <- ifelse(
  is.na(Viral_metadata_mothers$Lifestyle),NA,
  ifelse(Viral_metadata_mothers$Lifestyle == "Temperate", "Yes", "No")
)

# Convert to factor, keeping NAs
Viral_metadata_mothers$Lifestyle_binary <- factor(Viral_metadata_mothers$Lifestyle_binary,
  levels = c("No", "Yes"))

#-------------
# B.1. Association with persistence (Binary) adjusting for infant mean abundance (GLM)
#-------------

# Run the association using a GLM
glm_results_mothers <- run_glm_binary(
  metadata = Viral_metadata_mothers,
  abundance_col = "Mean_abundance_mothers",
  features = features
)

#-------------
# B.2. Association with persistence (Binary) adjusting for maternal mean abundance + viral host at genus level (GLMM)
#-------------

# Due to low numbers in mothers, we exclude Anti_RM_ral (n=12) and Anti_RecBCD_system (n=27) from maternal analysis
DGR_ADF_cols <- DGR_ADF_cols[!DGR_ADF_cols %in% c("Anti_RM_ral", "Anti_RecBCD_system")]

# Run the association using a binomial GLMM
glmm_results_mothers <- run_glmmTMB_binary(Viral_metadata_mothers,
                                           "Mean_abundance_mothers", DGR_ADF_cols)


#-------------
# B.3. Association with N mothers persistent correcting for viral host at genus level
#-------------

# Run the association using a Negative binomial regression with random effects
zinb_results_mothers <- run_glmmTMB_nb(Viral_metadata_mothers,
                                       "N_Mothers_Persistent", "Mean_abundance_mothers", DGR_ADF_cols)



#**************************************************************
# 2. Plots
#**************************************************************

############################################################
# A. Association of lifestyle with binary persistence
############################################################

# Encode Persistence as binary with 0/1
Viral_metadata_infants_plot <- Viral_metadata_infants
Viral_metadata_mothers_plot <- Viral_metadata_mothers

Viral_metadata_infants_plot$Persistent <-
  ifelse(Viral_metadata_infants_plot$Persistent == "Yes", 1, 0)
Viral_metadata_infants_plot$Persistent <-
  as.numeric(as.character(Viral_metadata_infants_plot$Persistent))

Viral_metadata_mothers_plot$Persistent <-
  ifelse(Viral_metadata_mothers_plot$Persistent == "Yes", 1, 0)
Viral_metadata_mothers_plot$Persistent <-
  as.numeric(as.character(Viral_metadata_mothers_plot$Persistent))

# Plot showing GLM association of lifestyle with vOTU persistence
pdf('11_PERSISTENCE//Plots/Lifestyle_binary_persistence_infants.pdf', width = 3.2, height = 3.5)
ggplot(Viral_metadata_infants_plot[!is.na(Viral_metadata_infants_plot$Lifestyle),],
       aes(x = log10(Mean_abundance_infants), y = Persistent, color = Lifestyle, fill = Lifestyle)) +
  geom_boxplot(aes(group = interaction(Persistent, Lifestyle)), position = position_dodge(width = 0.25),
               width = 0.15, alpha = 0.4, outlier.shape = NA, linewidth = 0.6, color = "black") +
  geom_smooth(method = "glm",method.args = list(family = "binomial"), se = TRUE, linewidth = 1.2, alpha = 0.25) +
  labs(x = expression(log[10]("Mean abundance")),
       y = "Probability of vOTU persistence", color = "Lifestyle", fill = "Lifestyle") +
  scale_color_manual(values = c("Virulent" = "#A3C293", "Temperate" = "#D7A67A")) +
  scale_fill_manual(values = c("Virulent" = "#A3C293", "Temperate" = "#D7A67A")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
dev.off()

# Plot showing GLM association of lifestyle with vOTU persistence
pdf('11_PERSISTENCE//Plots/Lifestyle_binary_persistence_mothers.pdf', width = 3.2, height = 3.5)
ggplot(Viral_metadata_mothers_plot[!is.na(Viral_metadata_mothers_plot$Lifestyle),],
       aes(x = log10(Mean_abundance_mothers), y = Persistent, color = Lifestyle, fill = Lifestyle)) +
  geom_boxplot(aes(group = interaction(Persistent, Lifestyle)), position = position_dodge(width = 0.25),
               width = 0.15, alpha = 0.4, outlier.shape = NA, linewidth = 0.6, color = "black") +
  geom_smooth(method = "glm",method.args = list(family = "binomial"), se = TRUE, linewidth = 1.2, alpha = 0.25) +
  labs(x = expression(log[10]("Mean abundance")),
       y = "Probability of vOTU persistence", color = "Lifestyle", fill = "Lifestyle") +
  scale_color_manual(values = c("Virulent" = "#A3C293", "Temperate" = "#D7A67A")) +
  scale_fill_manual(values = c("Virulent" = "#A3C293", "Temperate" = "#D7A67A")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
dev.off()

############################################################
# B. Association of ADS with persistence
############################################################

# Prepare infant data
zinb_plot_infants <- zinb_results_infants %>%
  filter(system != "DGR") %>%
  mutate(
    significance = ifelse(FDR <= 0.05, "FDR < 0.05", "NS"),
    cohort = "Infants"
  )

# Prepare mother data, removing systems with low counts in mothers
exclude_sys_mothers <- c("Anti_RM_ral", "Anti_RecBCD_system")
zinb_plot_mothers <- zinb_results_mothers %>%
  filter(system != "DGR", !system %in% exclude_sys_mothers) %>%
  mutate(
    significance = ifelse(FDR <= 0.05, "FDR < 0.05", "NS"),
    cohort = "Mothers"
  )

# Combine infants and mothers
zinb_plot_combined <- bind_rows(zinb_plot_infants, zinb_plot_mothers)

# Define system groups
zinb_plot_combined <- zinb_plot_combined %>%
  mutate(system_group = case_when(
    grepl("Anti_RM", system) ~ "Anti-RM",
    grepl("Antidefense", system) ~ "Defense system",
    TRUE ~ "Other systems"
  ))

zinb_plot_combined$system_group <- factor(
  zinb_plot_combined$system_group,
  levels = c("Defense system", "Anti-RM", "Other systems")
)

# Compute -log10(FDR) 
zinb_plot_combined <- zinb_plot_combined %>%
  mutate(size_fdr = -log10(FDR))

# Define color palette
col_low  <- "#D9D9D9"
col_mid  <- "#C9879E"
col_high <- "#321F3F"

# Order systems: Defense first, then Anti-RM (with Anti_RM_system first), then others
system_order <- c("Antidefense_system", "Anti_RM_system",     
  zinb_plot_combined %>% 
    filter(grepl("Anti_RM", system) & system != "Anti_RM_system") %>%
    arrange(desc(estimate)) %>%
    pull(system) %>% unique(),
  zinb_plot_combined %>%
    filter(!grepl("Anti_RM|Antidefense", system)) %>%
    arrange(desc(estimate)) %>%
    pull(system) %>% unique()
)
zinb_plot_combined$system <- factor(zinb_plot_combined$system, levels = system_order)

# Plot
pdf('11_PERSISTENCE/Plots/Combined_Infant_Mother_ZINB_Results.pdf', width = 12, height = 5)
ggplot(zinb_plot_combined, aes(x = system, y = cohort)) +
  geom_point(data = subset(zinb_plot_combined, FDR <= 0.05),
    aes(size = size_fdr), shape = 21, stroke = 1.1, fill = NA, color = "black",
    position = position_nudge(y = 0)) +
  geom_point(aes(color = z_value, size = size_fdr), shape = 16, alpha = 0.9,
    position = position_nudge(y = 0)) +
  scale_color_gradientn(colors = c(col_low, col_mid, col_high),
    values = scales::rescale(c(min(zinb_plot_combined$estimate),
                               0,max(zinb_plot_combined$estimate))),
    name = "Z-value") +
  scale_size_continuous(range = c(3, 8),name = "-log10(FDR)",
    breaks = c(1, 5, 10), limits = c(min(zinb_plot_combined$size_fdr, na.rm = TRUE),
               max(zinb_plot_combined$size_fdr, na.rm = TRUE)),
    labels = c("1", "5", "10")) +
  facet_grid(. ~ system_group, scales = "free_x", space = "free_x") +
  labs(x = "Anti-defense system", y = "Type") +
  theme_bw(base_size = 14) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.8),
    panel.spacing = unit(1.2, "lines"),
    axis.text.y = element_text(size = 14),
    axis.text.x = element_text(size = 14, angle = 45, hjust = 1),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 15),
    legend.position = "bottom"
  )
dev.off()


#****************
# Write results
#****************
write.table(glm_results_mothers,"11_PERSISTENCE/Results/GLM_Persistence_MOTHERS.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(glm_results_infants,"11_PERSISTENCE/Results/GLM_Persistence_INFANTS.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(glmm_results_mothers,"11_PERSISTENCE/Results/GLMM_Persistence_MOTHERS.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(glmm_results_infants,"11_PERSISTENCE/Results/GLMM_Persistence_INFANTS.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(zinb_results_mothers,"11_PERSISTENCE/Results/GLMM_ZING_Persistence_MOTHERS.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(zinb_results_infants,"11_PERSISTENCE/Results/GLMM_ZING_Persistence_INFANTS.txt",
            sep = "\t", row.names = F, quote = FALSE) 

