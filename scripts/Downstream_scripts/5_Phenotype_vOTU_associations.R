################################################################################
##### LL-NEXT: Associations between viral features and phenotypes
### Author(s):Asier Fernández-Pato
### Last updated: 14th December, 2025
################################################################################

#****************
# Load modules
#****************
library(lme4)
library(lmerTest)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(vegan)


#****************
# Define functions
#****************

# Function to test associations between vOTU features and phenotypes using mixed-effects models
# It includes the option to correct for the predicted bacterial host abundance
associate_features_with_phenotypes <- function(sample_metadata, subject_id_col,
                                               feature_data, phenotype_list, covariates = NULL, bacterial_host_data = NULL,
                                               adjust_for_bacterial_host = FALSE, random_effect = "NEXT_ID") {
  
  # Prepare metadata
  df <- sample_metadata
  row.names(df) <- df[[subject_id_col]]
  
  # If feature_data is a character vector, extract from metadata
  if (is.character(feature_data)) {
    feature_data <- df[, feature_data, drop = FALSE]
  }
  
  # Remove overlapping columns to avoid duplicates
  cols_to_remove <- intersect(names(df), colnames(feature_data))
  df_clean <- df[, !names(df) %in% cols_to_remove, drop = FALSE]
  
  # Merge metadata with viral features
  df <- cbind(df_clean, feature_data)
  
  # Base covariates
  base_covariates <- c("read_depth", "DNA_concentration_ng_ul", "Timepoint_categorical")
  all_covariates <- c(base_covariates, covariates)
  fixed_effects_null <- paste(all_covariates, collapse = " + ")
  
  # Loop over each feature
  features <- colnames(feature_data)
  results_all <- tibble()
  
  for (i in seq_along(features)) {
    feat <- features[i]
    cat("Feature", i, "/", length(features), ":", feat, "\n")
    feat_safe <- paste0("`", feat, "`")
    
    # Optional: Add bacterial host abundance as covariate
    bacterial_covariate <- NULL
    if (adjust_for_bacterial_host && !is.null(bacterial_host_data)) {
      if (feat %in% colnames(bacterial_host_data)) {
        df$host_abundance_tmp <- bacterial_host_data[[feat]] # add abundance of the corresponding host
        bacterial_covariate <- "host_abundance_tmp"
      } else {
        # Host abundance missing, skip this feature completely
        cat("Skipping", feat, "because host abundance is missing.\n")
        next
      }
    } else {
      df$host_abundance_tmp <- NA
    }
    
    # Test association with each phenotype
    for (phenotype in phenotype_list) {
      phenotype_safe <- paste0("`", phenotype, "`")
      
      # Keep only rows with non-missing phenotype and covariates
      non_missing <- !is.na(df[[phenotype]])
      if (!is.null(covariates)) {
        for (cov in covariates) {
          non_missing <- non_missing & !is.na(df[[cov]])
        }
      }
      if (adjust_for_bacterial_host && !is.null(bacterial_covariate)) {
        non_missing <- non_missing & !is.na(df$host_abundance_tmp)
      }
      df_sub <- df[non_missing, ]
      if (nrow(df_sub) < 3) next # skip if too few samples
      
      # Build model formulas
      fixed_effects <- fixed_effects_null
      if (!is.null(bacterial_covariate)) {
        fixed_effects <- paste(fixed_effects, "+", bacterial_covariate)
      }
      
      null_formula <- as.formula(paste(feat_safe, "~", fixed_effects, "+ (1|", random_effect, ")"))
      full_formula <- as.formula(paste(feat_safe, "~", fixed_effects, "+", phenotype_safe, "+ (1|", random_effect, ")"))
      
      # Fit models
      model_null <- try(lmer(null_formula, df_sub, REML = FALSE), silent = TRUE)
      if (inherits(model_null, "try-error")) next
      model_full <- try(lmer(full_formula, df_sub, REML = FALSE), silent = TRUE)
      if (inherits(model_full, "try-error")) next
      
      # Likelihood ratio test
      lrt <- anova(model_full, model_null)
      p_val <- lrt$`Pr(>Chisq)`[2]
      
      # Extract coefficients
      coef_info <- summary(model_full)$coefficients
      coef_table <- as.data.frame(coef_info)[grep(phenotype, rownames(coef_info)), ]
      
      # Store results
      temp_result <- coef_table %>%
        rownames_to_column("Term") %>%
        as_tibble() %>%
        mutate(
          P = p_val,
          Model = "Mixed",
          Feature = feat,
          Phenotype = phenotype
        )
      
      results_all <- bind_rows(results_all, temp_result)
    }
  }
  
  # Multiple testing correction
  results_all <- results_all %>%
    mutate(FDR = p.adjust(P, method = "BH"))
  
  return(results_all)
}

    
#****************
# Read input data
#****************

# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_analysis/VIRUSES/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_07052025.txt")
Abundance_table_infants <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_INFANTS_07052025.txt")
Abundance_table_mothers <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_MOTHERS_07052025.txt")
Abundance_table_infants_host <- read.delim("Abundance_table/Bacterial_genus_host_vOTU_abundance_INFANTS_17092025.txt")
Abundance_table_mothers_host <- read.delim("Abundance_table/Bacterial_genus_host_vOTU_abundance_MOTHERS_17092025.txt")
Abundance_table_bacteria <- read.delim("Abundance_table/Bacterial_abundance_GTDB_genus_28102025.txt", check.names = F)

Sample_metadata <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt")
Sample_metadata2 <- read.delim("Metadata_NEXT/LLNEXT_metadata_15_04_2024.txt")
Sample_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_07052025.txt")
Sample_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_07052025.txt")

Virus_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")
Virus_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_INFANTS_07052025.txt")
Virus_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_MOTHERS_07052025.txt")

# Make sure Phyllobacterium is not present in Abundance tables by host
# Abundance_table_infants_host <- Abundance_table_infants_host[, !colnames(Abundance_table_infants_host) %in% "Phyllobacterium"]
# Abundance_table_mothers_host <- Abundance_table_mothers_host[, !colnames(Abundance_table_mothers_host) %in% "Phyllobacterium"]

##************************************************************************
# 1. Section of phenotypes
##************************************************************************

pheno_selection <- read.delim("9_ASSOCIATION_PHENOTYPES/masterfile_cross_sectional_selection_AFTER_correlation_&_correction_09_05_2024.txt")
maternal_phenos <- na.omit(pheno_selection[pheno_selection$TCAM_mothers == 1, "variable_name"])
infant_phenos <- na.omit(pheno_selection[pheno_selection$TCAM_infants == 1, "variable_name"])
all_phenos <- union(maternal_phenos, infant_phenos)
all_phenos <- setdiff(all_phenos, c("infant_relations", "twin_pair")) # Exclude less interesting phenotypes

# Add phenotypes of interest to Sample metadata (infants)
phenotypes <- read.delim("9_ASSOCIATION_PHENOTYPES/masterfile_cross_sectional_2023_11_15.txt")
phenotypes_interest <- phenotypes[, all_phenos]

# Remove any overlapping columns before joining
cols_to_remove <- intersect(names(Sample_metadata_infants), names(phenotypes_interest))
Sample_metadata_infants <- Sample_metadata_infants[, !names(Sample_metadata_infants) %in% cols_to_remove, drop = FALSE]

# Join by matching NEXT_ID to next_id_infant
Sample_metadata_infants <- left_join(Sample_metadata_infants, phenotypes_interest,
  by = c("NEXT_ID" = "next_id_infant"))

# Add also bacterial shannon diversity
Sample_metadata2 <- Sample_metadata2[,c("NG_ID", "shannon")]
colnames(Sample_metadata2) <- c("NG_ID", "Bacterial_shannon")
Sample_metadata_infants <- left_join(Sample_metadata_infants, Sample_metadata2,
                                     by = "NG_ID")

# Generate a new variable birth_delivery_mode_CS
Sample_metadata_infants$birth_delivery_mode_CS <- ifelse(
  Sample_metadata_infants$birth_delivery_mode_simple == "post_labor_CS", "post_labor_CS",
  ifelse(Sample_metadata_infants$birth_delivery_mode_simple == "pre_labor_CS", "pre_labor_CS", NA)
)

# Set parity as a binary phenotype (no if 0, yes if > 0 )
Sample_metadata_infants <- Sample_metadata_infants %>%
  mutate(
    parity_num = as.numeric(as.character(mother_birthcard_parity)),
    mother_birthcard_parity = if_else(parity_num > 0, "yes", "no")
  )


##************************************************************************
# 2. Transformation of abundances and/or features/phenotypes
##************************************************************************

# Identify categorical and continuous variables
categorical_vars <- names(Sample_metadata_infants)[
  sapply(Sample_metadata_infants, function(x)
    is.character(x) || (is.numeric(x) && length(unique(na.omit(x))) < 10))
]

continuous_vars <- setdiff(names(Sample_metadata_infants), categorical_vars)

# Set as factors categorical variables
Sample_metadata_infants <- Sample_metadata_infants %>%
  mutate(across(all_of(categorical_vars), as.factor))

# Recode birth_delivery_mode_simple and birth_deliverybirthcard_mode_binary (VG as reference)
Sample_metadata_infants$birth_delivery_mode_simple <- factor(
  Sample_metadata_infants$birth_delivery_mode_simple,
  levels = c("VG", "post_labor_CS", "pre_labor_CS")
)

Sample_metadata_infants$birth_deliverybirthcard_mode_binary <- factor(
  Sample_metadata_infants$birth_deliverybirthcard_mode_binary,
  levels = c("VG", "CS")
)

################################
# Log-transform vOTU abundances
################################

# Get the list of vOTUs present in more than 5% samples
Abundance_table_infants <- t(Abundance_table_infants)
vOTU_nonZero  <- colSums(Abundance_table_infants>0) / nrow(Abundance_table_infants) 
vOTU_keep <- colnames(Abundance_table_infants)[vOTU_nonZero > 0.05]

# Log-transform viral abundance tables (including pseudocount)
pseudocount_vOTU_infants <- min(Abundance_table_infants[Abundance_table_infants > 0]) / 2
Abundance_table_infants_log <- Abundance_table_infants
Abundance_table_infants_log[Abundance_table_infants_log == 0] <- pseudocount_vOTU_infants
Abundance_table_infants_log <- log(Abundance_table_infants_log)

# Remove the vOTUs that are not prevalent
Abundance_table_infants_log_filtered <- as.data.frame(Abundance_table_infants_log[, vOTU_keep])

# Order abundance table according to Sample metadata
Abundance_table_infants_log_filtered <- Abundance_table_infants_log_filtered[match(Sample_metadata_infants$NG_ID,
                                                                                   rownames(Abundance_table_infants_log_filtered)), ]
#######################################
# Log-transform vOTU abundances by host
#######################################

# Get the list of host groups present in more than 5% samples
vOTU_nonZero_host <- colSums(Abundance_table_infants_host > 0) / nrow(Abundance_table_infants_host)
vOTU_keep_host <- colnames(Abundance_table_infants_host)[vOTU_nonZero_host > 0.05]

# Log-transform viral abundance tables (including pseudocount)
pseudocount_vOTU_infants_host <- min(Abundance_table_infants_host[Abundance_table_infants_host > 0]) / 2
Abundance_table_infants_host_log <- Abundance_table_infants_host
Abundance_table_infants_host_log[Abundance_table_infants_host_log == 0] <- pseudocount_vOTU_infants_host
Abundance_table_infants_host_log <- log(Abundance_table_infants_host_log)

# Subset to prevalent host-linked vOTUs
Abundance_table_infants_host_log_filtered <- as.data.frame(Abundance_table_infants_host_log[, vOTU_keep_host])

# Order abundance table according to Sample metadata
Abundance_table_infants_host_log_filtered  <- Abundance_table_infants_host_log_filtered [match(Sample_metadata_infants$NG_ID,
                                                                                   rownames(Abundance_table_infants_host_log_filtered)), ]

#######################################
# CLR-transform bacterial abundances
#######################################

# Subset infant from bacterial abundance table
Abundance_table_bacteria_infants <- Abundance_table_bacteria[rownames(Abundance_table_bacteria) %in% Sample_metadata_infants$NG_ID,]

# Get the list of genera to keep based on their vOTU prevalence
genus_keep <- vOTU_keep_host

# CLR-transform bacterial relative abundances
pseudocount_bact_infants <- min(Abundance_table_bacteria_infants[Abundance_table_bacteria_infants!=0])/2
Abundance_table_bacteria_infants_CLR <- decostand(Abundance_table_bacteria_infants, "clr", pseudocount=pseudocount_bact_infants)

# Subset to hosts of interests
genera_interest <- genus_keep[genus_keep %in% colnames(Abundance_table_bacteria_infants_CLR)]
Abundance_table_bacteria_infants_CLR_filtered <- Abundance_table_bacteria_infants_CLR[, genera_interest]

# Order abundance table according to Sample metadata
Abundance_table_bacteria_infants_CLR_filtered  <- Abundance_table_bacteria_infants_CLR_filtered [match(Sample_metadata_infants$NG_ID,
                                                                                               rownames(Abundance_table_bacteria_infants_CLR_filtered)), ]


##************************************************************************
# 3. Association with viral diversity and other viral features
##************************************************************************

# Select viral features of interest
Viral_features <- c("richness", "shannon", "relab_temperate_phages")
Phenotypes <- setdiff(colnames(phenotypes_interest), c("next_id_infant", "next_id_mother",
                                                       "FAMILY", "birth_delivery_mode_complex",
                                                        "infant_birthcard_feeding_mode_after_birth",
                                                       "mother_birthcard_duration_water_broken_min"))
Phenotypes <- c(Phenotypes, "birth_delivery_mode_CS")

##########################
# A. Unadjusted
##########################

# Run the association 
results_viral_features_unadj <- associate_features_with_phenotypes(
  sample_metadata = Sample_metadata_infants,
  subject_id_col = "NG_ID",
  feature_data = Viral_features,
  phenotype_list = Phenotypes
)

# Order by FDR
results_viral_features_unadj <- results_viral_features_unadj %>%
  arrange(FDR)

# Save results
write.table(results_viral_features_unadj, "9_ASSOCIATION_PHENOTYPES/RESULTS/Associations_Viral_features_Phenotypes_INFANTS_LOG_NO_adj.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


##########################
# B. Adjusted (Feeding + delivery)
##########################

# Run the association 
results_viral_features_adj <- associate_features_with_phenotypes(
  sample_metadata = Sample_metadata_infants,
  subject_id_col = "NG_ID",
  feature_data = Viral_features,
  phenotype_list = setdiff(Phenotypes, c("birth_deliverybirthcard_mode_binary", "birth_delivery_mode_simple",
                                         "infant_ffq_ever_never_breastfed")),
  covariates = c("birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed")
)

# Order by FDR
results_viral_features_adj <- results_viral_features_adj %>%
  arrange(FDR)

# Save results
write.table(results_viral_features_adj, "9_ASSOCIATION_PHENOTYPES/RESULTS/Associations_Viral_features_Phenotypes_INFANTS_LOG_adj.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


##########################
# C. Adjusted (Feeding + delivery + Bacterial shannon)
##########################

# Run the association 
results_viral_features_adj_bac <- associate_features_with_phenotypes(
  sample_metadata = Sample_metadata_infants,
  subject_id_col = "NG_ID",
  feature_data = Viral_features,
  phenotype_list = setdiff(Phenotypes, c("birth_deliverybirthcard_mode_binary", "birth_delivery_mode_simple",
                                         "infant_ffq_ever_never_breastfed")),
  covariates = c("birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed", "Bacterial_shannon")
)

# Order by FDR
results_viral_features_adj_bac <- results_viral_features_adj_bac %>%
  arrange(FDR)

# Save results
write.table(results_viral_features_adj_bac, "9_ASSOCIATION_PHENOTYPES/RESULTS/Associations_Viral_features_Phenotypes_INFANTS_LOG_adj_BAC.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


##************************************************************************
# 4. Association with vOTU abundances / presence
##************************************************************************

##########################
# A. Abundances: Not adjusted (Feeding + delivery)
##########################

# Run the association 
results_vOTUs_unadj <- associate_features_with_phenotypes(
  sample_metadata = Sample_metadata_infants,
  subject_id_col = "NG_ID",
  feature_data = Abundance_table_infants_log_filtered,
  phenotype_list = Phenotypes)

# Order results by FDR
results_vOTUs_unadj <- results_vOTUs_unadj %>%
  arrange(FDR)

# Save results
write.table(results_vOTUs_unadj, "9_ASSOCIATION_PHENOTYPES/RESULTS/Associations_vOTU_abundances_Phenotypes_INFANTS_LOG_unadj.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)


##########################
# B. Abundances: adjusted (Feeding + delivery)
##########################

# Run the association (using transformed abundances)
results_vOTUs_adj <- associate_features_with_phenotypes(
  sample_metadata = Sample_metadata_infants,
  subject_id_col = "NG_ID",
  feature_data = Abundance_table_infants_log_filtered,
  phenotype_list = setdiff(Phenotypes, c("birth_deliverybirthcard_mode_binary", "birth_delivery_mode_simple",
                                         "infant_ffq_ever_never_breastfed")),
  covariates = c("birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed")
)

# Order by FDR
results_vOTUs_adj <- results_vOTUs_adj %>%
  arrange(FDR)

# Save results
write.table(results_vOTUs_adj,"9_ASSOCIATION_PHENOTYPES/RESULTS/Associations_vOTU_abundances_Phenotypes_INFANTS_LOG_adj.txt",
  sep = "\t", quote = FALSE, row.names = FALSE)


##************************************************************************
# 5. Association with host-linked vOTU abundances/presence
##************************************************************************

##########################
# A. Abundances: Not adjusted (Feeding + delivery)
##########################

# Run the association 
results_vOTUs_host_unadj <- associate_features_with_phenotypes(
  sample_metadata = Sample_metadata_infants,
  subject_id_col = "NG_ID",
  feature_data = Abundance_table_infants_host_log_filtered,
  phenotype_list = Phenotypes)

# Order results by FDR
results_vOTUs_host_unadj <- results_vOTUs_host_unadj %>%
  arrange(FDR)

# Save results
write.table(results_vOTUs_host_unadj, "9_ASSOCIATION_PHENOTYPES/RESULTS/Associations_Host_vOTU_abundances_Phenotypes_INFANTS_LOG_unadj.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

##########################
# B. Abundances: Not adjusted (Feeding + delivery) but adjusted (Host)
##########################

# Run the association 
results_vOTUs_host_unadj_host <- associate_features_with_phenotypes(
  sample_metadata = Sample_metadata_infants,
  subject_id_col = "NG_ID",
  feature_data = Abundance_table_infants_host_log_filtered,
  phenotype_list = Phenotypes,
  bacterial_host_data = Abundance_table_bacteria_infants_CLR_filtered,
  adjust_for_bacterial_host = TRUE  
)

# Order results by FDR
results_vOTUs_host_unadj_host <- results_vOTUs_host_unadj_host %>%
  arrange(FDR)

# Save results
write.table(results_vOTUs_host_unadj_host, "9_ASSOCIATION_PHENOTYPES/RESULTS/Associations_Host_vOTU_abundances_Phenotypes_INFANTS_LOG_unadj_Host_adj.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

##########################
# C. Abundances: adjusted (Feeding + delivery)
##########################

# Run the association 
results_vOTUs_host_adj <- associate_features_with_phenotypes(
  sample_metadata = Sample_metadata_infants,
  subject_id_col = "NG_ID",
  feature_data = Abundance_table_infants_host_log_filtered,
  phenotype_list = setdiff(Phenotypes, c("birth_deliverybirthcard_mode_binary", "birth_delivery_mode_simple",
                                         "infant_ffq_ever_never_breastfed")),
  covariates = c("birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed")
)

# Order results by FDR
results_vOTUs_host_adj <- results_vOTUs_host_adj %>%
  arrange(FDR)

# Save results
write.table(results_vOTUs_host_adj, "9_ASSOCIATION_PHENOTYPES/RESULTS/Associations_Host_vOTU_abundances_Phenotypes_INFANTS_LOG_adj.txt",
  sep = "\t", quote = FALSE, row.names = FALSE)

##########################
# D. Abundances: adjusted (Feeding + delivery + Bacterial genus abundance)
##########################

results_vOTUs_host_adj_bact <- associate_features_with_phenotypes(
  sample_metadata = Sample_metadata_infants,
  subject_id_col = "NG_ID",
  feature_data = Abundance_table_infants_host_log_filtered,
  phenotype_list = setdiff(Phenotypes, c("birth_deliverybirthcard_mode_binary", "birth_delivery_mode_simple",
                                         "infant_ffq_ever_never_breastfed")),
  covariates = c("birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed"),
  bacterial_host_data = Abundance_table_bacteria_infants_CLR_filtered,
  adjust_for_bacterial_host = TRUE  
)

# Order results by FDR
results_vOTUs_host_adj_bact <- results_vOTUs_host_adj_bact %>%
  arrange(FDR)


##************************************************************************
# 6. Plots
##************************************************************************

##########################
# A. Associations with viral diversity
##########################

Sample_metadata_infants$birth_delivery_mode_simple <- factor(
  Sample_metadata_infants$birth_delivery_mode_simple,
  levels = c("VG", "pre_labor_CS", "post_labor_CS")
)

# Generate plots for infant delivery (3 levels) and feeding mode
pdf('9_ASSOCIATION_PHENOTYPES/Plots/Richness_Delivery.pdf', width=3, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$birth_delivery_mode_simple),], 
       aes(x = birth_delivery_mode_simple, y = richness, fill = birth_delivery_mode_simple)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Delivery mode', y = 'Viral Richness') +
  scale_fill_manual(values = c("VG" = "#B6D7A8", "post_labor_CS" = "#7FA6D9",
                               "pre_labor_CS" = "#A2C2E5")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none")
dev.off()

pdf('9_ASSOCIATION_PHENOTYPES/Plots/Shannon_Delivery.pdf', width=2.9, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$birth_delivery_mode_simple),], 
       aes(x = birth_delivery_mode_simple, y = shannon, fill = birth_delivery_mode_simple)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Delivery mode', y = 'Shannon Index') +
  scale_fill_manual(values = c("VG" = "#B6D7A8", "post_labor_CS" = "#7FA6D9",
                               "pre_labor_CS" = "#A2C2E5")) +
  scale_y_continuous(breaks = seq(0, 5, 1),limits = c(0, 5),
                     expand = expansion(mult = c(0.1, 0.2))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()

pdf('9_ASSOCIATION_PHENOTYPES/Plots/Richness_Feeding.pdf', width=2.9, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$infant_ffq_ever_never_breastfed),], 
       aes(x = infant_ffq_ever_never_breastfed, y = richness, fill = infant_ffq_ever_never_breastfed)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Feeding mode', y = 'Viral Richness') +
  scale_fill_manual(values = c("ever_BF" = "#A37BAA", "never_BF" = "#D8C9A9")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()

pdf('9_ASSOCIATION_PHENOTYPES/Plots/Shannon_Feeding.pdf', width=2.7, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$infant_ffq_ever_never_breastfed),], 
       aes(x = infant_ffq_ever_never_breastfed, y = shannon, fill = infant_ffq_ever_never_breastfed)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Feeding mode', y = 'Shannon Index') +
  scale_fill_manual(values = c("ever_BF" = "#A37BAA", "never_BF" = "#D8C9A9")) +
  scale_y_continuous(breaks = seq(0, 5, 1),limits = c(0, 5),
                     expand = expansion(mult = c(0.1, 0.2))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()

# Generate plots for maternal infections, net income and food allergy
pdf('9_ASSOCIATION_PHENOTYPES/Plots/Shannon_Fungal_inf.pdf', width=2.7, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$mother_deliveryhealth_preg_complaint_vaginal_fungal_inf),], 
       aes(x = mother_deliveryhealth_preg_complaint_vaginal_fungal_inf, y = shannon, fill = mother_deliveryhealth_preg_complaint_vaginal_fungal_inf)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Fungal infection', y = 'Shannon Index') +
  scale_fill_manual(values = c("yes" = "#A8B89F", "no" = "#C4BBAF")) +
  scale_y_continuous(breaks = seq(0, 5, 1),limits = c(0, 5),
                     expand = expansion(mult = c(0.1, 0.2))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()

pdf('9_ASSOCIATION_PHENOTYPES/Plots/Shannon_net_income.pdf', width=2.9, height=3.2)
Sample_metadata_infants$mother_income_net_p18 <- factor(Sample_metadata_infants$mother_income_net_p18)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$mother_income_net_p18),], 
       aes(x = mother_income_net_p18, y = shannon, fill = mother_income_net_p18)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Net income', y = 'Shannon Index') +
  scale_fill_manual(values = c("0" = "#9FB8C1", "1" = "#D9B38C", "2" = "#E6A6B0")) +
  scale_y_continuous(breaks = seq(0, 5, 1),limits = c(0, 5),
                     expand = expansion(mult = c(0.1, 0.2))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()

pdf('9_ASSOCIATION_PHENOTYPES/Plots/Richness_Allergy.pdf', width=2.8, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$infant_health_food_allergy),], 
       aes(x = infant_health_food_allergy, y = richness, fill = infant_health_food_allergy)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Food Allergy', y = 'Viral Richness') +
  scale_fill_manual(values = c("no" = "#4C72B0", "yes" = "#DD8452")) +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.2))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()

pdf('9_ASSOCIATION_PHENOTYPES/Plots/Shannon_Allergy.pdf', width=2.7, height=3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$infant_health_food_allergy),], 
       aes(x = infant_health_food_allergy, y = shannon, fill = infant_health_food_allergy)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Food Allergy', y = 'Shannon Index') +
  scale_fill_manual(values = c("no" = "#4C72B0", "yes" = "#DD8452")) +
  scale_y_continuous(breaks = seq(0, 5, 1),limits = c(0, 5),
                     expand = expansion(mult = c(0.1, 0.2))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()

##########################
# B. Associations with viral abundances
##########################

# Boxplot of Bacteroides phage abundances vs delivery mode
plot_Bacteroides <- Sample_metadata_infants %>%
  select(birth_delivery_mode_simple, Timepoint_categorical) %>%
  bind_cols(Bacteroides_phage_abundance = Abundance_table_infants_host_log_filtered$Bacteroides) %>%
  filter(!is.na(birth_delivery_mode_simple), !is.na(Timepoint_categorical), !is.na(Bacteroides_phage_abundance))

plot_Bacteroides$Timepoint_categorical <- factor(
  plot_Bacteroides$Timepoint_categorical,
  levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12"))

plot_Bacteroides$birth_delivery_mode_simple <- factor(
  plot_Bacteroides$birth_delivery_mode_simple,
  levels = c("VG", "pre_labor_CS", "post_labor_CS"))

pdf('9_ASSOCIATION_PHENOTYPES/Plots/Bacteroides_Delivery.pdf', width=7.2, height=2.8)
ggplot(plot_Bacteroides, aes(x = Timepoint_categorical, y = Bacteroides_phage_abundance,
                             fill = birth_delivery_mode_simple)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8),
              alpha = 0.4, color = "grey70", size = 1.8) +
  geom_boxplot(position = position_dodge(width = 0.8), width = 0.6,
               color = "black", outlier.color = NA, size = 0.8) +
  scale_fill_manual(values = c("VG" = "#B6D7A8", "post_labor_CS" = "#7FA6D9", "pre_labor_CS" = "#A2C2E5")) +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
  labs(x = "Infant Timepoint",y = "Bacteroides phage abundance (log RPKM)", fill = "Delivery mode") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )
dev.off()

# Generate a heatmap with associations with viral abundances (by host)
# Feeding and delivery mode: results from results_vOTUs_host_unadj
# Other phenotypes: results from results_vOTUs_host_adj
# 1. Combine adjusted and unadjusted results
association_heatmap_infant_host <- results_vOTUs_host_unadj[
  results_vOTUs_host_unadj$Phenotype %in% c("birth_deliverybirthcard_mode_binary",
                                            "infant_ffq_ever_never_breastfed"), ]

association_heatmap_infant_host <- rbind(association_heatmap_infant_host,
  results_vOTUs_host_adj[!results_vOTUs_host_adj$Phenotype %in%
      c("birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed"), ]
)

# 2. Set phenotype factor levels
phenotype_levels <- unique(association_heatmap_infant_host$Phenotype)
association_heatmap_infant_host$Phenotype <- factor(
  association_heatmap_infant_host$Phenotype,
  levels = phenotype_levels
)

# 3. Filter out redundant terms (multi-level phenotypes except for blood group)
remove_terms <- c("mother_birthcard_parity1","mother_birthcard_parity3",
  "mother_birthcard_parity4","mother_birthcard_parity7",
  "mother_income_net_p181","mother_education_p181",
  "mother_ffq_aMED_score1","mother_ffq_aMED_score2","mother_ffq_aMED_score3",
  "mother_ffq_aMED_score4","mother_ffq_aMED_score5","mother_ffq_aMED_score7",
  "mother_ffq_aMED_score8")

association_heatmap_infant_host <- association_heatmap_infant_host %>%
  filter(!Term %in% remove_terms)

# 4. Keep only phenotypes with at least one nominally significant association (P < 0.05)
signif_host_phenotypes <- association_heatmap_infant_host %>%
  group_by(Phenotype) %>%
  summarise(any_signif = any(P < 0.05)) %>%
  filter(any_signif) %>%
  pull(Phenotype)

association_heatmap_infant_host <- association_heatmap_infant_host %>%
  filter(Phenotype %in% signif_host_phenotypes)

# 5. Prepare data for heatmap
host_heatmap_long <- association_heatmap_infant_host %>%
  mutate(Host = Feature, sig = case_when(
      FDR < 0.05 ~ "*",
      P < 0.05 & FDR >= 0.05 ~ "~",
      TRUE ~ "")
  )

# 6. Add white gap between adjusted vs unadjusted phenotypes
special_phenos <- c("birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed")

# Create a spacer row (blank row)
spacer <- data.frame(Host = unique(host_heatmap_long$Host)[1],
                     Term = " ", `t value` = NA, sig = "", Phenotype = "spacer")

# Combine and order with the spacer in between
host_heatmap_long_gap <- host_heatmap_long %>%
  mutate(group = ifelse(Phenotype %in% special_phenos, "B", "A")) %>%
  bind_rows(spacer) %>%
  arrange(group, Phenotype)

# Set y-axis order so the gap is preserved
host_heatmap_long_gap$Term <- factor(host_heatmap_long_gap$Term,
  levels = c(unique(host_heatmap_long_gap$Term[host_heatmap_long_gap$group == "A" &
                                        host_heatmap_long_gap$Term != " "]), " ",
    unique(host_heatmap_long_gap$Term[host_heatmap_long_gap$group == "B"]))
)

pdf("9_ASSOCIATION_PHENOTYPES/Plots/vOTU_Host_Phenotype_association.pdf", width = 14, height = 12)
ggplot(host_heatmap_long_gap, aes(x = Host, y = Term)) +
  geom_tile(aes(fill = `t value`), colour = "lightgrey", linewidth = 0.5) +
  geom_text(aes(label = sig), color = "#333333", size = 5, vjust = 0.5) +
  scale_fill_gradient2(
    low = "#003F45", mid = "#F0F0F0", high = "#F48C42",
    midpoint = 0, name = "T-value", na.value = "white"
  ) +
  scale_y_discrete(drop = FALSE) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = 11),
    axis.title.x = element_blank(),
    panel.grid = element_blank(),
    axis.text.y = element_text(size = 12)
  ) +
  xlab("Bacterial host genus") +
  ylab("Phenotype (specific term)")
dev.off()


##************************************************************************
# 7. Write results
##************************************************************************
# Sample metadata
write.table(Sample_metadata_infants,"Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_phenos_07052025.txt",
            sep = "\t", row.names = F, quote = FALSE) 

# Phenotypes for association
write.table(Phenotypes,"9_ASSOCIATION_PHENOTYPES/Phenos_of_interest.txt",
            sep = "\t", row.names = F, quote = FALSE) 

