################################################################################
##### LL-NEXT: Associations between viral composition and phenotypes
### Author(s):Asier Fernández-Pato
### Last updated: 14th December, 2025
################################################################################

#****************
# Load modules
#****************
library(vegan)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(tempted)
library(scales)
library(ggExtra)
library(lme4)
library(lmerTest)
library(tidyverse)
library(foreach)


#****************
# Define functions
#****************

# Function to test associations between vOTU features and phenotypes using mixed-effects models
associate_features_with_phenotypes <- function(sample_metadata,subject_id_col,
                                               feature_data, phenotype_list, covariates = NULL, random_effect = "NEXT_ID") {
  
  # Prepare metadata with individual IDs (NG_IDs) as rownames
  df <- sample_metadata
  row.names(df) <- df[[subject_id_col]]
  
  # If feature_data is a character vector (name of variables), extract from metadata
  if (is.character(feature_data)) {
    feature_data <- df[, feature_data, drop = FALSE]
  }
  
  # Remove overlapping columns to avoid duplicates
  cols_to_remove <- intersect(names(df), colnames(feature_data))
  df_clean <- df[, !names(df) %in% cols_to_remove, drop = FALSE]
  
  # Merge metadata with feature data
  df <- cbind(df_clean, feature_data)
  
  # Generate the string for fixed effects
  base_covariates <- c("read_depth", "DNA_concentration_ng_ul", "Timepoint_categorical")
  all_covariates <- c(base_covariates, covariates)
  fixed_effects_null <- paste(all_covariates, collapse = " + ")
  
  # Features to test (can be vOTU abundances or other viral features)
  features <- colnames(feature_data)
  results_all <- tibble()
  
  # Loop over each feature
  for (i in seq_along(features)) {
    feat <- features[i]
    cat("Feature", i, "/", length(features), "\n")
    feat_safe <- paste0("`", feat, "`") # just in case the variable name has spaces
    
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
      df_sub <- df[non_missing, ]
      if (nrow(df_sub) < 3) next  # skip if there are too few samples
      
      # Build models
      null_formula <- as.formula(paste(feat_safe, "~", fixed_effects_null, "+ (1|", random_effect, ")"))
      full_formula <- as.formula(paste(feat_safe, "~", fixed_effects_null, "+", phenotype_safe, "+ (1|", random_effect, ")"))
      
      # Fit models 
      model_null <- try(lmer(null_formula, df_sub, REML = FALSE), silent = TRUE)
      if (inherits(model_null, "try-error")) next
      model_full <- try(lmer(full_formula, df_sub, REML = FALSE), silent = TRUE)
      if (inherits(model_full, "try-error")) next
      
      # Likelihood ratio test
      lrt <- anova(model_full, model_null)
      p_val <- lrt$`Pr(>Chisq)`[2]
      
      # Extract coefficient info
      coef_info <- summary(model_full)$coefficients
      coef_table <- as.data.frame(coef_info)[grep(phenotype, rownames(coef_info)), ]
      
      # Format results
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

# Function to run PERMANOVA analysis per timepoint
Call_adonis = function(phenotypes, Phenotype_list, Distance, perm=10000, cores=8, covariates="base") {
  Distance = as.matrix(Distance)
  adonis_results <- data.frame(matrix(ncol = 6, nrow = length(Phenotype_list)))
  colnames(adonis_results) <- c("Phenotype", "Df", "F", "R2", "p-value", "N")
  Phenotype_list = Phenotype_list[!is.na(Phenotype_list)]
  
  # Base formula 
  base_formula <- "read_depth + DNA_concentration_ng_ul"
  
  # Add covariates if required 
  if (covariates == "cor_delivery") {
    base_formula <- paste0(base_formula, " + birth_deliverybirthcard_mode_binary")
  } else if (covariates == "cor_delivery_feeding") {
    base_formula <- paste0(base_formula, " + birth_deliverybirthcard_mode_binary + infant_ffq_ever_never_breastfed")
  }
  
  adon <- foreach(i = 1:length(Phenotype_list), .combine = rbind) %do% {
    Pheno = Phenotype_list[i]
    print(Pheno)
    Values <- as.vector(as_vector(phenotypes[[Pheno]]))
    
    # Filter samples with non-missing phenotype values
    r <- row.names(phenotypes)[!is.na(Values)]
    phenos2 <- phenotypes[r, ]
    
    # Filter NAs in covariates too
    vars_needed <- c(all.vars(as.formula(paste0("~", base_formula))), Pheno)
    phenos2 <- phenos2 %>% dplyr::select(all_of(vars_needed)) %>% tidyr::drop_na()
    
    # Align distance matrix
    r <- rownames(phenos2)
    distmat_cleaned <- Distance[r, r]
    
    # Construct formula
    FR <- as.formula(paste0("distmat_cleaned ~ ", base_formula, " + ", Pheno))
    
    # Run adonis2
    ad1 <- adonis2(FR, phenos2, permutations = perm, parallel = cores, na.action = na.fail, by="margin")
    
    # Extract phenotype result
    pheno_row <- which(rownames(ad1) == Pheno)
    residual_row <- which(rownames(ad1) == "Residual")
    
    if (length(pheno_row) == 0) {
      warning(paste("Phenotype", Pheno, "not found in adonis2 result."))
      adonis_results[i, ] <- c(Pheno, NA, NA, NA, NA)
    } else {adonis_results[i, ] <- c(Pheno, ad1$Df[pheno_row], ad1$F[pheno_row],
        ad1$R2[pheno_row], ad1$'Pr(>F)'[pheno_row], ad1$Df[residual_row]
      )
    }
  }
  
  return(adonis_results %>% drop_na())
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

Sample_metadata <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt")
Sample_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_phenos_07052025.txt")
Sample_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_07052025.txt")

Virus_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")
Virus_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_INFANTS_07052025.txt")
Virus_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_MOTHERS_07052025.txt")


##************************************************************************
# 1. PERMANOVA analysis - Unadjusted for delivery and feeding mode
##************************************************************************

# Load phenotypes of interest for associations
Phenotypes <- read.delim("9_ASSOCIATION_PHENOTYPES/Phenos_of_interest.txt")
Phenotypes <- as.character(Phenotypes$x)

# Add feeding mode after birth also for this analysis
# Exclude "birth_delivery_mode_CS"  and "birth_delivery_mode_simple" 
Phenotypes <- c(Phenotypes, "infant_birthcard_feeding_mode_after_birth")      
Phenotypes <- setdiff(Phenotypes, c("birth_delivery_mode_CS", "birth_delivery_mode_simple")
                      )
# Set parity as a binary phenotype (no if 0, yes if > 0 )
Sample_metadata_infants <- Sample_metadata_infants %>%
  mutate(
    parity_num = as.numeric(as.character(mother_birthcard_parity)),
    mother_birthcard_parity = if_else(parity_num > 0, "yes", "no")
  )

# Estimate the statistical significance using PERMANOVA
# Use the viral composition grouped by host
dist_infants <- as.matrix(vegdist(Abundance_table_infants_host, method = "bray"))

# Running PERMANOVA with the base model (unadjusted for feeding and delivery)
ResultsAdonis_base <- tibble()

for (Timepoint in unique(Sample_metadata_infants$Timepoint_categorical)) {
  print(Timepoint)
  
  # Subset data (per timepoint) and define phenotypes to use
  data <- Sample_metadata_infants %>% filter(Timepoint_categorical == Timepoint)
  rownames(data) <- data$NG_ID
  PhenosTest <- Phenotypes
  
  # Run adonis
  temp <- Call_adonis(data, Phenotype_list = PhenosTest, perm = 10000,
    Distance = dist_infants,covariates = "base") %>%
    mutate(Timepoint = Timepoint)
  
  # Add to results
  ResultsAdonis_base <- bind_rows(ResultsAdonis_base, temp)
}

# Add FDR correction
ResultsAdonis_base <- ResultsAdonis_base %>%
  mutate(FDR = p.adjust(`p-value`, method = "BH")) %>%
  arrange(FDR)

# Write results
write.table(ResultsAdonis_base, "9_ASSOCIATION_PHENOTYPES/RESULTS/ADONIS_Associations_per_Timepoint_by_MARGIN_INFANTS_unadj.txt",
            sep = "\t", row.names = F)

##************************************************************************
# 2. PERMANOVA analysis - Adjusted for delivery and feeding mode
##************************************************************************

# Running PERMANOVA with the adjusted model
ResultsAdonis_adj <- tibble()

for (Timepoint in unique(Sample_metadata_infants$Timepoint_categorical)) {
  print(Timepoint)
  
  # Subset data (per timepoint) and define phenotypes to use
  data <- Sample_metadata_infants %>% filter(Timepoint_categorical == Timepoint)
  rownames(data) <- data$NG_ID
  PhenosTest <- setdiff(Phenotypes,c("birth_deliverybirthcard_mode_binary","infant_ffq_ever_never_breastfed",
                                     "birth_delivery_mode_simple", "infant_birthcard_feeding_mode_after_birth", 
                                     "birth_delivery_mode_CS"))
  
  # Run adonis
  temp <- Call_adonis(data, Phenotype_list = PhenosTest, perm = 10000,
    Distance = dist_infants, covariates = "cor_delivery_feeding") %>%
    mutate(Timepoint = Timepoint)
  
  # Add to results
  ResultsAdonis_adj <- bind_rows(ResultsAdonis_adj, temp)
}

# Add FDR correction
ResultsAdonis_adj  <- ResultsAdonis_adj %>%
  mutate(FDR = p.adjust(`p-value`, method = "BH")) %>%
  arrange(FDR)

# Write results
write.table(ResultsAdonis_adj, "9_ASSOCIATION_PHENOTYPES/RESULTS/ADONIS_Associations_per_Timepoint_by_MARGIN_INFANTS_adj.txt",
            sep = "\t", row.names = F)


##************************************************************************
# 3. TEMPTED analysis - Associations across timepoints
##************************************************************************

########################
# 3.A. Running TEMPTED
########################

# Add NG_ID as rownames in Sample metadata
Abundance_table_infants <- Abundance_table_infants_host
rownames(Sample_metadata_infants) <- Sample_metadata_infants$NG_ID

# Get subjects with more than 3 timepoint available (NEXT_IDs with more than 3 NG_ID)
subject_counts <- Sample_metadata_infants %>%
  dplyr::count(NEXT_ID, name = "n_timepoints")

subjects_with_multiple_timepoints <- subject_counts %>%
  dplyr::filter(n_timepoints > 3) %>%
  dplyr::pull(NEXT_ID)

# Select only them from Sample metadata and abundance table
Sample_metadata_infants_TEMPTED <- Sample_metadata_infants %>%
  dplyr::filter(NEXT_ID %in% subjects_with_multiple_timepoints)

Abundance_table_infants_TEMPTED <- Abundance_table_infants[rownames(Abundance_table_infants) %in% 
                                                                     Sample_metadata_infants_TEMPTED$NG_ID, ]

# Remove vOTUs that are no longer present (after the filtering)
zero_cols <- which(colSums(Abundance_table_infants_TEMPTED == 0) == nrow(Abundance_table_infants_TEMPTED))
Abundance_table_infants_TEMPTED <- Abundance_table_infants_TEMPTED[, -zero_cols, drop = FALSE]

# Log-transform viral abundance tables (including pseudocount)
pseudocount_vOTU_infants <- min(Abundance_table_infants_TEMPTED[Abundance_table_infants_TEMPTED > 0]) / 2
Abundance_table_infants_TEMPTED_log <- Abundance_table_infants_TEMPTED
Abundance_table_infants_TEMPTED_log[Abundance_table_infants_TEMPTED_log == 0] <- pseudocount_vOTU_infants
Abundance_table_infants_TEMPTED_log <- log(Abundance_table_infants_TEMPTED_log)

# Order abundance table according to Sample metadata
Abundance_table_infants_TEMPTED_log  <- Abundance_table_infants_TEMPTED_log [match(rownames(Sample_metadata_infants_TEMPTED),
                                                                                   rownames(Abundance_table_infants_TEMPTED_log )), ]

# Selected only phenotypes of interest from metadata
Sample_metadata_infants_TEMPTED <- Sample_metadata_infants_TEMPTED[,c("NEXT_ID", "exact_age_days_at_collection",Phenotypes )]

# Run TEMPTED
results_tempted <- tempted_all(
  featuretable =  Abundance_table_infants_TEMPTED_log,         
  timepoint =  Sample_metadata_infants_TEMPTED$exact_age_days_at_collection,  
  subjectID = Sample_metadata_infants_TEMPTED$NEXT_ID,       
  threshold = 0.95,                                
  transform = "none",                              
  r = 5,                                           
  smooth=1e-8,                                   
  pct_ratio = 0.05,                                 
  pct_aggregate = 1                               
)

# Save results
saveRDS(results_tempted, file = "9_ASSOCIATION_PHENOTYPES/RESULTS/Results_tempted.rds")


##########################################
# 3.B. Association with temporal trajectories
##########################################

# Retrieve subject trajectories per timepoint (in days)
subject_trajectories <- results_tempted$metafeature_aggregate
colnames(subject_trajectories) <- c("Value", "NEXT_ID", "exact_age_days_at_collection", "PC")

# Reformat to have 1 column for PC
subject_trajectories_wide <- subject_trajectories %>%
  pivot_wider(names_from = PC, values_from = Value)

# Add variables to Sample metadata infants
Sample_metadata_infants <- Sample_metadata_infants %>%
  left_join(subject_trajectories_wide,by = c("NEXT_ID", "exact_age_days_at_collection"))

# Select viral features of interest (PC1-5)
Viral_features <- c("PC1", "PC2", "PC3", "PC4", "PC5")

# Run the association (LMM unadjusted)
results_TEMPTED_PCs_unadj <- associate_features_with_phenotypes(
  sample_metadata = Sample_metadata_infants,
  subject_id_col = "NG_ID",
  feature_data = Viral_features,
  phenotype_list = Phenotypes
)

results_TEMPTED_PCs_unadj <- results_TEMPTED_PCs_unadj %>%
  arrange(FDR)

# Run the association (LMM adjusted)
results_TEMPTED_PCs_adj <- associate_features_with_phenotypes(
  sample_metadata = Sample_metadata_infants,
  subject_id_col = "NG_ID",
  feature_data = Viral_features,
  phenotype_list = setdiff(Phenotypes, c("birth_deliverybirthcard_mode_binary","infant_ffq_ever_never_breastfed",
                                         "birth_delivery_mode_simple", "infant_birthcard_feeding_mode_after_birth")),
  covariates = c("birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed")
)

results_TEMPTED_PCs_adj <- results_TEMPTED_PCs_adj %>%
  arrange(FDR)

# Save results
write.table(results_TEMPTED_PCs_unadj, "9_ASSOCIATION_PHENOTYPES/RESULTS/Associations_TEMPTED_PCs_Phenotypes_INFANTS_LOG_NO_adj.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(results_TEMPTED_PCs_adj, "9_ASSOCIATION_PHENOTYPES/RESULTS/Associations_TEMPTED_PCs_Phenotypes_INFANTS_LOG_adj.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)



##************************************************************************
# 4. Plots
##************************************************************************

##########################################
# A. PERMANOVA results
##########################################

# First, modify PERMANOVA results for plotting
ResultsAdonis_base_plot <- ResultsAdonis_base %>%
  mutate(p_value = as.numeric(`p-value`),
         R2 = as.numeric(R2),
         logP = -log10(p_value),
         Timepoint = factor(Timepoint, levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12"))
  ) %>%
  filter(p_value < 0.05)   

# Set order of phenotypes
phenotype_order <- c("mother_education_p18","mother_health_smoked_one_whole_year_p18",
  "mother_ffq_pattern_processed_foods","mother_ffq_pattern_meat_and_fat", "mother_ffq_aMED_score",
  "mother_birthcard_parity", "mother_deliveryhealth_preg_complaint_uti",
  "mother_deliveryhealth_preg_complaint_flu_rti", "mother_birthcardself_gestational_age_weeks",
  "infant_growth_birth_weight_kg", "infant_health_eczema_diagnosis_relaxed",
  "infant_health_eczema_diagnosis_strict", "infant_birthcard_feeding_mode_after_birth",
  "infant_deliverybirthcard_meconium_amniotic_fluid",
  "birth_deliverybirthcard_place_delivery_simple",
  "birth_deliverybirthcard_mode_binary")

ResultsAdonis_base_plot$Phenotype <- factor(ResultsAdonis_base_plot$Phenotype,
  levels = phenotype_order)
ResultsAdonis_base_plot <- ResultsAdonis_base_plot %>%
  filter(!is.na(Phenotype))

# Generate the plot
pdf("9_ASSOCIATION_PHENOTYPES/Plots/PERMANOVA_summary.pdf", width = 9, height = 8.5)
ggplot(ResultsAdonis_base_plot, aes(x = Timepoint, y = Phenotype, 
                              size = R2, color = logP, group = Phenotype)) +
  geom_line(color = "grey80", size = 0.4, alpha = 0.8) +  
  geom_point(alpha = 0.9) +
  scale_color_gradientn(colours = c("#E9E5F0", "#C3BED9", "#9B91C9", "#6E6BAF", "#3F3B73"),
                        name = expression(-log[10](p))) +
  scale_size_continuous(range = c(2, 8), name = expression(R^2)) +
  labs(x = "Timepoint", y = "Phenotype") +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text.x = element_text(size = 14, color = "black"),
    axis.text.y = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14))
dev.off()

pdf("9_ASSOCIATION_PHENOTYPES/Plots/PERMANOVA_summary2.pdf", width = 9, height = 7.9)
ggplot(ResultsAdonis_base_plot, aes(x = Timepoint, y = Phenotype)) +
  geom_tile(aes(fill = logP), color = "white", linewidth = 0.4, alpha = 0.85) +
  geom_point(aes(size = R2, color = logP), alpha = 0.95, stroke = 0.3) +
  scale_fill_gradientn(colours = c("#F6E9DC", "#DCC7E3", "#A999CA", "#6C5A9D", "#352C63"),
    name = expression(-log[10](p))) +
  scale_color_gradientn(colours = c("#EBD3BC", "#C7A7D8", "#8E74B7", "#59457E", "#241E4E"),
    name = expression(-log[10](p))) +
  scale_size_continuous(range = c(3, 9), name = expression(R^2)) +
  labs(x = "Timepoint", y = "Phenotype") +
  theme_minimal(base_size = 16) +
  theme(panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.75),
    axis.title.x = element_text(size = 15, margin = margin(t = 10)),
    axis.title.y = element_text(size = 15, margin = margin(r = 10)),
    axis.text.x = element_text(size = 13, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    legend.title = element_text(size = 14, face = "bold"),
    legend.text = element_text(size = 12),
    legend.key.height = unit(0.6, "cm"),
    legend.key.width = unit(0.6, "cm"),
    plot.margin = margin(10, 20, 10, 10)
  )
dev.off()

##########################################
# B. TEMPTED results (trajectory)
##########################################

# Generate subject trajectory plots for significant results 
subject_trajectories <- results_tempted$metafeature_aggregate

TEMPTED_delivery_grouping <- unique(Sample_metadata_infants_TEMPTED[,c("NEXT_ID",
                                                                       "birth_deliverybirthcard_mode_binary")])
TEMPTED_feeding_grouping <- unique(Sample_metadata_infants_TEMPTED[,c("NEXT_ID",
                                                                      "infant_birthcard_feeding_mode_after_birth")])
TEMPTED_parity_grouping <- unique(Sample_metadata_infants_TEMPTED[,c("NEXT_ID",
                                                                     "mother_birthcard_parity")])


pdf("9_ASSOCIATION_PHENOTYPES/Plots/TEMPTED_Trajectory_by_Delivery.pdf", width = 12.5, height = 3)
plot_metafeature(subject_trajectories, TEMPTED_delivery_grouping) + 
  xlab("Days of Life") +
  ylab("TEMPTED component") +
  scale_color_manual(values = c("VG" = "#B6D7A8", "CS" = "#A2C2E5")) +
  scale_fill_manual(values = c("VG" = "#B6D7A8", "CS" = "#A2C2E5")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
  )
dev.off()

TEMPTED_feeding_grouping$infant_birthcard_feeding_mode_after_birth <- factor(
  TEMPTED_feeding_grouping$infant_birthcard_feeding_mode_after_birth,
  levels = c("BF", "MF", "FF")
)

pdf("9_ASSOCIATION_PHENOTYPES/Plots/TEMPTED_PC5_Trajectory_by_Feeding.pdf", width = 13, height = 3)
plot_metafeature(subject_trajectories, TEMPTED_feeding_grouping) + 
  xlab("Days of Life") +
  ylab("TEMPTED component") +
  scale_color_manual(values = c("#A57BBE", "#D97C5A", "#D8C9A9")) +
  scale_fill_manual(values = c("#A57BBE", "#D97C5A", "#D8C9A9")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
  )
dev.off()

TEMPTED_parity_grouping$mother_birthcard_parity <- factor(
  TEMPTED_parity_grouping$mother_birthcard_parity,
  levels = c("no", "yes"))

pdf("9_ASSOCIATION_PHENOTYPES/Plots/TEMPTED_Trajectory_by_Parity.pdf", width = 12.5, height = 3)
plot_metafeature(subject_trajectories, TEMPTED_parity_grouping[, c("NEXT_ID", "mother_birthcard_parity")]) +
  xlab("Days of Life") +
  ylab("TEMPTED component") +
  scale_color_manual(values = c("no"  = "#3C5488","yes" = "#E64B35")) +
  scale_fill_manual(values = c("no"  = "#3C5488", "yes" = "#E64B35")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    )
dev.off()


pdf("9_ASSOCIATION_PHENOTYPES/Plots/TEMPTED_PC3_Temporal_Loadings_by_Delivery.pdf", width = 2.1, height = 3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$birth_deliverybirthcard_mode_binary),],
       aes(x = as.factor(birth_deliverybirthcard_mode_binary),
           y = PC3, fill = as.factor(birth_deliverybirthcard_mode_binary))) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = "Delivery Mode", y = "PC3 (TEMPTED component)") +
  scale_fill_manual(values = c("VG" = "#B6D7A8", "CS" = "#A2C2E5")) +
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

Sample_metadata_infants$infant_birthcard_feeding_mode_after_birth <- factor(
  Sample_metadata_infants$infant_birthcard_feeding_mode_after_birth,
  levels = c("BF", "MF", "FF")
)

pdf("9_ASSOCIATION_PHENOTYPES/Plots/TEMPTED_PC5_Temporal_Loadings_Feeding.pdf", width = 2.1, height = 3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$infant_birthcard_feeding_mode_after_birth),],
       aes(x = as.factor(infant_birthcard_feeding_mode_after_birth),
           y = PC5,
           fill = as.factor(infant_birthcard_feeding_mode_after_birth))) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = "Feeding mode", y = "PC5 (TEMPTED component)") +
  scale_fill_manual(values = c("#A57BBE", "#D97C5A", "#D8C9A9")) +
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

pdf("9_ASSOCIATION_PHENOTYPES/Plots/TEMPTED_PC2_Temporal_Loadings_Parity.pdf", width = 2.1, height = 3.2)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$mother_birthcard_parity),],
       aes(x = as.factor(mother_birthcard_parity),
           y = PC2,
           fill = as.factor(mother_birthcard_parity))) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = "Parity", y = "PC2 (TEMPTED component)") +
  scale_fill_manual(values = c("no"  = "#3C5488","yes" = "#E64B35")) +
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

##########################################
# C. TEMPTED results (subject loadings) 
##########################################

# Get subject loadings
subject_loadings <- results_tempted$A_hat

# Convert rownames to a column so we can join by NEXT_ID
subject_loadings <- as.data.frame(subject_loadings)
subject_loadings$NEXT_ID <- rownames(subject_loadings)

# Add the two phenotype columns from Sample_metadata_infants
subject_loadings <- subject_loadings %>%
  left_join(Sample_metadata_infants %>%
              dplyr::select(NEXT_ID, birth_deliverybirthcard_mode_binary, infant_birthcard_feeding_mode_after_birth),
            by = "NEXT_ID") %>% 
  distinct()

# Generate a 2-dimensional representation of infants according to phenotypes
pdf("9_ASSOCIATION_PHENOTYPES/Plots/TEMPTED_PC1_PC2_by_Delivery.pdf", width = 4, height = 4)
p <- ggplot(subject_loadings[!is.na(subject_loadings$birth_deliverybirthcard_mode_binary), ],
  aes(x = PC1, y = PC3, color = as.factor(birth_deliverybirthcard_mode_binary))) +
  geom_point(size = 2, alpha = 0.8) +
  labs(x = "Component 1", y = "Component 3") +
  scale_color_manual(values = c("VG" = "#B6D7A8", "CS" = "#A2C2E5")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )

ggMarginal(p, type = "density", groupColour = TRUE, groupFill = TRUE,
  alpha = 0.5,size = 4)
dev.off()

subject_loadings$infant_birthcard_feeding_mode_after_birth <- factor(
  subject_loadings$infant_birthcard_feeding_mode_after_birth,
  levels = c("BF", "MF", "FF")
)

pdf("9_ASSOCIATION_PHENOTYPES/Plots/TEMPTED_PC2_PC5_by_Feeding_mode.pdf", width = 4, height = 4)
p <- ggplot(subject_loadings[!is.na(subject_loadings$infant_birthcard_feeding_mode_after_birth), ],
  aes(x = PC2, y = PC5, color = as.factor(infant_birthcard_feeding_mode_after_birth))) +
  geom_point(size = 2, alpha = 0.8) + labs(x = "Component 2", y = "Component 5") +
  scale_color_manual(values =  c("#A57BBE", "#D97C5A", "#D8C9A9")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 13, color = "black"),
    axis.title = element_text(size = 15),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )

ggMarginal(p,type = "density", groupColour = TRUE, groupFill = TRUE,
  alpha = 0.5, size = 4)
dev.off()

##########################################
# D. TEMPTED results (other phenos) 
##########################################

# Generate plot with assocations results of phenotypes to TEMPTED PCs (after adjusting)
# Add direction of effect to results
phenos_df <- results_TEMPTED_PCs_adj[ results_TEMPTED_PCs_adj$FDR < 0.05,]
phenos_df <- phenos_df %>%
  mutate(Sign = ifelse(Estimate > 0, "Positive", "Negative"))

# Order for plotting
phenos_df$Phenotype <- factor(phenos_df$Phenotype,levels = rev(unique(phenos_df$Phenotype)))
phenos_df$Feature <- factor(phenos_df$Feature,levels = c("PC1","PC2","PC3","PC4","PC5"))

# Plot
pdf("9_ASSOCIATION_PHENOTYPES/Plots/TEMPTED_PC_Phenotype_Associations_FDR.pdf",
    width = 8.5, height = 3)
ggplot(phenos_df, aes(x = Feature, y = Phenotype)) + geom_point(aes(size = FDR, fill = Sign),
             shape = 21, color = "black", alpha = 0.8) +
  scale_fill_manual(values = c("Positive" = "#C27BA0", "Negative" = "#6FA8DC")) +
  scale_size_continuous(range = c(7, 3)) +
  labs(x = "Component", y = NULL, fill = "Direction", size = "FDR") +
  theme_bw() +
  coord_fixed() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(size = 12, color = "black"),
    axis.text.y = element_text(size = 12, color = "black"),
    axis.title.x = element_text(size = 13, margin = margin(t = 8)),
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 11),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75)
  )
dev.off()

pdf("9_ASSOCIATION_PHENOTYPES/Plots/TEMPTED_PC_Phenotype_Associations_FDR2.pdf",
    width = 8.5, height = 3)
ggplot(phenos_df, aes(x = Feature, y = Phenotype, fill = `t value`)) +
  geom_tile(color = "white", linewidth = 0.6) +
  scale_fill_gradient2(low = "#6FA8DC", mid = "white", high = "#C27BA0",
    midpoint = 0, name = "T-value") +
  theme_minimal(base_size = 14) +
  labs(x = "Component", y = NULL) +
  coord_fixed() +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_text(size = 12, color = "black", angle = 0, vjust = 0.5),
    axis.text.y = element_text(size = 12, color = "black"),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5)
  )
dev.off()

##*************
# Save output
#**************
write.table(Sample_metadata_infants,"Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 



