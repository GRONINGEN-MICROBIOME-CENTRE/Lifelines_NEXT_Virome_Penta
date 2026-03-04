################################################################################
##### LL-NEXT: Viral persistence vs human phenotypes
### Author(s): Asier Fernández-Pato
### Last updated: 11th December, 2025
################################################################################

#****************
# Load libraries
#****************
library(dplyr)
library(tidyr)
library(data.table)
library(ggplot2)
library(lmer)
library(lmerTest)
library(scales)


#****************
# Define functions
#****************

# Function to estimate the association between individual-level features and phenotypes of interest
test_phenotype_associations_lm <- function(data, outcome_var,phenotypes,
    covariates = c("clean_reads_FQ_2", "DNA_concentration_ng_ul")
) {
  
  results <- list()
  
  for (phenotype in phenotypes) {
    
    # Build RHS of models
    rhs_null <- c(covariates)
    rhs_full <- c(covariates, phenotype)
    
    # Build formulas
    formula_null <- as.formula(paste(outcome_var, "~", paste(rhs_null, collapse = " + ")))
    formula_full <- as.formula(paste(outcome_var, "~", paste(rhs_full, collapse = " + ")))
    
    # Filter out rows with missing phenotype (same as LMM code)
    df_sub <- data[!is.na(data[[phenotype]]), ]
    
    # Fit models on the filtered df_sub
    model_null <- lm(formula_null, data = df_sub)
    model_full <- lm(formula_full, data = df_sub)
    
    # Likelihood Ratio Test (LRT)
    lrt_p <- anova(model_null, model_full)$`Pr(>F)`[2]
    
    # Coefficients of full model
    coef_info <- summary(model_full)$coefficients
    
    # Extract rows for the phenotype term
    phenotype_rows <- grep(phenotype, rownames(coef_info), value = TRUE)
    
    for (term in phenotype_rows) {
      results[[length(results) + 1]] <- tibble(
        phenotype = phenotype,
        term = term,
        estimate = coef_info[term, "Estimate"],
        std_error = coef_info[term, "Std. Error"],
        p_value = lrt_p
      )
    }
  }
  
  bind_rows(results) %>%
    mutate(FDR = p.adjust(p_value, method = "BH"))
}

# Function to estimate the association between sample-level features and phenotypes of interest
test_phenotype_associations_lmm <- function(df, outcome_var, phenotypes, subject_col,
                                            covariates = c("clean_reads_FQ_2", "DNA_concentration_ng_ul", "Timepoint_categorical")) {
  results <- list()
  
  for (phenotype in phenotypes) {
    
    # Construct formulas
    rhs_null <- c(covariates)
    rhs_full <- c(covariates, phenotype)
    
    formula_null <- as.formula(paste(outcome_var, "~", paste(rhs_null, collapse = " + "), "+ (1|", subject_col, ")"))
    formula_full <- as.formula(paste(outcome_var, "~", paste(rhs_full, collapse = " + "), "+ (1|", subject_col, ")"))
    
    # Filter out rows with missing phenotype
    df_sub <- df[!is.na(df[[phenotype]]), ]
    
    # Fit models
    model_null <- lmer(formula_null, data = df_sub, REML = FALSE)
    model_full <- lmer(formula_full, data = df_sub, REML = FALSE)
    
    # LRT p-value
    p_val <- anova(model_full, model_null)["model_full", "Pr(>Chisq)"]
    
    # Coefficients
    coef_info <- summary(model_full)$coefficients
    phenotype_rows <- grep(phenotype, rownames(coef_info), value = TRUE)
    
    for (term in phenotype_rows) {
      res <- tibble(
        phenotype = phenotype,
        term = term,
        estimate = coef_info[term, "Estimate"],
        std_error = coef_info[term, "Std. Error"],
        p_value = p_val
      )
      results[[length(results) + 1]] <- res
    }
  }
  
  bind_rows(results) %>%
    mutate(FDR = p.adjust(p_value, method = "BH"))
}


#****************
# Load data
#****************

# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_Analysis/VIRUSES/")

#Load metadata and abundance table for the LL-NEXT samples 

# Read metadata tables 
Sample_metadata <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt")
Sample_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_phenos_07052025.txt")
Sample_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Sample_metadata_MOTHERS_phenos_15012025.txt")

# Read abundance table
Abundance_table <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_07052025.txt")
Abundance_table_infants <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_INFANTS_07052025.txt")
Abundance_table_mothers <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_MOTHERS_07052025.txt")

Viral_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")
Viral_metadata_infants <- Viral_metadata[Viral_metadata$Virus_ID %in% rownames(Abundance_table_infants),]
Viral_metadata_mothers <- Viral_metadata[Viral_metadata$Virus_ID %in% rownames(Abundance_table_mothers),]

# Load infant persistence data
infant_viral_persistence <- read.delim("11_PERSISTENCE/LLNEXT_Infant_Viral_Persistance_data.txt")
infant_complete_sampling <- as.vector(read.delim("11_PERSISTENCE/Infants_complete_sampling.txt", 
  header = FALSE)[,1])
  
#**************************************************************
# 1. Selection of infants with near-complete sampling
#**************************************************************

# For a fair comparison, associations are only done on a subset of infants with near-complete longitudinal sampling
# We subset infant_viral_persistence and Sample_metadata_infants accordingly
infant_viral_persistence_complete <- infant_viral_persistence[infant_viral_persistence$NEXT_ID %in% 
                                                                infant_complete_sampling, ]

Sample_metadata_infants_complete <- Sample_metadata_infants[Sample_metadata_infants$NEXT_ID %in%
                                                              infant_complete_sampling,]


#**************************************************************
# 2. Associations of human phenotypes with infant vOTU-level persistence
#**************************************************************

# Select a subset of 10 phenotypes of interest for this analysis
phenotypes_to_test <- c("birth_deliverybirthcard_mode_binary", "birth_deliverybirthcard_place_delivery_simple",
                        "family_pets_any", "infant_health_eczema_diagnosis_relaxed", "infant_health_eczema_diagnosis_strict", 
                        "mother_birthcard_parity", "infant_health_food_allergy", "infant_growth_birth_weight_kg",  
                        "mother_birthcard_age_at_delivery", "infant_ffq_ever_never_breastfed" )

# Identify categorical and continuous variables
categorical_vars <- names(Sample_metadata_infants_complete)[
  sapply(Sample_metadata_infants_complete, function(x)
    is.character(x) || (is.numeric(x) && length(unique(na.omit(x))) < 10))
]

continuous_vars <- setdiff(names(Sample_metadata_infants_complete), categorical_vars)

# Set as factors categorical variables
Sample_metadata_infants_complete <- Sample_metadata_infants_complete %>%
  mutate(across(all_of(categorical_vars), as.factor))

# Reorder birth_deliverybirthcard_place_delivery_simple
Sample_metadata_infants_complete$birth_deliverybirthcard_place_delivery_simple <- 
  factor(Sample_metadata_infants_complete$birth_deliverybirthcard_place_delivery_simple,
         levels = c("hospital","home")) 

# Set parity as a binary phenotype (no if 0, yes if > 0 )
Sample_metadata_infants_complete <- Sample_metadata_infants_complete %>%
  mutate(
    parity_num = as.numeric(as.character(mother_birthcard_parity)),
    mother_birthcard_parity = if_else(parity_num > 0, "yes", "no")
  )

###############
# A. Proportion
###############

# Add phenotypes of interest to infant_viral_persistence_complete
phenotypes_per_infant <- Sample_metadata_infants_complete %>%
  group_by(NEXT_ID) %>%
  dplyr::slice(1) %>%  # keep first only 1 (all same phenotypes)
  ungroup() %>%
  select(NEXT_ID, all_of(phenotypes_to_test),          
    clean_reads_FQ_2, DNA_concentration_ng_ul             
  )

infant_viral_persistence_complete <- infant_viral_persistence_complete %>%
  left_join(phenotypes_per_infant, by = "NEXT_ID")

# Reorder birth_deliverybirthcard_place_delivery_simple
infant_viral_persistence_complete$birth_deliverybirthcard_place_delivery_simple <- 
  factor(infant_viral_persistence_complete$birth_deliverybirthcard_place_delivery_simple,
         levels = c("hospital","home"))

# Run the association (LM - only 1 observation per infant)
lm_results_infant <- test_phenotype_associations_lm(
  data = infant_viral_persistence_complete,
  outcome_var = "Prop_Persistent_vOTUs",
  phenotypes = phenotypes_to_test,
  covariates = c("clean_reads_FQ_2", "DNA_concentration_ng_ul") 
)

# Run the association after correcting for feeding mode
lm_results_infant_adj <- test_phenotype_associations_lm(
  data = infant_viral_persistence_complete,
  outcome_var = "Prop_Persistent_vOTUs",
  phenotypes = setdiff(phenotypes_to_test, c("infant_ffq_ever_never_breastfed",
                                             "birth_deliverybirthcard_mode_binary")),
  covariates = c("clean_reads_FQ_2", "DNA_concentration_ng_ul",
                 "infant_ffq_ever_never_breastfed", "birth_deliverybirthcard_mode_binary") 
)

# Order results by FDR
lm_results_infant <- lm_results_infant %>%
  arrange(FDR)
lm_results_infant_adj <- lm_results_infant_adj %>%
  arrange(FDR)

###############
# B. Abundance
###############

# Test association with "Persistent_Abundance_Fraction"
lmm_results_infant <- test_phenotype_associations_lmm(
  df = Sample_metadata_infants_complete,          
  outcome_var = "Persistent_Abundance_Fraction", 
  phenotypes = phenotypes_to_test,
  subject_col = "NEXT_ID",
  covariates = c("clean_reads_FQ_2", "DNA_concentration_ng_ul", "Timepoint_categorical")
)

# Test association with "Persistent_Abundance_Fraction" after adjusting for feeding mode
lmm_results_infant_adj  <- test_phenotype_associations_lmm(
  df = Sample_metadata_infants_complete,          
  outcome_var = "Persistent_Abundance_Fraction", 
  phenotypes = setdiff(phenotypes_to_test, c("infant_ffq_ever_never_breastfed",
                                             "birth_deliverybirthcard_mode_binary")),
  subject_col = "NEXT_ID",
  covariates = c("clean_reads_FQ_2", "DNA_concentration_ng_ul", "Timepoint_categorical",
                 "infant_ffq_ever_never_breastfed", "birth_deliverybirthcard_mode_binary")
)

# Order results by FDR
lmm_results_infant <- lmm_results_infant %>%
  arrange(FDR)
lmm_results_infant_adj <- lmm_results_infant_adj %>%
  arrange(FDR)

###############
# C. PLOTS 
###############

# Infant persistent proportion vs delivery place
pdf('11_PERSISTENCE/Plots/Prop_persistent_vOTUs_delivery_place.pdf', width=2.8, height=3.2)
ggplot(infant_viral_persistence_complete[!is.na(infant_viral_persistence_complete$birth_deliverybirthcard_place_delivery_simple),], 
       aes(x = birth_deliverybirthcard_place_delivery_simple, y = Prop_Persistent_vOTUs, fill = birth_deliverybirthcard_place_delivery_simple)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Delivery place', y = '% Persistent vOTUs') +
  scale_fill_manual(values = c("#D8B05B", "#6B7AA1")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, by = 0.2),
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

# Infant relab persistent vOTUs vs feeding mode
pdf('11_PERSISTENCE/Plots/Relab_persistent_vOTUs_delivery_place.pdf', width=2.9, height=3.2)
ggplot(Sample_metadata_infants_complete[!is.na(Sample_metadata_infants_complete$birth_deliverybirthcard_place_delivery_simple),], 
       aes(x = birth_deliverybirthcard_place_delivery_simple, y = Persistent_Abundance_Fraction, fill = birth_deliverybirthcard_place_delivery_simple)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Delivery place', y = 'Persistent vOTU abundance') +
  scale_fill_manual(values = c("#D8B05B", "#6B7AA1")) +
  scale_y_continuous(labels = percent_format(accuracy = 1), breaks = seq(0, 1, by = 0.25), limits = c(0, 1),
                     expand = expansion(mult = c(0.1, 0.25))
  ) +
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

#****************
# Write results
#****************
write.table(lm_results_infant,"11_PERSISTENCE/Results/Prop_Persistent_vOTUs_association_INFANTS.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(lm_results_infant_adj,"11_PERSISTENCE/Results/Prop_Persistent_vOTUs_association_adj_INFANTS.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(lmm_results_infant,"11_PERSISTENCE/Results/Persistent_Abundance_Fraction_association_INFANTS.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(lmm_results_infant_adj,"11_PERSISTENCE/Results/Persistent_Abundance_Fraction_association_adj_INFANTS.txt", sep = "\t", row.names = F, quote = FALSE) 


