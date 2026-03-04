################################################################################
##### LL-NEXT:  DGR analysis -  General characterization
### Author(s): Asier Fernández-Pato
### Last updated: 13th December, 2025
################################################################################

#****************
# Load libraries
#****************
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(Biostrings)
library(ggplot2)
library(scales)


#****************
# Define functions
#****************
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


#***********************************************
# 1. DGRs: General characterization
#***********************************************

# Read DGR identification results
DGR_results <- read.delim("5_DGR_ANALYSIS/All_RT_1k_DGR_detection_filtered.tsv", check.names = F)
DGR_target_results <- read.delim("5_DGR_ANALYSIS/All_RT_1k_DGR_target_detection_filtered.tsv", check.names = F)

#####################
# --- Filtering ----
#####################

# Filter those DGRs with less than 75% of Adenine mutation rate
DGR_results_out <-DGR_results[DGR_results$`Max A/T bias` < 0.75,]
DGR_results <- DGR_results[DGR_results$`Max A/T bias` >= 0.75,]

# Filter those DGRs also from the targets results 
# For this, add A/T bias to the targets file (matching to A/T bias in DGR_results also possible - highly similar)
DGR_target_results$`Max A/T bias` <- sapply(strsplit(DGR_target_results$`Adjusted bias (A;T;C;G)`, ";"), function(x) {
  counts <- as.numeric(x)
  A <- counts[1]
  T <- counts[2]
  max(A, T) / sum(counts)
})

DGR_target_results <- DGR_target_results[DGR_target_results$`Max A/T bias` >= 0.75,]
DGR_target_results$Virus_ID <- sub("_[^_]*$", "", sub(".*:", "", DGR_target_results[[1]]))
colnames(DGR_target_results) [1] <- "DGR_code"

# Filter out systems with different RTs linked to same target/VR region. Keep the systems with the closest RT to the target gene
# For this, first select cases of different RTs linked to same target/VR region
dup_candidates <- DGR_target_results[
  duplicated(DGR_target_results[, c("Target", "VR start", "VR end")]) |
    duplicated(DGR_target_results[, c("Target", "VR start", "VR end")], fromLast = TRUE),
]

# Estimate distance to target gene
dup_candidates$RT_distance <- abs(dup_candidates$`Target start` - dup_candidates$`VR start`)

# Keep only the closest RT for each group of duplicates
filtered_dups <- dup_candidates %>%
  group_by(Target, `VR start`, `VR end`) %>%
  slice_min(RT_distance, with_ties = FALSE) %>%
  ungroup()

# Generate the filtered dataset (#2548 vOTUs with DGRs)
DGR_target_results <- rbind(
  DGR_target_results[!(DGR_target_results$DGR_code %in% dup_candidates$DGR_code), ],
  filtered_dups[, colnames(DGR_target_results)]
)

# Note: We see some RTs targeting different genes at the same coordinates (both strands). We keep these cases.

###########################
# --- Characterization ----
###########################

# Retrieve the DGR+ phages and targets
DGR_phages <- Viral_metadata[Viral_metadata$DGR == "Yes", "Virus_ID"] #2548
DGR_targets <- unique(DGR_target_results$Target) #2930

# Estimate total number of DGRs
length(unique(DGR_target_results$DGR_code)) #2,587

# Estimate number of vOTUs with >1 DGR system (vOTUs with multiple RTs) (n=39 with 2)
phages_with_multiple_DGRs <- DGR_target_results %>%
  group_by(Virus_ID) %>%
  summarise(n_RT_genes = n_distinct(`DGR_code`)) %>%
  filter(n_RT_genes > 1)

# Estimate number of DGRs with > 1 target (same RT, multiple different targets)
dgr_target_counts <- aggregate(Target ~ DGR_code, data = DGR_target_results, function(x) length(unique(x)))
dgrs_with_multiple_targets <- dgr_target_counts[dgr_target_counts$Target > 1, ] # n= 336 (6 with 3 targets)
virus_ids_with_multi_target_DGRs <- unique(DGR_target_results$Virus_ID[DGR_target_results$DGR_code %in% 
                                                                         dgrs_with_multiple_targets$DGR_code])

# Estimate their relative abundance and proportion per sample
Sample_metadata$relab_DGR <- 100*(colSums(Abundance_table[rownames(Abundance_table) %in% DGR_phages,]) / 
                                    colSums(Abundance_table))
Sample_metadata$proportion_DGR <- 100*(colSums(Abundance_table[rownames(Abundance_table) %in% DGR_phages, ] > 0) /
                                         colSums(Abundance_table > 0))

# Add to Sample metadata infants and mothers
Sample_metadata_infants <- Sample_metadata_infants %>%
  left_join(Sample_metadata %>% select(NG_ID, relab_DGR, proportion_DGR),
            by = "NG_ID")
Sample_metadata_mothers <- Sample_metadata_mothers %>%
  left_join(Sample_metadata %>% select(NG_ID, relab_DGR, proportion_DGR),
            by = "NG_ID")

#***********************************************
# 2. Annotation of target genes
#***********************************************

# Read target annotation results (PHROGs) (done on 1,800,152 proteins)
HMMER_result_PHROGs <- fread("5_DGR_ANALYSIS/HMMER_processed_table.tsv", check.names = FALSE)
#colnames(HMMER_result_PHROGs) [colnames(HMMER_result_PHROGs) == "Description"] <- "Protein_Info"

# Read HMMER metadata
metadata <- read.delim("5_DGR_ANALYSIS/phrog_annot_v4.tsv", check.names = F)
metadata$phrog <- paste0("phrog_", metadata$phrog)

# Combine HMMER results and metadata to get the annotations for each protein (order by E-value)
Protein_annotations <- left_join(HMMER_result_PHROGs, metadata,  by = c("Match" = "phrog"))
Protein_annotations <- Protein_annotations[order(Protein_annotations$E_value), ]
Protein_annotations <- Protein_annotations[, c("Protein_ID","Match", "E_value", "Description", "annot", "category")]
colnames(Protein_annotations) <- c("Protein_ID","Match", "E_value", "Protein_Info", "DESCRIPTION", "CATEGORY")

# Retrieve annotations corresponding to target genes 
# Keep only top annotation (min E-value)
Annotation_DGR_targets <- Protein_annotations %>%
  filter(Protein_ID %in% DGR_targets) %>%
  group_by(Protein_ID) %>%
  slice_min(order_by = E_value, with_ties = FALSE) %>%
  ungroup()
table(Annotation_DGR_targets$CATEGORY)

# Add annotation column to DGR target results
DGR_target_results <- DGR_target_results %>%
  left_join(
    Annotation_DGR_targets %>%
      select(Protein_ID, DESCRIPTION, CATEGORY),
    by = c("Target" = "Protein_ID")
  ) %>%
  dplyr::rename(
    Annotation = DESCRIPTION,
    DGR_target_category = CATEGORY
  )

DGR_target_results$`Target product` <- NULL

# Explore if any DGR targets are annotated as ani-defense systems
antidefense_systems <- read.delim("6_ANTIDEFENSE_ANALYSIS/Present_vOTUs_proteins_defense_finder_systems.tsv", 
                                  header = T, check.names = F)
antidefense_systems <- distinct(antidefense_systems) #Get unique systems
all_antidefense_proteins <- unlist(strsplit(antidefense_systems$protein_in_syst, ","))
length(which(DGR_targets %in% all_antidefense_proteins)) #0

# Plot target annotation results (Annotation_DGR_targets$CATEGORY)
category_counts <- Annotation_DGR_targets %>%
  count(CATEGORY) %>%
  mutate(Percent = 100 * n / sum(n))

# Define colors (expandable palette)
category_colors <- c("#80CBC4","#A6CEE3", "#FDBF6F", "#D7CCC8")

pdf("5_DGR_ANALYSIS/Plots/Functional_categories_piechart.pdf", width = 4.2, height = 4.2)
ggplot(category_counts, aes(x = "", y = Percent, fill = CATEGORY)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = category_colors) +
  labs(x = NULL, y = NULL, fill = "Functional category") +
  theme_void(base_size = 14) +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 14),
    legend.text  = element_text(size = 12)
  )
dev.off()

#***********************************************
# 3. Association with shared vOTUs
#***********************************************

# Load list of shared vOTUs
shared_vOTUs <- read.delim("10_STRAIN_TRANSMISSION/shared_vOTUs.txt")
shared_vOTUs <- shared_vOTUs$x

# Add to infant viral metadata
Viral_metadata_infants$Sharing <- ifelse(Viral_metadata_infants$Virus_ID %in% shared_vOTUs, 1, 0)
Viral_metadata_infants$Sharing <-
  factor(ifelse(Viral_metadata_infants$Virus_ID %in% shared_vOTUs, "Yes", "No"))

# Test if maternally shared vOTUs are more likely to encode DGRs
# Test if both adjusting and not for bacterial host
Viral_metadata_infants$DGR <- factor(Viral_metadata_infants$DGR,levels = c("No", "Yes"))

DGR_maternal_sharing_association <- glm(DGR ~ Sharing + Mean_abundance_mothers,
  family = binomial, data = Viral_metadata_infants)
summary(DGR_maternal_sharing_association)

DGR_maternal_sharing_association_adj <- glmer(DGR ~ Sharing + Mean_abundance_mothers + (1 | Bacterial_genus_host),
                                        family = binomial, data = Viral_metadata_infants)
summary(DGR_maternal_sharing_association_adj)

# Generate tables with the number of Shared vOTUs and DGR_+ vOTUs per host
table_DGR_host <- table(Viral_metadata_infants$Bacterial_genus_host[Viral_metadata_infants$DGR == "Yes"])
table_sharing_host <- table(Viral_metadata_infants$Bacterial_genus_host[Viral_metadata_infants$Sharing == "Yes"])
prop_Bacteroides_DGR <- table_DGR_host["Bacteroides"] / sum(table_DGR_host)
prop_Bacteroides_Sharing <- table_sharing_host["Bacteroides"] / sum(table_sharing_host)


# Test if maternally shared vOTUs are more likely to persist in the infant gut (GLMM)
Viral_metadata_infants$Persistent <- factor(Viral_metadata_infants$Persistent,levels = c("No", "Yes"))

Persistence_maternal_sharing_association <- glmer( Persistent ~ Sharing + Mean_abundance_infants + (1 | Bacterial_genus_host),
  data = Viral_metadata_infants,family = binomial)

summary_coef <- coef(summary(Persistence_maternal_sharing_association))
p_value_sharing <- summary_coef["SharingYes", "Pr(>|z|)"]

################
# --- Plots ----
################

Viral_metadata_infants_plot <- Viral_metadata_infants %>%
  mutate(
    Persistent = case_when(
      Persistent == "Yes" ~ "Persistent",
      Persistent == "No"  ~ "Non-Persistent",
      TRUE ~ Persistent
    ),
    Sharing = case_when(
      Sharing == "Yes" ~ "Shared",
      Sharing == "No"  ~ "Not Shared",
      TRUE ~ Sharing
    )
  )
  
# Generate barplots for the association between DGR presence / Persistence and maternal sharing
pdf('11_PERSISTENCE/Plots/vOTU_persistence_by_sharing.pdf', width = 2.6, height = 3.6)
ggplot(Viral_metadata_infants_plot, aes(x = factor(Sharing, labels = c("Not Shared", "Shared")), fill = Persistent)) +
  geom_bar(position = "fill", width = 0.9) +
  geom_text(stat = "count", aes(label = ..count..), 
            position = position_fill(vjust = 0.5), color = "black", size = 4.5) +
  scale_fill_manual(values = c("Persistent" = "#457B9D", "Non-Persistent" = "#DADADA")) +
  labs(x = "Maternal sharing", y = "Proportion of vOTUs (%)", fill = "Phage persistence") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.15),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(add = c(0.05, 0.05))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()

pdf('11_PERSISTENCE/Plots/vOTU_DGR_by_sharing.pdf', width = 2.6, height = 3.6)
ggplot(Viral_metadata_infants_plot, aes(x = factor(Sharing, labels = c("Not Shared", "Shared")), fill = DGR)) +
  geom_bar(position = "fill", width = 0.9) +
  geom_text(stat = "count", aes(label = ..count..), 
            position = position_fill(vjust = 0.5), color = "black", size = 4.5) +
  scale_fill_manual(values = c("Yes" = "#457B9D", "No" = "#DADADA")) +
  labs(x = "Maternal sharing", y = "Proportion of vOTUs (%)", fill = "DGR presence") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.15),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(add = c(0.05, 0.05))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()

#***********************************************
# 4. Estimate number and proportion of DGR+ vOTUs
#***********************************************

DGR_votus_infants <- Viral_metadata_infants[Viral_metadata_infants$DGR == "Yes", "Virus_ID"]
DGR_votus_mothers <- Viral_metadata_mothers[Viral_metadata_mothers$DGR == "Yes", "Virus_ID"]

# Subset abundance table to DGR+ vOTUs and set to long format
# Infants
Ab_DGR_infants <- Abundance_table_infants[rownames(Abundance_table_infants) %in% DGR_votus_infants, ]
Ab_DGR_long_infants <- Ab_DGR_infants %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Virus_ID") %>%
  pivot_longer(
    cols = -Virus_ID,
    names_to = "NG_ID",
    values_to = "abundance"
  )

# Mothers
Ab_DGR_mothers <- Abundance_table_mothers[rownames(Abundance_table_mothers) %in% DGR_votus_mothers, ]
Ab_DGR_long_mothers <- Ab_DGR_mothers %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Virus_ID") %>%
  pivot_longer(
    cols = -Virus_ID,
    names_to = "NG_ID",
    values_to = "abundance"
  )

# Count number of DGR vOTUs present (>0) per sample
DGR_per_sample_infants <- Ab_DGR_long_infants %>%
  group_by(NG_ID) %>%
  summarise(Number_DGRs = sum(abundance > 0))

DGR_per_sample_mothers <- Ab_DGR_long_mothers %>%
  group_by(NG_ID) %>%
  summarise(Number_DGRs = sum(abundance > 0))

# Add to Sample_metadata_infants
Sample_metadata_infants <- Sample_metadata_infants %>%
  left_join(DGR_per_sample_infants, by = "NG_ID") %>%
  mutate(Number_DGRs = replace_na(Number_DGRs, 0)) 

Sample_metadata_mothers <- Sample_metadata_mothers %>%
  left_join(DGR_per_sample_mothers, by = "NG_ID") %>%
  mutate(Number_DGRs = replace_na(Number_DGRs, 0))  

# Estimate also the proportion of DGR+ vOTUs per sample
Sample_metadata_infants$Proportion_DGR_vOTUs <- 100*(Sample_metadata_infants$Number_DGRs / Sample_metadata_infants$richness)
Sample_metadata_mothers$Proportion_DGR_vOTUs <- 100*(Sample_metadata_mothers$Number_DGRs / Sample_metadata_mothers$richness)


#***********************************************
# 4. Association of proportion of DGR+ vOTUs with human phenotypes
#***********************************************

# Run associations with time (in mothers and infants) using a LMM
# Set W2 as reference and P12 as reference
Sample_metadata_infants$Timepoint_categorical <-
  relevel(as.factor(Sample_metadata_infants$Timepoint_categorical),
          ref = "W2")
Sample_metadata_mothers$Timepoint_categorical <-
  relevel(as.factor(Sample_metadata_mothers$Timepoint_categorical),
          ref = "P12")

# A) Infants
MM_time_DGR_presence_infant <- lmer(Proportion_DGR_vOTUs ~ exact_age_months_at_collection + DNA_concentration_ng_ul + read_depth + richness + (1|NEXT_ID),
                              REML = F, data = Sample_metadata_infants)
summary(MM_time_DGR_presence_infant)  # 3.862e-01  4.142e-02  2.643e+03   9.324  < 2e-16 ***

MM_time_DGR_presence_infant_cat <- lmer(Proportion_DGR_vOTUs ~ Timepoint_categorical + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID),
                                        REML = F, data = Sample_metadata_infants)
summary(MM_time_DGR_presence_infant_cat)  

# B) Mothers
MM_time_DGR_presence_mother <- lmer(Proportion_DGR_vOTUs  ~ exact_age_months_at_collection + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID),
                                REML = F, data = Sample_metadata_mothers)
summary(MM_time_DGR_presence_mother) #exact_age_months_at_collection -1.375e-02  4.637e-03  7.201e+02  -2.964 0.003134 **

MM_time_DGR_presence_mother_cat <- lmer(Proportion_DGR_vOTUs  ~ Timepoint_categorical + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID),
                                    REML = F, data = Sample_metadata_mothers)
summary(MM_time_DGR_presence_mother_cat) 


# Run associations with a subset of 10 infant phenotypes
phenotypes_to_test <- c("birth_deliverybirthcard_mode_binary", "birth_deliverybirthcard_place_delivery_simple",
                        "family_pets_any", "infant_health_eczema_diagnosis_relaxed", "infant_health_eczema_diagnosis_strict", 
                        "mother_birthcard_parity", "infant_health_food_allergy", "infant_growth_birth_weight_kg",  
                        "mother_birthcard_age_at_delivery", "infant_ffq_ever_never_breastfed" )

# Identify categorical and continuous variables
categorical_vars <- names(Sample_metadata_infants)[
  sapply(Sample_metadata_infants, function(x)
    is.character(x) || (is.numeric(x) && length(unique(na.omit(x))) < 10))
]

continuous_vars <- setdiff(names(Sample_metadata_infants), categorical_vars)

# Set as factors categorical variables
Sample_metadata_infants <- Sample_metadata_infants %>%
  mutate(across(all_of(categorical_vars), as.factor))

# Reorder birth_deliverybirthcard_place_delivery_simple
Sample_metadata_infants$birth_deliverybirthcard_place_delivery_simple <- 
  factor(Sample_metadata_infants$birth_deliverybirthcard_place_delivery_simple,
         levels = c("hospital","home"))

# Set parity as a binary phenotype (no if 0, yes if > 0 )
Sample_metadata_infants <- Sample_metadata_infants %>%
  mutate(
    parity_num = as.numeric(as.character(mother_birthcard_parity)),
    mother_birthcard_parity = if_else(parity_num > 0, "yes", "no")
  )

# Test association with the proportion of DGR+ vOTUs
lmm_results_infant_DGRs <- test_phenotype_associations_lmm(
  df = Sample_metadata_infants,          
  outcome_var = "Proportion_DGR_vOTUs", 
  phenotypes = phenotypes_to_test,
  subject_col = "NEXT_ID",
  covariates = c("clean_reads_FQ_2", "DNA_concentration_ng_ul", "Timepoint_categorical")
)

# Test association after adjusting for feeding mode and delivery mode
lmm_results_infant_DGR_adj  <- test_phenotype_associations_lmm(
  df = Sample_metadata_infants,          
  outcome_var = "Proportion_DGR_vOTUs", 
  phenotypes = setdiff(phenotypes_to_test, c("infant_ffq_ever_never_breastfed",
                                             "birth_deliverybirthcard_mode_binary")),
  subject_col = "NEXT_ID",
  covariates = c("clean_reads_FQ_2", "DNA_concentration_ng_ul", "Timepoint_categorical",
                 "infant_ffq_ever_never_breastfed", "birth_deliverybirthcard_mode_binary")
)

# Order results by FDR
lmm_results_infant_DGRs <- lmm_results_infant_DGRs %>%
  arrange(FDR)
lmm_results_infant_DGR_adj <- lmm_results_infant_DGR_adj %>%
  arrange(FDR)

################
# --- Plots ----
################

# A) --- Time ----

Sample_metadata_infants$Timepoint_categorical <- factor(Sample_metadata_infants$Timepoint_categorical, 
                                                levels=c("W2", "M1", "M2", "M3","M6","M9", "M12"))
Sample_metadata_mothers$Timepoint_categorical <- factor(Sample_metadata_mothers$Timepoint_categorical, 
                                                        levels=c("P12", "P28", "B", "M3"))


pdf('5_DGR_ANALYSIS/Plots/DGR_Presence_Infants.pdf', width=3.6, height=3.1)
ggplot(Sample_metadata_infants, aes(x = Timepoint_categorical, y = Proportion_DGR_vOTUs, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Timepoint', y = 'Proportion DGR+ vOTUs') +
  scale_fill_manual(values = c("Infant" = "#66A6AD")) +
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

pdf('5_DGR_ANALYSIS/Plots/DGR_Presence_Mothers.pdf', width=3.6, height=3.1)
ggplot(Sample_metadata_mothers, aes(x = Timepoint_categorical, y = Proportion_DGR_vOTUs, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Timepoint', y = 'Proportion DGR+ vOTUs') +
  scale_fill_manual(values = c("Mother" = "#8E7CA6")) +
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


# B) --- Delivery mode and allergy ----

pdf('5_DGR_ANALYSIS/Plots/DGR_Presence_Infants_Delivery.pdf', width=2, height=3.2)
ggplot(Sample_metadata_infants %>% filter(!is.na(delivery_mode)),
       aes(x = delivery_mode, y = Proportion_DGR_vOTUs, fill = delivery_mode)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Delivery', y = 'Proportion DGR+ vOTUs') +
  scale_fill_manual(values = c("VG" = "#B6D7A8", "CS" = "#A2C2E5")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()

pdf('5_DGR_ANALYSIS/Plots/DGR_Presence_Infants_Delivery_time.pdf', width = 5.5, height = 3.9)
ggplot(Sample_metadata_infants %>% filter(!is.na(delivery_mode)),  aes(x = Timepoint_categorical,
           y = Proportion_DGR_vOTUs, fill = delivery_mode)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  geom_boxplot(width = 0.7, color = "black", outlier.color = NA, size = 0.9,
               position = position_dodge(width = 0.8)) +
  labs(x = 'Timepoint', y = 'Proportion DGR+ vOTUs') +
  scale_fill_manual(values = c("VG" = "#B6D7A8", "CS" = "#A2C2E5")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()


pdf('5_DGR_ANALYSIS/Plots/DGR_Presence_Infants_Allergy.pdf', width=2, height=3.2)
ggplot(Sample_metadata_infants %>% filter(!is.na(infant_health_food_allergy)),
       aes(x = infant_health_food_allergy, y = Proportion_DGR_vOTUs, fill = infant_health_food_allergy)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Food Allergy', y = 'Proportion DGR+ vOTUs') +
  scale_fill_manual(values = c("no" = "#4C72B0", "yes" = "#DD8452")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.5),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()


pdf('5_DGR_ANALYSIS/Plots/DGR_Presence_Infants_Allergy_time.pdf', width = 5.5, height = 3.9)
ggplot(Sample_metadata_infants %>% filter(!is.na(infant_health_food_allergy)),  aes(x = Timepoint_categorical,
                                                                       y = Proportion_DGR_vOTUs, fill = infant_health_food_allergy)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, 
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  geom_boxplot(width = 0.7, color = "black", outlier.color = NA, size = 0.9,
               position = position_dodge(width = 0.8)) +
  labs(x = 'Timepoint', y = 'Proportion DGR+ vOTUs') +
  scale_fill_manual(values = c("no" = "#4C72B0", "yes" = "#DD8452")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 18),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()


#***********************************************
# 5. Write results
#***********************************************
write.table(DGR_target_results,"5_DGR_ANALYSIS/Results/DGR_target_results.txt",
sep = "\t", row.names = F, quote = FALSE) 
write.table(Annotation_DGR_targets,"5_DGR_ANALYSIS/Results/Annotation_DGR_targets.txt",
sep = "\t", row.names = F, quote = FALSE) 
write.table(lmm_results_infant_DGRs,"5_DGR_ANALYSIS/Results/Proportion_DGR_vOTUs_association_INFANTS.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(lmm_results_infant_DGR_adj,"5_DGR_ANALYSIS/Results/Proportion_DGR_vOTUs_association_adj_INFANTS.txt",
            sep = "\t", row.names = F, quote = FALSE) 

