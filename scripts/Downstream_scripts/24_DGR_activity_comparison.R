################################################################################
##### LL-NEXT:  DGR analysis - Comparison in mothers vs infants and shared vOTUs
### Author(s): Asier Fernández-Pato
### Last updated: 18th December, 2025
################################################################################

#****************
# Load libraries
#****************
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(scales)
library(lme4)
library(lmerTest)


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

# Read DGR activity results
dgr_activity_table_infants <- read.delim("5_DGR_ANALYSIS/INFANTS/DGR_activity_table.txt")
dgr_activity_table_mothers <- read.delim("5_DGR_ANALYSIS/MOTHERS/DGR_activity_table_MOTHERS.txt")


#*#***********************************************
# 1. Comparison of DGR activity (NT changes) in mothers and infants
#***********************************************

# Add binary activity variable (here NT changes will be used)
dgr_activity_table_mothers$long_AA_div_binary <- ifelse(
  dgr_activity_table_mothers$long_AA_diversification == "yes", 1, 0
)

dgr_activity_table_infants$long_AA_div_binary <- ifelse(
  dgr_activity_table_infants$long_AA_diversification == "yes", 1, 0
)

dgr_activity_table_mothers$long_NT_div_binary <- ifelse(
  dgr_activity_table_mothers$long_NT_diversification == "yes", 1, 0
)

dgr_activity_table_infants$long_NT_div_binary <- ifelse(
  dgr_activity_table_infants$long_NT_diversification == "yes", 1, 0
)


# Add mother/infant information
dgr_activity_table_mothers$group <- "mother"
dgr_activity_table_infants$group <- "infant"

# Generate a combined table with DGR activity results in mothers and infants
df_mothers <- dgr_activity_table_mothers %>%
  select(DGR_code,VR_ID, vOTU_ID, NEXT_ID, long_NT_div_binary, long_AA_div_binary, group,
         pregnancy_AA_change,post_pregnancy_AA_change,
         DGR_target_ann_description, DGR_target_ann_category)

df_infants <- dgr_activity_table_infants %>%
  select(DGR_code,VR_ID, vOTU_ID, NEXT_ID, long_NT_div_binary, num_var_VR_positions, long_AA_div_binary, group,
  first_month_AA_change, early_AA_change, late_AA_change, early_to_late_AA_change,
  DGR_target_ann_description, DGR_target_ann_category)

dgr_activity_df <- bind_rows(df_mothers, df_infants)

# Run GLMM to test whether whether the probability of a DGR being active differs between infants and mothers
# We account for repeated measurements at the vOTU and infant level
GLMM_DGR_activity_infant_mother <- glmer(long_NT_div_binary ~ group + (1 | vOTU_ID) + (1 | NEXT_ID),
  data = dgr_activity_df, family = binomial)
summary_GLMM_DGR_activity_infant_mother <- summary(GLMM_DGR_activity_infant_mother)
summary_GLMM_DGR_activity_infant_mother$coefficients #groupmother -1.079944  0.1260200 -8.569623 1.038247e-17

###########
# PLOTS
###########

# Get the table of active DGRs based on nucleotide variability in mothers and infants
active_VRs_mothers <- dgr_activity_table_mothers[dgr_activity_table_mothers$long_NT_diversification == "yes",]

active_DGRs_mothers <- active_VRs_mothers %>% # n=1472 DGRs
  group_by(DGR_code) %>%
  filter(VR_ID %in% unique(VR_ID)[1]) %>%
  ungroup() 
nrow(active_DGRs_mothers)

active_VRs_infants <- dgr_activity_table_infants[dgr_activity_table_infants$long_NT_diversification == "yes",]

active_DGRs_infants <- active_VRs_infants %>% # n=415 DGRs
  group_by(DGR_code) %>%
  filter(VR_ID %in% unique(VR_ID)[1]) %>%
  ungroup() 
nrow(active_DGRs_infants)

# Knowing the total number of DGRs explored and the active ones (NT changes), build a table
n_total_DGRs_mothers <- 3340
n_total_DGRs_infants <- 628
n_active_DGRs_mothers <- nrow(active_DGRs_mothers)
n_active_DGRs_infants <- nrow(active_DGRs_infants)

df_prop <- data.frame(
  Group = c("Infants", "Mothers"),
  Active = c(n_active_DGRs_infants, n_active_DGRs_mothers),
  Total = c(n_total_DGRs_infants, n_total_DGRs_mothers)
)

df_prop$Inactive <- df_prop$Total - df_prop$Active

# Convert to long format
df_long <- tidyr::pivot_longer(
  df_prop,
  cols = c("Active", "Inactive"),
  names_to = "Status",
  values_to = "Count"
)

# Estimate proportions
df_long$Proportion <- df_long$Count / rep(df_prop$Total, each = 2)
df_long$Status <- factor(df_long$Status, levels = c("Inactive", "Active"))

fill_colors <- c("Active" = "#457B9D", "Inactive" = "#DADADA")

pdf("5_DGR_ANALYSIS/Plots/Active_DGRs_infants_vs_mothers.pdf", width = 2.3, height = 3.6)
ggplot(df_long, aes(x = Group, y = Proportion, fill = Status)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_text(data = df_long[df_long$Status == "Active", ],
    aes(label = paste0(round(Proportion * 100, 1), "%\n(",Count, "/", Total, ")")),
    position = position_stack(vjust = 0.5), color = "black", size = 4.5) +
  scale_fill_manual(values = fill_colors) +
  scale_y_continuous(labels = scales::percent,limits = c(0, 1.15), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(add = c(0.05, 0.05))) +
  labs(x = "", y = "Active DGRs (%)", fill = "") +
  theme_bw(base_size = 15) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
    )
dev.off()


#***********************************************
# 2. DGR activity (AA changes) in shared vOTUs
#***********************************************

# Read the list of shared vOTUs
shared_vOTUs <- read.delim("10_STRAIN_TRANSMISSION/shared_vOTUs.txt")
shared_vOTUs <- shared_vOTUs$x

# Select DGR+ vOTUs among shared ones (n=99)
DGR_vOTUs <- Viral_metadata[Viral_metadata$DGR == "Yes", "Virus_ID", drop = TRUE]
DGR_shared_vOTUs <- shared_vOTUs[shared_vOTUs %in% DGR_vOTUs]

# Get number of DGRs among shared vOTUs
df_infants_DGRs <- df_infants %>% 
  group_by(DGR_code) %>%
  filter(VR_ID %in% unique(VR_ID)[1]) %>%
  ungroup() 

df_infants_DGRs$shared <- ifelse(df_infants_DGRs$vOTU_ID %in% DGR_shared_vOTUs, "yes", "no")
table(df_infants_DGRs$shared) # 201 DGRs in shared vOTUs

# GLMM: Test if shared vOTUs are likely more active in the infant gut than non-shared ones
# We account for repeated measurements at the vOTU and infant level
df_infants$shared <- ifelse(df_infants$vOTU_ID %in% DGR_shared_vOTUs, "yes", "no")
GLMM_DGR_activity_shared_infant <- glmer(long_AA_div_binary ~ shared + (1 | vOTU_ID) + (1 | NEXT_ID),
                                         data = df_infants, family = binomial)
summary_GLMM_DGR_activity_shared_infant <- summary(GLMM_DGR_activity_shared_infant)
summary_GLMM_DGR_activity_shared_infant$coefficients

# Estimate proportion of DGRs from shared and non-shared vOTUs that show AA changes
shared_DGRs <- df_infants_DGRs[df_infants_DGRs$shared == "yes", ]
n_shared <- nrow(shared_DGRs)
n_AA_changes_shared <- sum(shared_DGRs$long_AA_div_binary == "1", na.rm = TRUE)
prop_AA_changes_shared <- n_AA_changes_shared / n_shared #67.6

non_shared_DGRs <- df_infants_DGRs[df_infants_DGRs$shared == "no", ]
n_non_shared <- nrow(non_shared_DGRs)
n_AA_changes_non_shared <- sum(non_shared_DGRs$long_AA_div_binary == "1", na.rm = TRUE)
prop_AA_changes_non_shared <- n_AA_changes_non_shared / n_non_shared #64.6

# Next, subset results for shared vOTUs and merge infant and mothers
df_infants_shared <- df_infants %>%
  filter(vOTU_ID %in% shared_vOTUs)
df_mothers_shared <- df_mothers %>%
  filter(vOTU_ID %in% shared_vOTUs)

df_infants_shared$Type <- "Infant"
df_mothers_shared$Type <- "Mother"
df_mothers_shared$shared <- "yes"
  
df_infants_shared <- df_infants_shared[,c("DGR_code","VR_ID", "vOTU_ID", "NEXT_ID",
                                "Type", "shared", "long_AA_div_binary", 
                                "DGR_target_ann_description", "DGR_target_ann_category")]
df_mothers_shared <- df_mothers_shared[,c("DGR_code","VR_ID", "vOTU_ID", "NEXT_ID",
                                          "Type", "shared", "long_AA_div_binary", 
                                          "DGR_target_ann_description", "DGR_target_ann_category")]

df_shared_merged <- bind_rows(df_infants_shared, df_mothers_shared)

# GLMM: Test if the probability of DGR activity (AA diversification) differs between mothers and infants among shared vOTUs
df_shared_merged$Type <- relevel(factor(df_shared_merged$Type), ref = "Mother")
GLMM_DGR_activity_shared_mother_vs_infant <- glmer(long_AA_div_binary ~ Type + (1 | vOTU_ID) + (1 | NEXT_ID),
                                         data = df_shared_merged, family = binomial)
summary_GLMM_DGR_activity_shared_mother_vs_infant <- summary(GLMM_DGR_activity_shared_mother_vs_infant )
summary_GLMM_DGR_activity_shared_mother_vs_infant$coefficients #TypeInfant  0.8534795  0.2103729 4.056984 4.971049e-05

# Format both df_infants and df_shared_merged for supplementary tables
df_infants_sup <- df_infants[,c("DGR_code","VR_ID", "vOTU_ID", "NEXT_ID",
                                "group", "shared", "long_AA_div_binary", 
                                "DGR_target_ann_description", "DGR_target_ann_category")]

df_shared_merged_sup <- df_shared_merged[,c("DGR_code","VR_ID", "vOTU_ID", "NEXT_ID",
                                            "Type", "shared","long_AA_div_binary",
                                            "DGR_target_ann_description", "DGR_target_ann_category")]

lookup_family <- Sample_metadata %>%
  distinct(NEXT_ID, FAMILY)

df_infants_sup <- df_infants_sup %>%
  left_join(lookup_family, by = "NEXT_ID") %>%
  select(-NEXT_ID)

df_shared_merged_sup <- df_shared_merged_sup %>%
  left_join(lookup_family, by = "NEXT_ID") %>%
  select(-NEXT_ID)

#***********************************************
# 3. DGR activity patterns in infants/mothers of shared vOTUs
#***********************************************

# Add FAMILY IDs to infants and mothers
sample_unique <- Sample_metadata %>%
  distinct(NEXT_ID, FAMILY)

df_infants_shared <- df_infants_shared %>%
  left_join(sample_unique, by = "NEXT_ID") %>%
  dplyr::rename(FAM_ID_infant = FAMILY)

df_mothers_shared <- df_mothers_shared %>%
  left_join(sample_unique, by = "NEXT_ID") %>%
  dplyr::rename(FAM_ID_mother = FAMILY)

# Create unique IDs for each VR-FAMILY 
df_infants_shared$Unique_ID <- paste(df_infants_shared$FAM_ID_infant, df_infants_shared$VR_ID, sep = "_")
df_mothers_shared$Unique_ID <- paste(df_mothers_shared$FAM_ID_mother, df_mothers_shared$VR_ID, sep = "_")

# Generate one combined df with infant and maternal pairs in the same row
# Keep only those Unique_IDs present in df_infants_shared
# Note that in 3 cases DGRs were shared to both siblings of a family (but in all not persistent in the mother - no problem for merging just based on FAM)
paired_DGR_activity_shared <- df_infants_shared %>%
  inner_join(df_mothers_shared,
             by = "Unique_ID",
             suffix = c("_infant", "_mother"))

# For simplicity, remove comparisons (n=3) including the second pregnancy (if 1st pregnancy also present), so each pair is considered once
paired_DGR_activity_shared <- paired_DGR_activity_shared %>%
  filter(!str_ends(NEXT_ID_mother, "_2"))
paired_DGR_activity_shared$Pair_ID <- paste(paired_DGR_activity_shared$NEXT_ID_mother,
                                            paired_DGR_activity_shared$NEXT_ID_infant, sep = "_")

# Count number of VRs, DGRs and pairs assessed
n_VRs_pairs_shared <- nrow(paired_DGR_activity_shared) #102
DGR_pairs_shared <- paired_DGR_activity_shared %>% 
  group_by(DGR_code_infant) %>%
  filter(VR_ID_infant %in% unique(VR_ID_infant)[1]) %>%
  ungroup() 
n_DGRs_pairs_shared <- nrow(DGR_pairs_shared) #93
n_pairs_shared <- length(unique(paired_DGR_activity_shared$Pair_ID)) #45

# Summarize diversification patterns per mother-infant pair (for each VR)
summary_DGR_activity_shared_pairs <- paired_DGR_activity_shared %>%
  mutate(
    pattern = case_when(
      long_AA_div_binary_mother == 1 & long_AA_div_binary_infant == 1 ~ "Both mother and infant",
      long_AA_div_binary_mother == 1 & long_AA_div_binary_infant == 0 ~ "Mother only",
      long_AA_div_binary_mother == 0 & long_AA_div_binary_infant == 1 ~ "Infant only",
      long_AA_div_binary_mother == 0 & long_AA_div_binary_infant == 0 ~ "Neither"
    )
  ) %>%
  count(pattern) %>%
  mutate(percent = n / sum(n) * 100)

summary_DGR_activity_shared_pairs

#****************
# Write results
#****************
write.table(df_infants_sup,"5_DGR_ANALYSIS/Results/Shared_vOTUs_DGR_activity_infants.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(df_shared_merged_sup,"5_DGR_ANALYSIS/Results/Shared_vOTUs_DGR_activity_infants_mothers.txt",
            sep = "\t", row.names = F, quote = FALSE) 

