################################################################################
##### LL-NEXT:  DGR analysis -  Activity in MOTHERS (II)
### Author(s): Asier Fernández-Pato
### Last updated: 18th December, 2025
################################################################################

#****************
# Load libraries
#****************
library(ggplot2)
library(ggalluvial)
library(dplyr)
library(tidyr)
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
dgr_activity_table <- read.delim("5_DGR_ANALYSIS/MOTHERS/DGR_activity_table_MOTHERS.txt")

#***********************************************
# 1. Estimate number of active DGRs and mutation bias in codons
#***********************************************

# Estimate the number of VRs (and DGRs) that are potentially active (based on per-site nucleotide diversity) in the VR
length(which(dgr_activity_table$P_value_VR_var_adj < 0.05)) # 1801 VRs potentially active
active_VRs_per_site_var <- dgr_activity_table[!is.na(dgr_activity_table$P_value_VR_var_adj) &
                                                dgr_activity_table$P_value_VR_var_adj < 0.05,]
active_DGRs_per_site_var <- active_VRs_per_site_var %>% # n=1609 DGRs
  group_by(DGR_code) %>%
  filter(VR_ID %in% unique(VR_ID)[1]) %>%
  ungroup() 
nrow(active_DGRs_per_site_var)

# Indicate which of the active DGRs also show temporal VR diversification of NT and AA
table(active_VRs_per_site_var$long_NT_diversification) # 1509
table(active_DGRs_per_site_var$long_NT_diversification) # 1353 
table(active_VRs_per_site_var$long_AA_diversification) # 1483 
table(active_DGRs_per_site_var$long_AA_diversification) # 1333

# Indicate which of the non-active DGRs (Wilcox test) show temporal VR diversification of NT and AA
non_active_VRs_per_site_var <- dgr_activity_table[
  !is.na(dgr_activity_table$P_value_VR_var_adj) &
    dgr_activity_table$P_value_VR_var_adj >= 0.05, ]

non_active_DGRs_per_site_var <- non_active_VRs_per_site_var %>%
  group_by(DGR_code) %>%
  filter(VR_ID %in% unique(VR_ID)[1]) %>%
  ungroup()

table(non_active_VRs_per_site_var$long_NT_diversification) #139
table(non_active_DGRs_per_site_var$long_NT_diversification) #118
table(non_active_VRs_per_site_var$long_AA_diversification) #136
table(non_active_DGRs_per_site_var$long_AA_diversification) #115

# Estimate % of codons with temporal NT changes that resulted in AA substitutions (88.1%)
fraction_codons_with_aa_substitutions <- 100 * (
  sum(dgr_activity_table$codons_with_AA_changes, na.rm = TRUE) /
    sum(dgr_activity_table$codons_with_NT_changes, na.rm = TRUE)
)

# Statistical test: GLMM (accounting for multiple vOTUs and infants with repeated VRs)
# Count of mutations in pos1+2 vs pos3
dgr_activity_table$mut_pos1_2 <- dgr_activity_table$codon_pos1_changes + dgr_activity_table$codon_pos2_changes
dgr_activity_table$mut_pos3 <- dgr_activity_table$codon_pos3_changes

# Fit GLMM
GLMM_base_mutation_bias <- glmer(cbind(mut_pos1_2, mut_pos3) ~ 1 + (1 | vOTU_ID) + (1 | NEXT_ID),
  data = dgr_activity_table, family = binomial)

GLMM_base_mutation_bias_results <- summary(GLMM_base_mutation_bias) # p = 2.019502e-308
GLMM_base_mutation_bias_results$coefficients["(Intercept)", "Pr(>|z|)"]

dgr_activity_table$mut_pos1_2 <- NULL
dgr_activity_table$mut_pos3 <- NULL


############
# PLOTS
############

# A) Plot showing the distribution of the number of mothers in which of the tested DGRs is persistent (n=229)
vOTU_persistence <- read.delim("5_DGR_ANALYSIS/MOTHERS/vOTU_persistence_mother_per_vOTU.txt")

vOTU_persistence <- vOTU_persistence %>%
  left_join(
    Viral_metadata_mothers %>% 
      select(Virus_ID, Lifestyle, Bacterial_genus_host),
    by = c("vOTU" = "Virus_ID")
  )

pdf('5_DGR_ANALYSIS/Plots/Persistent_vOTU_distribution_mothers.pdf', width = 3.9, height = 3.9)
ggplot(vOTU_persistence, aes(x = n_persistent_mothers)) +
  stat_ecdf(geom = "step", color = "#1F4E4E", linewidth = 1.2) +
  geom_vline(xintercept = median(vOTU_persistence$n_persistent_mothers),
             linetype = "dashed", color = "#4B4FC5", linewidth = 1) +
  labs(x = "Number of mothers (persistent)", y = "Cumulative fraction of vOTUs") +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.02))) + 
  theme_bw(base_size = 16) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 16)
  )
dev.off()

# B) Sankey plot showing DGR activity
sankey_data <- data.frame(
  total = c("Total DGRs", "Total DGRs", "Total DGRs", "Total DGRs"),
  hypervariable = c("Hypervariable", "Hypervariable",
                    "Not hypervariable", "Not hypervariable"),
  temporal = c("Temporal diversification", "No temporal diversification",
               "Temporal diversification", "No temporal diversification"),
  count = c(1353, 256, 118, 1613) 
)

#HP and diversi 1353 
#HP and not diversity 256
#Not HP and divesi 118
#Not hypervariable not diversification 1613

# Convert to factors and reorder levels
sankey_data <- sankey_data %>%
  mutate(
    total = factor(total, levels = c("Total DGRs")),
    hypervariable = factor(hypervariable, levels = c("Hypervariable", "Not hypervariable")),
    temporal = factor(temporal, levels = c("No temporal diversification", "Temporal diversification"))
  )

# Generate the plot
pdf('5_DGR_ANALYSIS/Plots/DGR_activity_Sankey_MOTHERS.pdf', width = 6.5, height = 3.9)
ggplot(sankey_data, aes(axis1 = total, axis2 = hypervariable, axis3 = temporal, y = count)) +
  geom_alluvium(aes(fill = temporal), width = 1/12) +
  geom_stratum(width = 1/9, fill = "#E8E8E8", color = "black") +
  geom_text(stat = "stratum", aes(label = after_stat(stratum)), size = 4) +
  scale_fill_manual(values = c(
    "Temporal diversification" = "#4B4FC5",
    "No temporal diversification" = "#D3D3D3"
  )) +
  labs(y = "Number of DGRs", x = "") +
  theme_minimal(base_size = 14) +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )
dev.off()
