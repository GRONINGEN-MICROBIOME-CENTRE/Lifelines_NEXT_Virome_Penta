################################################################################
##### LL-NEXT:  DGR analysis -  Activity in INFANTS (II)
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
dgr_activity_table <- read.delim("5_DGR_ANALYSIS/INFANTS/DGR_activity_table.txt")

#***********************************************
# 1. Estimate number of active DGRs and mutation bias in codons
#***********************************************

# Estimate the number of VRs (and DGRs) that are potentially active (based on per-site nucleotide diversity) in the VR
length(which(dgr_activity_table$P_value_VR_var_adj < 0.05)) # 475 VRs potentially active
active_VRs_per_site_var <- dgr_activity_table[!is.na(dgr_activity_table$P_value_VR_var_adj) &
                                                dgr_activity_table$P_value_VR_var_adj < 0.05,]
active_DGRs_per_site_var <- active_VRs_per_site_var %>% # n=439 DGRs
  group_by(DGR_code) %>%
  filter(VR_ID %in% unique(VR_ID)[1]) %>%
  ungroup() 
nrow(active_DGRs_per_site_var)

# Indicate which of the active DGRs (Wilcox test) also show temporal VR diversification of NT and AA
table(active_VRs_per_site_var$long_NT_diversification) # 425
table(active_DGRs_per_site_var$long_NT_diversification) # 392
table(active_VRs_per_site_var$long_AA_diversification) # 425
table(active_DGRs_per_site_var$long_AA_diversification) # 392

# Indicate which of the non-active DGRs (Wilcox test) show temporal VR diversification of NT and AA
non_active_VRs_per_site_var <- dgr_activity_table[
  !is.na(dgr_activity_table$P_value_VR_var_adj) &
    dgr_activity_table$P_value_VR_var_adj >= 0.05, ]

non_active_DGRs_per_site_var <- non_active_VRs_per_site_var %>%
  group_by(DGR_code) %>%
  filter(VR_ID %in% unique(VR_ID)[1]) %>%
  ungroup()

table(non_active_VRs_per_site_var$long_NT_diversification) #28
table(non_active_DGRs_per_site_var$long_NT_diversification) #26
table(non_active_VRs_per_site_var$long_AA_diversification) #28
table(non_active_DGRs_per_site_var$long_AA_diversification) #26

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

GLMM_base_mutation_bias_results <- summary(GLMM_base_mutation_bias) # p = 1.293621e-89
GLMM_base_mutation_bias_results $coefficients["(Intercept)", "Pr(>|z|)"]

dgr_activity_table$mut_pos1_2 <- NULL
dgr_activity_table$mut_pos3 <- NULL

#***********************************************
# 2. Association with infant phenotypes
#***********************************************

# Use stricter definition of activity (longitudinal diversification of nucleotides == yes)
# As VRs are assessed sometimes across same vOTUs and infants, a GLMM should be used

# First, convert DGR activity to binary 0/1
dgr_activity_table <- dgr_activity_table %>%
  mutate(long_NT_div_binary = ifelse(long_NT_diversification == "yes", 1, 0))

# Add viral features
dgr_activity_table <- dgr_activity_table %>%
  left_join(
    Viral_metadata_infants %>%
      select(Virus_ID, Lifestyle, Bacterial_genus_host),
    by = c("vOTU_ID" = "Virus_ID")
  ) %>%
  # Create a binary host variable
  mutate(Host_Bacteroides = ifelse(Bacterial_genus_host == "Bacteroides", "Bacteroides", "Non-Bacteroides"))

# Add human phenotypes
phenotypes_to_test <- c("birth_deliverybirthcard_mode_binary", "infant_ffq_ever_never_breastfed")

dgr_activity_table <- dgr_activity_table %>%
  left_join(
    Sample_metadata_infants %>%
      select(NEXT_ID, all_of(phenotypes_to_test)),
    by = "NEXT_ID"
  ) %>%
  distinct()

# Convert character phenotypes to factors
phenotypes_to_test <- c("Lifestyle", phenotypes_to_test)

for (pheno in phenotypes_to_test) {
  if (is.character(dgr_activity_table[[pheno]])) {
    dgr_activity_table[[pheno]] <- factor(dgr_activity_table[[pheno]])
  }
}

# Initialize results table
GLM_DGR_activity_phenos <- data.frame(
  Variable = character(),
  Level = character(),
  Estimate = numeric(),
  Std_Error = numeric(),
  z_value = numeric(),
  p_value = numeric(),
  stringsAsFactors = FALSE
)

for (i in phenotypes_to_test) {
  
  formula_glmm <- as.formula(paste0("long_NT_div_binary ~ ", i, " + (1 | vOTU_ID) + (1 | NEXT_ID)"))
  fit <- glmer(formula_glmm, family = binomial, data = dgr_activity_table,
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5)))
  sum_fit <- summary(fit)
  
  # Extract coefficients for all levels (excluding intercept)
  coef_table <- sum_fit$coefficients[-1, , drop = FALSE]
  
  for (lvl in rownames(coef_table)) {
    GLM_DGR_activity_phenos <- rbind(
      GLM_DGR_activity_phenos,
      data.frame(
        Variable = i,
        Level = lvl,
        Estimate = coef_table[lvl, "Estimate"],
        Std_Error = coef_table[lvl, "Std. Error"],
        z_value = coef_table[lvl, "z value"],
        p_value = coef_table[lvl, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      )
    )
  }
}

# Adjust p-values for multiple testing (FDR)
GLM_DGR_activity_phenos$FDR <- p.adjust(GLM_DGR_activity_phenos$p_value, method = "fdr")

#***********************************************
# 3. PLOTS
#***********************************************

# A) Plot showing the distribution of the number of infants in which of the tested DGRs is persistent (n=229)
vOTU_persistence <- read.delim("5_DGR_ANALYSIS/INFANTS/vOTU_persistence.txt")

vOTU_persistence <- vOTU_persistence %>%
  left_join(
    Viral_metadata_infants %>% 
      select(Virus_ID, Lifestyle, Bacterial_genus_host),
    by = c("vOTU" = "Virus_ID")
  )

pdf('5_DGR_ANALYSIS/Plots/Persistent_vOTU_distribution_infants.pdf', width = 4.1, height = 3.9)
ggplot(vOTU_persistence, aes(x = n_persistent_infants)) +
  geom_histogram(binwidth = 1, fill = "#D3D3D3", color = "black", alpha = 0.8) +
  geom_vline(xintercept = median(vOTU_persistence$n_persistent_infants), 
             color = "#4B4FC5", size = 1.2, linetype = "dashed") +
  labs(x = "Number of infants (persistent)", y = "Number of vOTUs") +
  annotate("text",
    x = median(vOTU_persistence$n_persistent_infants), 
    y = max(table(vOTU_persistence$n_persistent_infants)) * 0.9,
    label = paste0("Median = ", median(vOTU_persistence$n_persistent_infants)),
    hjust = -0.1, color = "#4B4FC5", size = 5
  ) +
  coord_cartesian(ylim = c(0, max(table(vOTU_persistence$n_persistent_infants)) * 1.1),
                  xlim = c(0, max(vOTU_persistence$n_persistent_infants) + 1)) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
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
  count = c(392, 47, 26, 163)  # 50 + 113 = 163 for Not hypervariable, no diversification
)

# Convert to factors and reorder levels
sankey_data <- sankey_data %>%
  mutate(
    total = factor(total, levels = c("Total DGRs")),
    hypervariable = factor(hypervariable, levels = c("Hypervariable", "Not hypervariable")),
    temporal = factor(temporal, levels = c("No temporal diversification", "Temporal diversification"))
  )

# Generate the plot
pdf('5_DGR_ANALYSIS/Plots/DGR_activity_Sankey_INFANTS.pdf', width = 8, height = 3.9)
ggplot(sankey_data, aes(axis1 = total, axis2 = hypervariable, axis3 = temporal, y = count)) +
  geom_alluvium(aes(fill = temporal), width = 1/12) +
  geom_stratum(width = 0.3, fill = "#E8E8E8", color = "black") +
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

# C) Plot codon position bias and non-synonimous changes
# Summarize counts per codon position
codon_counts <- dgr_activity_table %>%
  summarise(
    pos1 = sum(codon_pos1_changes, na.rm = TRUE),
    pos2 = sum(codon_pos2_changes, na.rm = TRUE),
    pos3 = sum(codon_pos3_changes, na.rm = TRUE)
  ) %>%
  pivot_longer(everything(),
               names_to = "Codon_position",
               values_to = "Mutations") %>%
  mutate(Proportion = Mutations / sum(Mutations))

# Choose colors
colors <- c("pos1" = "#1F4E4E","pos2" = "#76B39D", "pos3" = "#A0D6C7")

pdf("5_DGR_ANALYSIS/Plots/Barplot_codon_position_bias.pdf", width = 2.5, height = 3.5)
ggplot(codon_counts, aes(x = Codon_position, y = Proportion, fill = Codon_position)) +
  geom_col(width = 0.8) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(labels = scales::percent, expand = expansion(mult = c(0.05, 0.25))) +
  labs(x = "Codon position", y = "Proportion of temporal mutations") +
  theme_bw(base_size = 14) +
  theme(
    panel.grid = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "none"
  )
dev.off()

aa_fraction <- fraction_codons_with_aa_substitutions
fractions <- data.frame(
  Type = c("AA substitutions", "Silent/NT only"),
  Percent = c(aa_fraction, 100 - aa_fraction)
)

aa_colors <- c("AA substitutions" = "#80CBC4", 
               "Silent/NT only"   = "#D7CCC8")

pdf("5_DGR_ANALYSIS/Plots/AA_substitution_fraction_pie.pdf", width = 2, height = 2)
ggplot(fractions, aes(x = "", y = Percent, fill = Type)) +
  geom_col(width = 1, color = NA, alpha = 0.85) + 
  coord_polar(theta = "y") +
  geom_text(aes(label = paste0(round(Percent, 1), "%")),
            position = position_stack(vjust = 0.5),
            size = 4.5, color = "black") +
  scale_fill_manual(values = aa_colors) +
  theme_void() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 13, face = "plain"),
    legend.position = "none"
  )
dev.off()


#****************
# Write results
#****************
write.table(GLM_DGR_activity_phenos,"5_DGR_ANALYSIS/Results/GLM_DGR_activity_phenos.txt", sep = "\t", row.names = F, quote = FALSE) 


