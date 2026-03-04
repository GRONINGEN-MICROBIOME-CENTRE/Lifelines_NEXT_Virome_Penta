################################################################################
##### LL-NEXT: Estimation of vOTU persistence
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


#****************
# Define functions
#****************
# Define function to estimate the abundance explained by persistent vOTUs per sample
compute_persistent_abundance <- function(sample_id, abundance_table, metadata, persistent_list) {
  next_id <- metadata$NEXT_ID[metadata$NG_ID == sample_id][1] # [1] added because Sample_metadata_mothers can have 2 entries for the same sample (if multiple pregnancies)
  
  if (length(next_id) == 0 || is.na(next_id) || !(next_id %in% names(persistent_list))) {
    return(NA)
  }
  
  persistent_votus <- persistent_list[[next_id]]
  persistent_votus <- intersect(persistent_votus, rownames(abundance_table))
  
  total_abundance <- sum(abundance_table[, sample_id], na.rm = TRUE)
  persistent_abundance <- sum(abundance_table[persistent_votus, sample_id], na.rm = TRUE)
  
  if (total_abundance == 0) return(0)
  return(persistent_abundance / total_abundance)
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

Viral_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")
Viral_metadata_infants <- Viral_metadata[Viral_metadata$Virus_ID %in% rownames(Abundance_table_infants),]
Viral_metadata_mothers <- Viral_metadata[Viral_metadata$Virus_ID %in% rownames(Abundance_table_mothers),]
  

#*******************************
# 1. Estimate persistence in infants and mothers
#*******************************

# First, we identify persistent vOTUs in infants and in mothers to relate it to phage traits
# Persistent vOTU: 
# -- Infants: vOTU detected in at least 3 timepoints of the same infant, with at least 1 early (W2,M1,M2,M3) and one late (M6, M9, M12) timepoints
# -- Mothers: vOTU detected in at least 1 pregnancy (P12,P28,B) and 1 post-pregnancy timepoint (M3) of the same mother
# Next, we estimate the proportion of persistent vOTUs for each infant and mother (viral stability) to relate it with human phenotypes
# -- We estimate this proportion also per vOTU host (top 5) and DGR+ phages.

###############
# A. INFANTS
###############

#---------------------------------------
# A.1. First, identify persistent vOTUs
#---------------------------------------

# Identify infant vOTUs
infant_vOTUs <- rownames(Abundance_table_infants)[rowSums(Abundance_table_infants != 0) > 0]

persistent_vOTUs <- c()
num_infants_persistent <- c()

for (vOTU in infant_vOTUs) {
  # Get samples where the vOTU is present and retrieve NEXT IDs (infants)
  samples_present <- colnames(Abundance_table_infants)[which(Abundance_table_infants[vOTU, ] != 0)]
  infants_present <- unique(Sample_metadata$NEXT_ID[Sample_metadata$NG_ID %in% samples_present])
  
  # Count how many infants have this vOTU persistently
  persistent_count <- sum(sapply(infants_present, function(infant) {
    timepoints_present <- unique(Sample_metadata$Timepoint[Sample_metadata$NG_ID %in% samples_present & Sample_metadata$NEXT_ID == infant])
    
    # Check if there are at least 3 timepoints
    has_enough_timepoints <- length(timepoints_present) >= 3
    
    # Check if there is at least one early and one late timepoint
    has_early_timepoint <- any(timepoints_present %in% c("W2", "M1", "M2", "M3"))
    has_late_timepoint <- any(timepoints_present %in% c("M6", "M9", "M12"))
    
    # Persistent condition: At least 3 timepoints, with one early and one late timepoint
    has_enough_timepoints && has_early_timepoint && has_late_timepoint
  }))
  
  # If vOTU is persistent in at least one infant, store it
  if (persistent_count > 0) {
    persistent_vOTUs <- c(persistent_vOTUs, vOTU)
    num_infants_persistent <- c(num_infants_persistent, persistent_count)
  }
}

# Add persistency results to to Viral metadata (n=2,238)
Viral_metadata_infants$Persistent <- ifelse(
  Viral_metadata_infants$Virus_ID %in% persistent_vOTUs,
  "Yes", "No"
)

# Create a named vector mapping vOTUs to their persistence count
vOTU_infant_counts <- setNames(num_infants_persistent, persistent_vOTUs)

# Add the number of infants in which a vOTU is persistent as a new column to Viral metadata
Viral_metadata_infants$N_Infants_Persistent <- vOTU_infant_counts[Viral_metadata_infants$Virus_ID]

# Estimate the median/range infants in which persistent vOTUs were "persistent"
Viral_metadata_infants_pers <- Viral_metadata_infants[Viral_metadata_infants$Persistent == "Yes",]
median(Viral_metadata_infants_pers$N_Infants_Persistent, na.rm = TRUE)
range(Viral_metadata_infants_pers$N_Infants_Persistent, na.rm = TRUE)
quantile(Viral_metadata_infants_pers$N_Infants_Persistent, probs = c(0.25, 0.75), na.rm = TRUE)

# Estimate the number of vOTUs persisting in at least 10 infants
length(which(Viral_metadata_infants$N_Infants_Persistent >=10))

# Replace NA values with 0 (for non-persistent vOTUs)
Viral_metadata_infants$N_Infants_Persistent[is.na(Viral_metadata_infants$N_Infants_Persistent)] <- 0


#---------------------------------------------------------
# A.2. Estimate proportion of persistent vOTUs 
#---------------------------------------------------------

# Convert abundance to binary (present/absent)
abundance_binary <- Abundance_table_infants != 0

# Merge with metadata to get timepoints per sample
sample_timepoints <- Sample_metadata[, c("NG_ID", "NEXT_ID", "Timepoint_categorical")] 
abundance_with_timepoints <- merge(
  data.frame(NG_ID = colnames(abundance_binary), t(abundance_binary), check.names = FALSE),
  sample_timepoints, 
  by = "NG_ID"
)

# Aggregate timepoints per vOTU per infant across all samples
timepoints_per_votu_infant <- abundance_with_timepoints %>%
  pivot_longer(cols = -c(NG_ID, NEXT_ID, Timepoint_categorical), names_to = "vOTU", values_to = "Present") %>%
  filter(Present) %>%
  group_by(NEXT_ID, vOTU) %>%
  summarise(Timepoints = list(unique(Timepoint_categorical)), .groups = "drop")

# Define persistence criteria
early_timepoints <- c("W2", "M1", "M2", "M3")
late_timepoints <- c("M6", "M9", "M12")

# Calculate overall, genus-specific, and DGR-mediated persistence metrics
infant_viral_persistence <- timepoints_per_votu_infant %>%
  left_join(Viral_metadata_infants[, c("Virus_ID", "Bacterial_genus_host", "DGR")], 
            by = c("vOTU" = "Virus_ID")) %>%
  mutate(Is_Persistent = sapply(Timepoints, function(tp) {
    length(tp) >= 3 && any(tp %in% early_timepoints) && any(tp %in% late_timepoints)
  })) %>%
  group_by(NEXT_ID) %>%
  summarise(
    Total_vOTUs = n_distinct(vOTU),
    Num_Persistent_vOTUs = sum(Is_Persistent),
    Prop_Persistent_vOTUs = ifelse(Total_vOTUs > 0, Num_Persistent_vOTUs / Total_vOTUs, 0)) %>%
  distinct(NEXT_ID, .keep_all = TRUE)

# Add "Prop_Persistent_vOTUs" to Sample metadata (Note that all samples from the same individual will have the same value)
Sample_metadata_infants <- Sample_metadata_infants %>%
  left_join(infant_viral_persistence[, c("NEXT_ID", "Prop_Persistent_vOTUs")], by = "NEXT_ID")

#---------------------------------------------------------
# A.3. Estimate relative abundance of persistent vOTUs per sample
#---------------------------------------------------------

# Define early and late timepoints
early_tp <- c("W2","M1","M2","M3")
late_tp <- c("M6","M9","M12")

# For this, first we retrieve the list of persistent vOTUs per infant
persistent_votus_by_infant <- timepoints_per_votu_infant %>%
  mutate(Is_Persistent = sapply(Timepoints, function(tp) {
    length(tp) >= 3 &&
      any(tp %in% early_tp) &&
      any(tp %in% late_tp)
  })) %>%
  filter(Is_Persistent) %>%
  group_by(NEXT_ID) %>%
  summarise(Persistent_vOTUs = list(unique(vOTU)), .groups = "drop")

persistent_votu_list_infants <- setNames(persistent_votus_by_infant$Persistent_vOTUs,
                                         persistent_votus_by_infant$NEXT_ID)


# Add to Sample metadata
Sample_metadata_infants$Persistent_Abundance_Fraction <- sapply(
  Sample_metadata_infants$NG_ID,
  compute_persistent_abundance,
  abundance_table = Abundance_table_infants,
  metadata = Sample_metadata_infants,
  persistent_list = persistent_votu_list_infants  
)

Sample_metadata_infants$Persistent_Abundance_Fraction[is.na(Sample_metadata_infants$Persistent_Abundance_Fraction)] <- 0


###############
# B. MOTHERS
###############

#---------------------------------------
# B.1. First, identify persistent vOTUs
#---------------------------------------

# Identify mother vOTUs
mother_vOTUs <- rownames(Abundance_table_mothers)[rowSums(Abundance_table_mothers != 0) > 0]

# Initialize storage for persistent vOTUs and their mother counts
persistent_vOTUs_mothers <- c()
num_mothers_persistent <- c()

for (vOTU in mother_vOTUs) {
  # Get samples where the vOTU is present and retrieve NEXT IDs (mothers)
  samples_present <- colnames(Abundance_table_mothers)[which(Abundance_table_mothers[vOTU, ] != 0)]
  mothers_present <- unique(Sample_metadata$NEXT_ID[Sample_metadata$NG_ID %in% samples_present])
  
  # Count how many mothers have this vOTU persistently
  persistent_count <- sum(sapply(mothers_present, function(mother) {
    timepoints_present <- unique(Sample_metadata$Timepoint[Sample_metadata$NG_ID %in% samples_present & Sample_metadata$NEXT_ID == mother])
    has_pregnancy_timepoint <- any(timepoints_present %in% c("P12", "P28"))
    has_post_pregnancy_timepoint <- any(timepoints_present %in% c("M3"))
    has_pregnancy_timepoint && has_post_pregnancy_timepoint  # TRUE (1) if persistent, FALSE (0) otherwise
  }))
  
  # If vOTU is persistent in at least one mother, store it
  if (persistent_count > 0) {
    persistent_vOTUs_mothers <- c(persistent_vOTUs_mothers, vOTU)
    num_mothers_persistent <- c(num_mothers_persistent, persistent_count)
  }
}

# Add persistency results to to Viral metadata
Viral_metadata_mothers$Persistent <- ifelse(
  Viral_metadata_mothers$Virus_ID %in% persistent_vOTUs_mothers,
  "Yes", "No"
)

# Create a named vector mapping vOTUs to their persistence count
vOTU_mother_counts <- setNames(num_mothers_persistent, persistent_vOTUs_mothers)

# Add the number of mothers in which a vOTU is persistent as a new column to Viral metadata
Viral_metadata_mothers$N_Mothers_Persistent <- vOTU_mother_counts[Viral_metadata_mothers$Virus_ID]

# Estimate the median/range mothers in which persistent vOTUs were "persistent"
Viral_metadata_mothers_pers <- Viral_metadata_mothers[Viral_metadata_mothers$Persistent == "Yes",]
median(Viral_metadata_mothers_pers$N_Mothers_Persistent, na.rm = TRUE)
range(Viral_metadata_mothers_pers$N_Mothers_Persistent, na.rm = TRUE)
quantile(Viral_metadata_mothers_pers$N_Mothers_Persistent, probs = c(0.25, 0.75), na.rm = TRUE)

# Estimate the number of vOTUs persisting in at least 10 mothers
length(which(Viral_metadata_mothers$N_Mothers_Persistent >=10))

# Replace NA values with 0 (for non-persistent vOTUs)
Viral_metadata_mothers$N_Mothers_Persistent[is.na(Viral_metadata_mothers$N_Mothers_Persistent)] <- 0

#---------------------------------------------------------
# A.2. Estimate proportion of persistent vOTUs
#---------------------------------------------------------

# Convert abundance to binary (present/absent) for mothers
abundance_binary_mothers <- Abundance_table_mothers != 0

# Merge with metadata to get timepoints per sample
sample_timepoints <- Sample_metadata[, c("NG_ID", "NEXT_ID", "Timepoint_categorical")]
abundance_with_timepoints_mothers <- merge(
  data.frame(NG_ID = colnames(abundance_binary_mothers), t(abundance_binary_mothers), check.names = FALSE),
  sample_timepoints, 
  by = "NG_ID"
)

# Aggregate timepoints per vOTU per mother across all samples
timepoints_per_votu_mother <- abundance_with_timepoints_mothers %>%
  pivot_longer(cols = -c(NG_ID, NEXT_ID, Timepoint_categorical), names_to = "vOTU", values_to = "Present") %>%
  filter(Present) %>%
  group_by(NEXT_ID, vOTU) %>%
  summarise(Timepoints = list(unique(Timepoint_categorical)), .groups = "drop")

# Define persistence criteria for mothers
pregnancy_timepoints <- c("P12", "P28")
post_pregnancy_timepoints <- c("M3")

# Calculate overall and genus-specific persistence metrics (one row per mother)
mother_viral_persistence <- timepoints_per_votu_mother %>%
  left_join(Viral_metadata_infants[, c("Virus_ID", "Bacterial_genus_host", "DGR")], 
            by = c("vOTU" = "Virus_ID")) %>%
  mutate(Is_Persistent = sapply(Timepoints, function(tp) {
    any(tp %in% pregnancy_timepoints) && any(tp %in% post_pregnancy_timepoints)
  })) %>%
  group_by(NEXT_ID) %>%
  summarise(
    Total_vOTUs = n_distinct(vOTU),
    Num_Persistent_vOTUs = sum(Is_Persistent),
    Prop_Persistent_vOTUs = ifelse(Total_vOTUs > 0, Num_Persistent_vOTUs / Total_vOTUs, 0)) %>%
  distinct(NEXT_ID, .keep_all = TRUE)

# Add "Prop_Persistent_vOTUs" to Sample metadata (Note that all samples from the same individual will have the same value)
Sample_metadata_mothers <- Sample_metadata_mothers %>%
  left_join(mother_viral_persistence[, c("NEXT_ID", "Prop_Persistent_vOTUs")], by = "NEXT_ID")

#---------------------------------------------------------
# A.3. Estimate relative abundance of persistent vOTUs per sample
#---------------------------------------------------------

# For this, first we retrieve the list of persistent vOTUs per mother 
persistent_votus_by_mother <- timepoints_per_votu_mother %>%
  mutate(Is_Persistent = sapply(Timepoints, function(tp) {
    any(tp %in% pregnancy_timepoints) && any(tp %in% post_pregnancy_timepoints)
  })) %>%
  filter(Is_Persistent) %>%
  group_by(NEXT_ID) %>%
  summarise(Persistent_vOTUs = list(unique(vOTU)), .groups = "drop")

persistent_votu_list_mothers <- setNames(persistent_votus_by_mother$Persistent_vOTUs, persistent_votus_by_mother$NEXT_ID)

# Add to Sample_metadata
Sample_metadata_mothers$Persistent_Abundance_Fraction <- sapply(
  Sample_metadata_mothers$NG_ID,
  compute_persistent_abundance,
  abundance_table = Abundance_table_mothers,
  metadata = Sample_metadata_mothers,
  persistent_list = persistent_votu_list_mothers
)

Sample_metadata_mothers$Persistent_Abundance_Fraction[is.na(Sample_metadata_mothers$Persistent_Abundance_Fraction)] <- 0

# We add the Persistent_Abundance_Fraction to the general Sample_metadata
Persistent_Abundance_Fraction_mothers <- Sample_metadata_mothers %>%
  select(NG_ID, Persistent_Abundance_Fraction)
Persistent_Abundance_Fraction_infants <- Sample_metadata_infants %>%
  select(NG_ID, Persistent_Abundance_Fraction)
persistent_abundance_all <- bind_rows(Persistent_Abundance_Fraction_mothers, Persistent_Abundance_Fraction_infants)

Sample_metadata <- Sample_metadata %>%
  left_join(persistent_abundance_all, by = "NG_ID")

#*******************************
# 2. Statistical analysis
#*******************************

########################################################
# 2.1. Subsetting mothers and infants (complete sampling)
########################################################

# Note statistical analysis was performed on a subset of mothers (n=174) and infants (n=166) with (almost) complete timepoints.
# Note: For the proportion of persistent vOTUs per individual, no NA values are present (if no persistent vOTUs -> 0)
# For the relative abundance in a sample explained by persistent vOTUs:
# 0 means: the individual has persistent vOTUs, but none are present in that sample 
# Also 0 means (after conversion of NAs to 0):  the individual had no persistent vOTUs at all

# Select subset of mothers and infants with complete timepoints available
all_infant_timepoints <- c( "M1", "M2", "M3","M6", "M9", "M12")
all_mother_timepoints <- c("P12", "P28", "M3")

infants_with_all_timepoints <- Sample_metadata_infants %>%
  filter(Timepoint_categorical %in% all_infant_timepoints) %>%
  group_by(NEXT_ID) %>%
  summarise(n_timepoints = n_distinct(Timepoint_categorical)) %>%
  filter(n_timepoints == length(all_infant_timepoints)) %>%
  pull(NEXT_ID)

mothers_with_all_timepoints <- Sample_metadata_mothers %>%
  filter(Timepoint_categorical %in% all_mother_timepoints) %>%
  group_by(NEXT_ID) %>%
  summarise(n_timepoints = n_distinct(Timepoint_categorical)) %>%
  filter(n_timepoints == length(all_mother_timepoints)) %>%
  pull(NEXT_ID)

persistence_infant_selection <- infant_viral_persistence[infant_viral_persistence$NEXT_ID %in% infants_with_all_timepoints,]
persistence_mother_selection <- mother_viral_persistence[mother_viral_persistence$NEXT_ID %in% mothers_with_all_timepoints,]

NG_IDs_mothers_with_all_timepoints <- Sample_metadata_mothers[Sample_metadata_mothers$NEXT_ID %in% mothers_with_all_timepoints,"NG_ID"]
NG_IDs_infants_with_all_timepoints <- Sample_metadata_infants[Sample_metadata_infants$NEXT_ID %in% infants_with_all_timepoints,"NG_ID"]

########################################################
# 2.2. Test proportion of persistent vOTUs
########################################################

Personal_persistent_virome_prop_mother_infant <- wilcox.test(persistence_mother_selection$Prop_Persistent_vOTUs, persistence_infant_selection$Prop_Persistent_vOTUs)
Personal_persistent_virome_prop_mother_infant$p.value #4.0817e-31
median(persistence_mother_selection$Prop_Persistent_vOTUs) # 0.343
median(persistence_infant_selection$Prop_Persistent_vOTUs) # 0.11

########################################################
# 2.3. Test relative abundance of persistent vOTUs (per sample)
########################################################

# Subset the relative abundance is samples belonging to indivduals from the subset
Persistent_Abundance_Fraction_mothers_test <- Persistent_Abundance_Fraction_mothers %>%
  filter(NG_ID %in% NG_IDs_mothers_with_all_timepoints) %>%
  left_join(Sample_metadata_mothers %>% select(NG_ID, NEXT_ID), 
            by = c("NG_ID" = "NG_ID")) %>%
  mutate(Type = "Mother")

Persistent_Abundance_Fraction_infants_test <- Persistent_Abundance_Fraction_infants %>%
  filter(NG_ID %in% NG_IDs_infants_with_all_timepoints) %>%
  left_join(Sample_metadata_infants %>% select(NG_ID, NEXT_ID), 
            by = c("NG_ID" = "NG_ID")) %>%
  mutate(Type = "Infant")

# Combine both
Persistent_Abundance_Fraction_all_test <- bind_rows(Persistent_Abundance_Fraction_mothers_test,
  Persistent_Abundance_Fraction_infants_test)

# Perform the test (using LMM)
model_persistence_abundance <- lmer(Persistent_Abundance_Fraction ~ Type + (1 | NEXT_ID),
  data = Persistent_Abundance_Fraction_all_test)
summary(model_persistence_abundance)

median_summary <- Persistent_Abundance_Fraction_all_test %>%
  group_by(Type) %>%
  summarise(
    median_Persistent_Abundance_Fraction = median(Persistent_Abundance_Fraction, na.rm = TRUE),
    n = n()
  )

###############
# D. PLOTS 
###############

#---------------------------------------------------------
# 1. Barplot comparing the proportion of persistent phages in mothers and infants
#---------------------------------------------------------

infants_summary <- Viral_metadata_infants %>%
  group_by(Persistent) %>%
  summarise(Count = n()) %>%
  mutate(Group = "Infants", Proportion = Count / sum(Count))

mothers_summary <- Viral_metadata_mothers %>%
  group_by(Persistent) %>%
  summarise(Count = n()) %>%
  mutate(Group = "Mothers", Proportion = Count / sum(Count))

number_both <- bind_rows(infants_summary, mothers_summary)
number_both$Persistent <- ifelse(number_both$Persistent == "Yes", "Persistent", "Non-Persistent")

pdf('11_PERSISTENCE/Plots/vOTU_persistence_infants_mothers.pdf', width = 2.6, height = 3.6)
ggplot(number_both, aes(x = Group, y = Proportion, fill = Persistent)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_text(aes(label = Count),
            position = position_stack(vjust = 0.5), color = "black", size = 4.5) +
  scale_fill_manual(values = c("Persistent" = "#457B9D", "Non-Persistent" = "#DADADA")) +
  labs(x = "Group", y = "Proportion of vOTUs (%)", fill = "Phage persistence") +
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

#---------------------------------------------------------
# 2. Boxplot comparing proportion of persistent voTUs in mothers and infants 
#---------------------------------------------------------

# Combine mother and infant datasets
combined_persistence <- rbind(
  data.frame(Group = "Mother", Prop_Persistent_vOTUs = persistence_mother_selection$Prop_Persistent_vOTUs),
  data.frame(Group = "Infant", Prop_Persistent_vOTUs = persistence_infant_selection$Prop_Persistent_vOTUs)
)

pdf('11_PERSISTENCE/Plots/Proportion_persistent_vOTUs_per_individual_mothers_infants.pdf', width = 2.7, height = 3.2)
ggplot(combined_persistence, aes(x = Group, y = Prop_Persistent_vOTUs, fill = Group)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Group', y = '% Persistent vOTUs') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 0.9), expand = expansion(mult = c(0.05, 0.02))) +
  scale_fill_manual(values = c("Infant" = "#66A6AD", "Mother" = "#8E7CA6")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()


#---------------------------------------------------------
# 3. Boxplot comparing the relative abundance of persistent vOTUs in mothers and infants 
#---------------------------------------------------------

pdf('11_PERSISTENCE/Plots/Persistent_abundance_fraction_vOTUs_per_individual_mothers_infants.pdf', width = 2.8, height = 3.2)
ggplot(Persistent_Abundance_Fraction_all_test, aes(x = Type, y = Persistent_Abundance_Fraction, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Group', y = 'Persistent abundance fraction') +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.25),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),expand = expansion(mult = c(0.05, 0.02))) +
  scale_fill_manual(values = c("Infant" = "#66A6AD", "Mother" = "#8E7CA6")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14, color = "black"),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none"
  )
dev.off()


#****************
# Write results
#****************

# Infants with near-complete longitudinal sampling
write.table(infants_with_all_timepoints,"11_PERSISTENCE/Infants_complete_sampling.txt", sep = "\t", row.names = F, quote = FALSE) 

# Viral metadata
write.table(Viral_metadata,"Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Viral_metadata_infants,"Metadata_NEXT/LLNEXT_Viral_Metadata_INFANTS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Viral_metadata_mothers,"Metadata_NEXT/LLNEXT_Viral_Metadata_MOTHERS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(infant_viral_persistence,"11_PERSISTENCE/LLNEXT_Infant_Viral_Persistance_data.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(mother_viral_persistence,"11_PERSISTENCE/LLNEXT_Mother_Viral_Persistance_data.txt", sep = "\t", row.names = F, quote = FALSE) 

# Sample metadata
write.table(Sample_metadata,"Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Sample_metadata_infants,"Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Sample_metadata_mothers,"Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 

