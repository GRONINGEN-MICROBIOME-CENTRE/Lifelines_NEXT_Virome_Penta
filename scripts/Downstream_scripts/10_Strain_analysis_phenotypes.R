################################################################################
##### LL-NEXT: Viral strain transmission - Phenotypes
### Author(s): Asier Fernández-Pato
### Last updated: 5th December, 2025
################################################################################

#****************
# Load libraries
#****************
library(dplyr)
library(data.table)
library(scales)
library(stringr)
library(tidyr)
library(lme4)
library(lmerTest)
library(ggplot2)

# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_Analysis/VIRUSES/")


##************************************************************************
# 1. Load metadata, abundance table and inStrain results for the LL-NEXT samples 
#*************************************************************************
# Read metadata tables 
Sample_metadata <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt")
Sample_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_07052025.txt")
Sample_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_phenos_07052025.txt")

# Read abundance tables
Abundance_table <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_07052025.txt")
Abundance_table_infants <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_INFANTS_07052025.txt")
Abundance_table_mothers <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_MOTHERS_07052025.txt")

# Read viral metadata tables
Viral_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")
Viral_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_INFANTS_07052025.txt")
Viral_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_MOTHERS_07052025.txt")

# Read processed inStrain results 
inStrain_results <- read.delim("10_STRAIN_TRANSMISSION/inStrain_results_processed.txt")

#*******************************
# 2. Analysis of AGE
#*******************************

#******************************************************************
# 2.A. Check how mother-infant strain sharing (in pairs) changes over time 
#******************************************************************

#################################################################################
# PLOT: proportion of Mother-infant pairs that have strain sharing per timepoint
#################################################################################

#Note: Here we are exploring sharing per pair (1st and 2nd pregnancies and twins are considered as separated pairs)

# Get the number of infants per family and timepoint
infants_per_timepoint <- Sample_metadata %>%
  filter(Type == "Infant") %>%
  group_by(Timepoint_categorical, FAMILY) %>%
  summarize(Infant_NEXT_IDs = n_distinct(NEXT_ID), .groups = "drop")

# Get the number of mothers per family
mothers_per_family <- Sample_metadata %>%
  filter(Type == "Mother") %>%
  group_by(FAMILY) %>%
  summarize(Mother_NEXT_IDs = n_distinct(NEXT_ID), .groups = "drop")

# Merge infant and mother data by FAMILY
pairs_per_timepoint <- infants_per_timepoint %>%
  left_join(mothers_per_family, by = "FAMILY") %>%
  filter(!is.na(Mother_NEXT_IDs)) %>%
  mutate(
    # Calculate the number of pairs for each family at each timepoint
    Total_Mother_Infant_Pairs = Infant_NEXT_IDs * Mother_NEXT_IDs
  )

# Summarize the total number of pairs per timepoint
all_pairs_per_timepoint <- pairs_per_timepoint %>%
  group_by(Timepoint_categorical) %>%
  summarize(Total_Pairs = sum(Total_Mother_Infant_Pairs), .groups = "drop")

# Calculate families with inStrain results per timepoint
compared_pairs_per_timepoint <- inStrain_results %>%
  filter(Mother_Infant_pair == "Pair") %>%
  group_by(Infant_timepoint) %>%
  summarize(Compared_Pairs = length(unique(Pair_ID))) 

# Merge the two datasets
family_sharing_per_timepoint <- all_pairs_per_timepoint  %>%
  left_join(compared_pairs_per_timepoint, by = c("Timepoint_categorical" = "Infant_timepoint")) %>%  
  mutate(Compared_Pairs = ifelse(is.na(Compared_Pairs), 0, Compared_Pairs))  

# Count the number of mother-infant pairs with strain sharing per infant timepoint 
sharing_per_timepoint <- inStrain_results %>%
  filter(Mother_Infant_pair == "Pair" & Mother_infant_sharing == "Yes") %>%
  group_by(Infant_timepoint) %>%
  summarize(Sharing_Families = n_distinct(Pair_ID))

# Add results to families_per_timepoint and estimate proportion
family_sharing_per_timepoint <- family_sharing_per_timepoint %>%
  left_join(sharing_per_timepoint, by = c("Timepoint_categorical" = "Infant_timepoint")) %>%
  mutate(Proportion_Sharing = Sharing_Families / Total_Pairs)

# Ensure Timepoint_categorical is a factor with the correct order
timepoint_order <- c("W2", paste0("M", 1:12))  
family_sharing_per_timepoint$Timepoint_categorical <- factor(family_sharing_per_timepoint$Timepoint_categorical, levels = timepoint_order)

# Generate the bar plot (proportion of families per timepoint with strain sharing)
pdf('10_STRAIN_TRANSMISSION/Plots/Proportion_families_with_viral_strain_sharing_per_timepoint.pdf', width = 4.1, height = 3.9)
timepoint_colors <- c("#4B4FC5", "#3B81B3", "#4CA6B1", "#88CFA4", "#C1D97F", "#E8C26D", "#F4A219")

ggplot(family_sharing_per_timepoint, aes(x = Timepoint_categorical, y = Proportion_Sharing, fill = Timepoint_categorical)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_text(aes(label = paste0("n=", Sharing_Families)), 
            position = position_stack(vjust = 0.5), color = "black", size = 4.55) +
  scale_fill_manual(values = timepoint_colors) +
  labs(x = "Timepoint", y = "Proportion of families with strain sharing (%)", fill = "Timepoint") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,0.52)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()


######################################################################################
# PLOT: Family-specific sharing rate: Proportion of shared infant vOTUs per timepoint 
######################################################################################

# Generate inStrain results df only for pairs
inStrain_results_pairs <-inStrain_results[inStrain_results$Mother_Infant_pair == "Pair",]

# Count the number of pairs in which strains coould be profiled and compared
length(unique(inStrain_results_pairs$Pair_ID)) # 432
  
# Generate a dataframe with the number of vOTUs shared per pair (Pair_ID) and per timepoint (related pairs)
vOTU_sharing_per_family <- inStrain_results_pairs %>%
  mutate(FAM_ID = str_extract(Pair_ID, "^[^_]+")) %>% 
  group_by(Infant_timepoint, Pair_ID) %>%  
  summarize(
    Num_vOTUs = n_distinct(Virus_ID),
    Num_vOTUs_shared = n_distinct(Virus_ID[Mother_infant_sharing == "Yes"]),  # Unique Virus_IDs where sharing is "Yes"
    .groups = "drop" 
  ) 

# Identify the total number of vOTUs present in the infant sample for each pair and timepoint
# Note: For some infants, multiple samples of the same timepoint are available (we consider them together)
vOTU_sharing_per_family <- vOTU_sharing_per_family %>%
  rowwise() %>%
  mutate(
    Num_vOTUs_present_infant = {
      # Extract the timepoint, Pair ID and FAM_ID for the current row
      timepoint <- Infant_timepoint
      pair <- Pair_ID
      family <- sub("_.*", "", pair)
      infant_NEXT_ID <- inStrain_results_pairs %>%
        filter(Pair_ID == pair) %>%
        pull(Infant_NEXT_ID) %>%
        unique()
      
      # Find the SAMPLE_IDs corresponding to Infant NEXT_ID in Sample_metadata_infants
      # Same infants have more than 1 sample at the same timepoint
      ng_ids <- Sample_metadata_infants %>%
        filter(Timepoint_categorical == timepoint & NEXT_ID == infant_NEXT_ID) %>%
        pull(NG_ID)
      
      # Calculate the total number of vOTUs present for each infant and each timepoint
      total_vOTUs <- length(unique(unlist(sapply(ng_ids, function(ng_id) {
        which(Abundance_table_infants[[ng_id]] > 0)
      }))))
      
      # Return the total number of vOTUs present
      total_vOTUs
    }
  )

# Add additional phenotypes to vOTU_sharing_per_family 
phenotypes <- c("Feeding_mode", "Delivery_mode", "Delivery_place",
                "Infant_eczema", "Infant_food_allergy")

for (phenotype in phenotypes) {
  vOTU_sharing_per_family[[phenotype]] <- apply(vOTU_sharing_per_family, 1, function(row) {
    pair <- row["Pair_ID"]
    
    # Get phenotype values from inStrain_results_pairs
    subset_row <- inStrain_results_pairs[inStrain_results_pairs$Pair_ID == pair, ]
    
    # Extract phenotype values
    values <- unique(subset_row[[phenotype]])

    # Remove NA
    values <- values[!is.na(values)]
    
    if (length(values) == 0) {
      return(NA)
    } else {
      return(values[1])  # take the first non-NA value
    }
  })
}

# Estimate the proportion of infant vOTUs shared per family
vOTU_sharing_per_family$Proportion_sharing <- vOTU_sharing_per_family$Num_vOTUs_shared /vOTU_sharing_per_family$Num_vOTUs_present_infant
vOTU_sharing_per_family$Infant_timepoint <- factor(vOTU_sharing_per_family$Infant_timepoint, levels = timepoint_order)

# Get average sharing rate at W2
W2_sharing_rate <- as.numeric(unlist(vOTU_sharing_per_family[vOTU_sharing_per_family$Infant_timepoint == "W2",
                                                             "Proportion_sharing"]))
mean(W2_sharing_rate) * 100 # 7.856%

# Generate the plot:
# The Family sharing rate including only those families for which strain comparisons results are available
pdf('10_STRAIN_TRANSMISSION/Plots/Family_sharing_rate_families_with_strain_results_per_timepoint.pdf', width = 4.2, height = 3.9)
ggplot(vOTU_sharing_per_family, aes(x = Infant_timepoint, y = Proportion_sharing, fill = Infant_timepoint)) +
  geom_jitter(aes(size = Num_vOTUs_shared), alpha = 0.4, color = "lightgrey", width = 0.2) +  
  geom_boxplot(width = 0.6, color = "black", outlier.shape = NA, size = 0.9) +
  scale_fill_manual(values = timepoint_colors) +
  scale_size_continuous(name = "Number of vOTUs", range = c(1, 4)) +
  labs(x = "Timepoint", y = "vOTU sharing rate (%)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
dev.off()

###################
# STATISTICAL TEST  
###################

# A) Proportion of families with viral strain sharing over time
# Here all 666 families are included (as detection of sharing is being compared and not sharing itself)
family_sharing <- inStrain_results_pairs %>%
  group_by(Infant_timepoint, Infant_FAM_ID) %>%
  summarise(
    Any_Sharing = ifelse(any(Mother_infant_sharing == "Yes"), 1, 0),
    .groups = "drop"
  )

full_family_sharing <- pairs_per_timepoint %>%
  left_join(family_sharing,
            by = c("Timepoint_categorical" = "Infant_timepoint",
                   "FAMILY" = "Infant_FAM_ID")) %>%
  mutate(
    Any_Sharing = replace_na(Any_Sharing, 0) # add 0 to families without sharing
  )

# Test significance using a GLMM
full_family_sharing$Timepoint_categorical <- factor(
  full_family_sharing$Timepoint_categorical,
  levels = c("W2", "M1", "M12", "M2", "M3", "M6", "M9") # W2 first
)

GLMM_family_sharing_time <- glmer(
  Any_Sharing ~ Timepoint_categorical + (1 | FAMILY),
  data = full_family_sharing,
  family = binomial
)

summary(GLMM_family_sharing_time)

# B) vOTU sharing rate over time
# Add FAM_ID and timepoint in months to vOTU_sharing_per_family
vOTU_sharing_per_family <- vOTU_sharing_per_family %>%
  mutate(FAM_ID = str_extract(Pair_ID, "^FAM\\d+"))

vOTU_sharing_per_family <- vOTU_sharing_per_family %>%
  left_join(
    inStrain_results %>%
      select(Pair_ID, Infant_timepoint, Infant_timepoint_months) %>%
      distinct(),
    by = c("Pair_ID", "Infant_timepoint")
  )

# Test significance using a GLMM
GLMM_vOTU_sharing_time <- glmer(
  cbind(Num_vOTUs_shared, Num_vOTUs_present_infant - Num_vOTUs_shared) ~ 
    Infant_timepoint_months + (1 | FAM_ID),
  data = vOTU_sharing_per_family,
  family = binomial
)

summary(GLMM_vOTU_sharing_time)


#*******************************
# 3. Analysis of DELIVERY MODE
#*******************************

#******************************************************************
# 3.A. Check if mother-infant pairs with CS/VG have more similar strains
#******************************************************************

########################
# STATISTICAL ANALYSIS
########################

# Check if CS/VG born pairs have more similar strains
inStrain_results_pairs$Delivery_mode <- factor(inStrain_results_pairs$Delivery_mode,
                                                                    levels = c("VG", "CS"))

# Perform statistical test (use permutation-based test due to lack of independence between infants)
num_permutations <- 10000
permuted_stats <- numeric(num_permutations)
set.seed(123) 

# Extract the observed test statistic: difference in means between VG and CS delivery modes
observed_stat <- mean(inStrain_results_pairs$popANI[inStrain_results_pairs$Delivery_mode == "VG"], na.rm = TRUE) - 
  mean(inStrain_results_pairs$popANI[inStrain_results_pairs$Delivery_mode == "CS"], na.rm = TRUE)

for (i in 1:num_permutations) {
  permuted_labels <- sample(inStrain_results_pairs$Delivery_mode)
  permuted_stat <- mean(inStrain_results_pairs$popANI[permuted_labels == "VG"], na.rm = TRUE) - 
    mean(inStrain_results_pairs$popANI[permuted_labels == "CS"], na.rm = TRUE)
  permuted_stats[i] <- permuted_stat
}

# Calculate the p-value for a two-sided test
# We compare the absolute value of the observed test statistic to the absolute values of the permuted statistics
p_value <- mean(abs(permuted_stats) >= abs(observed_stat))


########
# PLOT
########

# Generate boxplot of strain similarity in pairs vs non-pairs
pdf('10_STRAIN_TRANSMISSION/Plots/Mother_infant_related_distances_delivery.pdf', width = 2.8, height = 3.2)
ggplot(inStrain_results_pairs[!is.na(inStrain_results_pairs$Delivery_mode),], 
       aes(x = Delivery_mode, y = popANI, fill = Delivery_mode)) +
  geom_jitter(alpha = 0.4, aes(color = Delivery_mode), size = 1.8, width = 0.2) + 
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) + 
  scale_fill_manual(values = c("VG" = "#B6D7A8", "CS" = "#A2C2E5")) +
  scale_color_manual(values = c("VG" = "#B6D7A8", "CS" = "#A2C2E5")) +
  labs(x = "Delivery mode", y = "Strain similarity (popANI)", fill = "Delivery mode") +
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

#******************************************************************
# 3.B. Check if pairs with CS/VG have more strain sharing
#******************************************************************

########################
# STATISTICAL ANALYSIS
########################

# Count potential sharing events in CS/VG pairs
# Build contingency table and test significance
mother_infant_sharing_delivery_mode <- na.omit(inStrain_results_pairs %>%
                                                 group_by(Mother_infant_comparison_ID, Delivery_mode) %>%
                                                 summarise(Mother_infant_sharing= ifelse(any(Mother_infant_sharing== "Yes"), "Yes", "No"),
                                                           .groups = "drop"))

contingency_table_delivery_mode <- table(mother_infant_sharing_delivery_mode$Delivery_mode, mother_infant_sharing_delivery_mode$Mother_infant_sharing)

# Perform chi-square test
chi2_test <- chisq.test(contingency_table_delivery_mode)
chi2_test
chi2_test$p.value # 0.0001715

########
# PLOT
########

# Transform contingency table into data frame for ggplot
contingency_table_delivery_mode_plot <- as.data.frame(as.table(contingency_table_delivery_mode))
colnames(contingency_table_delivery_mode_plot) <- c("Delivery_mode", "Mother_infant_sharing", "Count")

# Calculate proportions for stacked bar plot
contingency_table_delivery_mode_plot <- contingency_table_delivery_mode_plot %>%
  group_by(Delivery_mode) %>%
  mutate(Proportion = Count / sum(Count),
         Total = sum(Count))

# Generate stacked bar plot with proportions and counts
pdf('10_STRAIN_TRANSMISSION/Plots/Mother_infant_sharing_delivery.pdf', width = 2.5, height = 3.6)
ggplot(contingency_table_delivery_mode_plot, 
       aes(x = Delivery_mode, y = Proportion, fill = Mother_infant_sharing)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_text(aes(label = paste0("n=", Count)), 
            position = position_stack(vjust = 0.5), 
            color = "black", size = 4.55) +
  scale_fill_manual(values = c("No" = "#E0E0E0", "Yes" = "#A0A0A0")) +
  labs(x = "Infant Delivery Mode", y = "Proportion of sharing events (%)", fill = "Sharing") +
  scale_y_continuous(labels = scales::percent,limits = c(0, 1.15), 
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(add = c(0.05, 0.05))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()

#******************************************************************
# 3.C Check how mother-infant strain sharing (in pairs) changes with delivery mode
#******************************************************************

#################################################################################
# PLOT: proportion of Mother-infant pairs that have strain sharing per timepoint
#################################################################################

# Estimate the number of families with CS or VG per timepoint
family_sharing_per_timepoint <- family_sharing_per_timepoint %>%
  mutate(
    All_Families_CS = sapply(Timepoint_categorical, function(tp) {
      sum(
        Sample_metadata_infants$delivery_mode == "CS" &
          as.character(Sample_metadata_infants$Timepoint_categorical) == as.character(tp) &
          Sample_metadata_infants$FAMILY %in% infants_per_timepoint$FAMILY[infants_per_timepoint$Timepoint_categorical == tp],
        na.rm = TRUE
      )
    }),
    All_Families_VG = sapply(Timepoint_categorical, function(tp) {
      sum(
        Sample_metadata_infants$delivery_mode == "VG" &
          as.character(Sample_metadata_infants$Timepoint_categorical) == as.character(tp) &
          Sample_metadata_infants$FAMILY %in% infants_per_timepoint$FAMILY[infants_per_timepoint$Timepoint_categorical == tp],
        na.rm = TRUE
      )
    })
  )

# Compute the number of sharingfamilies for CS
transmission_per_timepoint_cs <- inStrain_results_pairs %>%
  filter(Mother_infant_sharing== "Yes" & Delivery_mode == "CS") %>%
  group_by(Infant_timepoint) %>%
  summarize(Transmission_Families_CS = n_distinct(Pair_ID))

# Compute the number of sharing families for VG
transmission_per_timepoint_vg <- inStrain_results_pairs %>%
  filter(Mother_infant_sharing== "Yes" & Delivery_mode == "VG") %>%
  group_by(Infant_timepoint) %>%
  summarize(Transmission_Families_VG = n_distinct(Pair_ID))

# Merge the results with the families_per_timepoint table
family_sharing_per_timepoint  <- family_sharing_per_timepoint  %>%
  left_join(transmission_per_timepoint_cs, by = c("Timepoint_categorical" = "Infant_timepoint")) %>%
  left_join(transmission_per_timepoint_vg, by = c("Timepoint_categorical" = "Infant_timepoint")) %>%
  mutate(Proportion_Transmission_CS = Transmission_Families_CS / All_Families_CS) %>%
  mutate(Proportion_Transmission_VG = Transmission_Families_VG / All_Families_VG)

# Generate the bar plot (proportion of families per timepoint with strain sharing by delivery mode)
family_sharing_per_timepoint$Timepoint_categorical <- factor(family_sharing_per_timepoint$Timepoint_categorical,
                                                             levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12"))

family_sharing_plot <- family_sharing_per_timepoint %>%
  pivot_longer(
    cols = starts_with("Proportion_Transmission_"), 
    names_to = "Delivery_Mode", 
    values_to = "Proportion"
  ) %>%
  mutate(
    Delivery_Mode = recode(Delivery_Mode, 
                           "Proportion_Transmission_CS" = "CS", 
                           "Proportion_Transmission_VG" = "VG"),
    Delivery_Mode = factor(Delivery_Mode, levels = c("VG", "CS"))
  )

pdf('10_STRAIN_TRANSMISSION/Plots/Proportion_families_with_viral_strain_sharing_per_timepoint_delivery.pdf', width = 3.6, height = 3.2)
ggplot(family_sharing_plot, aes(x = Timepoint_categorical, y = Proportion, fill = Delivery_Mode)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.9) +
  geom_text(aes(label = paste0("n=", ifelse(Delivery_Mode == "CS", Transmission_Families_CS, Transmission_Families_VG))),
            position = position_dodge(width = 0.8), size = 4, vjust = -0.4, color = "black") +
  scale_fill_manual(values = c("VG" = "#B6D7A8", "CS" = "#A2C2E5")) +
  labs(x = "Timepoint", y = "Families with sharing (%)", fill = "Delivery Mode") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.6)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "top",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
dev.off()

######################################################################################
# PLOT: Family-specific sharing rate: Proportion of shared infant vOTUs per timepoint 
######################################################################################

# Estimate the proportion of infant vOTUs shared per family and delivery mode
vOTU_sharing_per_family$Infant_timepoint <- factor(vOTU_sharing_per_family$Infant_timepoint,
                                                             levels = timepoint_order)

pdf('10_STRAIN_TRANSMISSION/Plots/Family_sharing_rate_delivery_per_timepoint.pdf', width = 5.5, height = 3.9)
ggplot(vOTU_sharing_per_family[!is.na(vOTU_sharing_per_family$Delivery_mode),], 
       aes(x = Infant_timepoint, y = Proportion_sharing, fill = Delivery_mode)) +
  geom_jitter(aes(size = Num_vOTUs_shared), alpha = 0.4, color = "lightgrey",
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +  
  geom_boxplot(aes(fill = Delivery_mode), width = 0.8, color = "black", outlier.shape = NA, size = 0.6) +
  scale_fill_manual(values = c("VG" = "#B6D7A8", "CS" = "#A2C2E5")) + 
  scale_size_continuous(name = "Number of vOTUs", range = c(1, 4)) +
  labs(x = "Timepoint", y = "vOTU sharing rate (%)", fill = "Delivery Mode") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.5)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "bottom",
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )
dev.off()

###################
# STATISTICAL TEST  
###################

# A) Proportion of families with viral strain sharing over time

# Add delivery mode to full_family_sharing
full_family_sharing <- full_family_sharing %>%
  left_join(
    Sample_metadata_infants %>%
      select(FAMILY, delivery_mode) %>%
      distinct(),
    by = "FAMILY"
  )

full_family_sharing$delivery_mode <- factor(
  full_family_sharing$delivery_mode,
  levels = c("CS", "VG")
)

# GLMM testing effect of delivery mode
GLMM_family_sharing_DM <- glmer(
  Any_Sharing ~ delivery_mode + (1 | FAMILY),
  data = full_family_sharing,
  family = binomial
)

summary(GLMM_family_sharing_DM) # delivery_modeVG   1.9085     0.3645   5.235 1.65e-07 ***

# B) vOTU sharing rate over time
vOTU_sharing_per_family$Delivery_mode <- factor(
  vOTU_sharing_per_family$Delivery_mode,
  levels = c("CS", "VG")
)

# GLMM testing effect of delivery mode
GLMM_vOTU_sharing_DM <- glmer(
  cbind(Num_vOTUs_shared, Num_vOTUs_present_infant - Num_vOTUs_shared) ~ 
    Delivery_mode + (1 | FAM_ID),
  data = vOTU_sharing_per_family,
  family = binomial
)

summary(GLMM_vOTU_sharing_DM) #Delivery_modeVG   0.3556     0.1665   2.136   0.0327 * 


#*******************************
# 4. Analysis of DELIVERY PLACE
#*******************************

#******************************************************************
# 4.A. Check if mother-infant pairs delivered at home/hospital have more similar strains
#******************************************************************

########################
# STATISTICAL ANALYSIS
########################

inStrain_results_pairs$Delivery_place <- factor(
  inStrain_results_pairs$Delivery_place,
  levels = c("hospital", "home")
)

num_permutations <- 10000
permuted_stats <- numeric(num_permutations)
set.seed(123)

observed_stat <- mean(inStrain_results_pairs$popANI[inStrain_results_pairs$Delivery_place == "home"], na.rm = TRUE) -
  mean(inStrain_results_pairs$popANI[inStrain_results_pairs$Delivery_place == "hospital"], na.rm = TRUE)

for (i in 1:num_permutations) {
  permuted_labels <- sample(inStrain_results_pairs$Delivery_place)
  permuted_stat <- mean(inStrain_results_pairs$popANI[permuted_labels == "home"], na.rm = TRUE) -
    mean(inStrain_results_pairs$popANI[permuted_labels == "hospital"], na.rm = TRUE)
  permuted_stats[i] <- permuted_stat
}

p_value <- mean(abs(permuted_stats) >= abs(observed_stat))


########
# PLOT
########

pdf('10_STRAIN_TRANSMISSION/Plots/Mother_infant_related_distances_delivery_place.pdf', width = 2.8, height = 3.2)
ggplot(inStrain_results_pairs[!is.na(inStrain_results_pairs$Delivery_place),],
       aes(x = Delivery_place, y = popANI, fill = Delivery_place)) +
  geom_jitter(alpha = 0.4, aes(color = Delivery_place), size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  scale_fill_manual(values = c("Hospital" = "#A2C2E5", "Home" = "#B6D7A8")) +
  scale_color_manual(values = c("Hospital" = "#A2C2E5", "Home" = "#B6D7A8")) +
  labs(x = "Delivery place", y = "Strain similarity (popANI)", fill = "Delivery place") +
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

#******************************************************************
# 4.B. Check if pairs delivered at home/hospital have more strain sharing
#******************************************************************

########################
# STATISTICAL ANALYSIS
########################

mother_infant_sharing_delivery_place <- na.omit(
  inStrain_results_pairs %>%
    group_by(Mother_infant_comparison_ID, Delivery_place) %>%
    summarise(
      Mother_infant_sharing = ifelse(any(Mother_infant_sharing == "Yes"), "Yes", "No"),
      .groups = "drop"
    )
)

contingency_table_delivery_place <- table(
  mother_infant_sharing_delivery_place$Delivery_place,
  mother_infant_sharing_delivery_place$Mother_infant_sharing
)

chi2_test <- chisq.test(contingency_table_delivery_place)
chi2_test
chi2_test$p.value


########
# PLOT
########

contingency_table_delivery_place_plot <- as.data.frame(as.table(contingency_table_delivery_place))
colnames(contingency_table_delivery_place_plot) <- c("Delivery_place", "Mother_infant_sharing", "Count")

contingency_table_delivery_place_plot <- contingency_table_delivery_place_plot %>%
  group_by(Delivery_place) %>%
  mutate(
    Proportion = Count / sum(Count),
    Total = sum(Count)
  )

pdf('10_STRAIN_TRANSMISSION/Plots/Mother_infant_sharing_delivery_place.pdf', width = 2.5, height = 3.6)
ggplot(contingency_table_delivery_place_plot,
       aes(x = Delivery_place, y = Proportion, fill = Mother_infant_sharing)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_text(aes(label = paste0("n=", Count)),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4.55) +
  scale_fill_manual(values = c("No" = "#E0E0E0", "Yes" = "#A0A0A0")) +
  labs(x = "Infant delivery place", y = "Proportion of sharing events (%)", fill = "Sharing") +
  scale_y_continuous(labels = scales::percent, limits = c(0, 1.15),
                     breaks = c(0, 0.25, 0.5, 0.75, 1),
                     expand = expansion(add = c(0.05, 0.05))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "top",
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14)
  )
dev.off()

#******************************************************************
# 4.C Check how mother-infant strain sharing (in pairs) changes with delivery place
#******************************************************************

#################################################################################
# PLOT: proportion of Mother-infant pairs that have strain sharing per timepoint
#################################################################################

family_sharing_per_timepoint <- family_sharing_per_timepoint %>%
  mutate(
    All_Families_Hospital = sapply(Timepoint_categorical, function(tp) {
      sum(
        Sample_metadata_infants$birth_deliverybirthcard_place_delivery_simple == "hospital" &
          as.character(Sample_metadata_infants$Timepoint_categorical) == as.character(tp) &
          Sample_metadata_infants$FAMILY %in% infants_per_timepoint$FAMILY[
            infants_per_timepoint$Timepoint_categorical == tp
          ],
        na.rm = TRUE
      )
    }),
    All_Families_Home = sapply(Timepoint_categorical, function(tp) {
      sum(
        Sample_metadata_infants$birth_deliverybirthcard_place_delivery_simple == "home" &
          as.character(Sample_metadata_infants$Timepoint_categorical) == as.character(tp) &
          Sample_metadata_infants$FAMILY %in% infants_per_timepoint$FAMILY[
            infants_per_timepoint$Timepoint_categorical == tp
          ],
        na.rm = TRUE
      )
    })
  )

transmission_per_timepoint_hospital <- inStrain_results_pairs %>%
  filter(Mother_infant_sharing == "Yes" & Delivery_place == "hospital") %>%
  group_by(Infant_timepoint) %>%
  summarize(Transmission_Families_Hospital = n_distinct(Pair_ID))

transmission_per_timepoint_home <- inStrain_results_pairs %>%
  filter(Mother_infant_sharing == "Yes" & Delivery_place == "home") %>%
  group_by(Infant_timepoint) %>%
  summarize(Transmission_Families_Home = n_distinct(Pair_ID))

family_sharing_per_timepoint <- family_sharing_per_timepoint %>%
  left_join(transmission_per_timepoint_hospital,
            by = c("Timepoint_categorical" = "Infant_timepoint")) %>%
  left_join(transmission_per_timepoint_home,
            by = c("Timepoint_categorical" = "Infant_timepoint")) %>%
  mutate(Proportion_Transmission_Hospital = Transmission_Families_Hospital / All_Families_Hospital) %>%
  mutate(Proportion_Transmission_Home = Transmission_Families_Home / All_Families_Home)

family_sharing_plot <- family_sharing_per_timepoint %>%
  pivot_longer(
    cols = starts_with("Proportion_Transmission_"),
    names_to = "Delivery_place",
    values_to = "Proportion"
  ) %>%
  mutate(
    Delivery_place = recode(
      Delivery_place,
      "Proportion_Transmission_Hospital" = "hospital",
      "Proportion_Transmission_Home" = "home"
    ),
    Delivery_place = factor(Delivery_place, levels = c("hospital", "home"))
  )

pdf('10_STRAIN_TRANSMISSION/Plots/Proportion_families_with_viral_strain_sharing_per_timepoint_delivery_place.pdf',
    width = 3.6, height = 3.2)
ggplot(family_sharing_plot %>% dplyr::filter(!is.na(Delivery_place)),
  aes(x = Timepoint_categorical, y = Proportion, fill = Delivery_place)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.9) +
  scale_fill_manual(values = c("hospital" = "#A2C2E5", "home" = "#B6D7A8")) +
  labs(x = "Timepoint", y = "Families with sharing (%)", fill = "Delivery place") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0, 0.7)) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    legend.position = "top"
  )
dev.off()

###################
# STATISTICAL TEST
###################

full_family_sharing <- full_family_sharing %>%
  left_join(
    Sample_metadata_infants %>%
      select(FAMILY, birth_deliverybirthcard_place_delivery_simple) %>%
      distinct(),
    by = "FAMILY"
  )

full_family_sharing$birth_deliverybirthcard_place_delivery_simple <- factor(
  full_family_sharing$birth_deliverybirthcard_place_delivery_simple,
  levels = c("hospital", "home")
)

# Change variable name
full_family_sharing <- full_family_sharing %>%
  dplyr::rename(Delivery_place = birth_deliverybirthcard_place_delivery_simple)

# Run GLMM correcting for delivery mode
GLMM_family_sharing_DP <- glmer(
  Any_Sharing ~ Delivery_place + delivery_mode + (1 | FAMILY),
  data = full_family_sharing,
  family = binomial
)

summary(GLMM_family_sharing_DP)


vOTU_sharing_per_family$Delivery_place <- factor(
  vOTU_sharing_per_family$Delivery_place,
  levels = c("hospital", "home")
)

GLMM_vOTU_sharing_DP <- glmer(
  cbind(Num_vOTUs_shared, Num_vOTUs_present_infant - Num_vOTUs_shared) ~
    Delivery_place + Delivery_mode + (1 | FAM_ID),
  data = vOTU_sharing_per_family,
  family = binomial
)

summary(GLMM_vOTU_sharing_DP)
