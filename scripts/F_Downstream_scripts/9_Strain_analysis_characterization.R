################################################################################
##### LL-NEXT: Viral strain transmission - Initial analysis
### Author(s): Asier Fernández-Pato
### Last updated: 5th December, 2025
################################################################################


#****************
# Load libraries
#****************
library(data.table)
library(dplyr)
library(scales)
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


##************************************************************************
# 2. Process in inStrain results: Add metadata 
#*************************************************************************

# After mapping the vOTUs_inStrain genomes to the LL-NEXT samples and running inStrain, we load the results
inStrain_results1 <- fread("10_STRAIN_TRANSMISSION/all_InStrain_comparisons.tsv")
inStrain_results2 <- fread("10_STRAIN_TRANSMISSION/new_InStrain_comparisons.tsv", header = T)
colnames(inStrain_results1) <- colnames(inStrain_results2)

# Combine both inStrain results files (removing rows from first file that also exist in second file (based on first 3 cols))
cols <- c("scaffold", "name1", "name2")
inStrain_results1_filtered <- inStrain_results1[!inStrain_results2, on = cols]
inStrain_results <- rbind(inStrain_results1_filtered, inStrain_results2)

# Filter those comparisons with a minimum of 75% of the genome compared
inStrain_results <- inStrain_results[inStrain_results$percent_genome_compared >=0.75,]

# Filter comparisons involving potential contaminants: LLNEXT_25294 and LLNEXT_79800
contaminants <- c("LLNEXT_25294", "LLNEXT_79800")
inStrain_results <- inStrain_results[!inStrain_results$scaffold %in% contaminants,]

# Format the results table
inStrain_results$name1 <- gsub("\\.sorted\\.bam$", "", inStrain_results$name1)
inStrain_results$name2 <- gsub("\\.sorted\\.bam$", "", inStrain_results$name2)
inStrain_results <- inStrain_results[, c("scaffold","name1","name2","popANI")]
colnames(inStrain_results) <- c("Virus_ID", "Sample1", "Sample2", "popANI")

# Add sample metadata
inStrain_results <- inStrain_results %>%
  mutate(
    NEXT_ID1 = Sample_metadata$NEXT_ID[match(Sample1, Sample_metadata$NG_ID)],
    NEXT_ID2 = Sample_metadata$NEXT_ID[match(Sample2, Sample_metadata$NG_ID)],
    FAM_ID1 = Sample_metadata$FAMILY[match(Sample1, Sample_metadata$NG_ID)],
    FAM_ID2 = Sample_metadata$FAMILY[match(Sample2, Sample_metadata$NG_ID)],
    Type1 = Sample_metadata$Type[match(Sample1, Sample_metadata$NG_ID)],
    Type2 = Sample_metadata$Type[match(Sample2, Sample_metadata$NG_ID)],
    Timepoint1 = Sample_metadata$Timepoint_categorical[match(Sample1, Sample_metadata$NG_ID)],
    Timepoint2 = Sample_metadata$Timepoint_categorical[match(Sample2, Sample_metadata$NG_ID)],
    Timepoint_days1 = Sample_metadata$exact_age_days_at_collection[match(Sample1, Sample_metadata$NG_ID)],
    Timepoint_days2 = Sample_metadata$exact_age_days_at_collection[match(Sample2, Sample_metadata$NG_ID)],
    Timepoint_months1 = Sample_metadata$exact_age_months_at_collection[match(Sample1, Sample_metadata$NG_ID)],
    Timepoint_months2 = Sample_metadata$exact_age_months_at_collection[match(Sample2, Sample_metadata$NG_ID)],
    Host = Viral_metadata$Bacterial_genus_host[match(Virus_ID, Viral_metadata$Virus_ID)]
  )

# Determine comparison type
inStrain_results <- inStrain_results %>%
  mutate(
    Comparison = case_when(
      Type1 == "Mother" & Type2 == "Mother" ~ "Mother-Mother",
      Type1 == "Infant" & Type2 == "Infant" ~ "Infant-Infant",
      TRUE ~ "Mother-Infant"
    ),
    Mother_Infant_pair = case_when(
      Comparison == "Mother-Infant" & FAM_ID1 == FAM_ID2 ~ "Pair",
      Comparison == "Mother-Infant" & FAM_ID1 != FAM_ID2 ~ "Not pair",
      TRUE ~ NA_character_
    )
  )

# Create individual IDs
# Each infant will have a different individual_ID (twins have a different NEXT_ID)
# Each mother will have a different individual_ID (same mother at different pregnancy have different NEXT_ID (_2))
inStrain_results <- inStrain_results %>%
  mutate(
    Individual_ID1 = paste(FAM_ID1, NEXT_ID1, tolower(Type1), sep = "_"),
    Individual_ID2 = paste(FAM_ID2, NEXT_ID2, tolower(Type2), sep = "_"),
    Pair_ID = paste(pmin(Individual_ID1, Individual_ID2),
                    pmax(Individual_ID1, Individual_ID2),
                    sep = "_")
  )

# Determine within/between comparisons
inStrain_results <- inStrain_results %>%
  mutate(
    Within_Between = case_when(
      Individual_ID1 == Individual_ID2 ~ "Within",
      gsub("_2$", "", NEXT_ID1) == gsub("_2$", "", NEXT_ID2) ~ "Within_diff_pregnancy",
      TRUE ~ "Between"
    )
  )

# Mother-infant comparison ID
inStrain_results <- inStrain_results %>%
  mutate(
    Mother_infant_comparison_ID = ifelse(
      Comparison == "Mother-Infant",
      ifelse(grepl("mother", Individual_ID1),
             paste(Virus_ID, Individual_ID1, Individual_ID2, sep = "_"),
             paste(Virus_ID, Individual_ID2, Individual_ID1, sep = "_")),
      NA
    )
  )

# Add mother/infant-specific metadata
inStrain_results <- inStrain_results %>%
  mutate(
    Infant_sample = ifelse(Type1 == "Infant", Sample1, Sample2),
    Mother_sample = ifelse(Type1 == "Infant", Sample2, Sample1),
    Infant_timepoint = ifelse(Type1 == "Infant", Timepoint1, Timepoint2),
    Mother_timepoint = ifelse(Type1 == "Infant", Timepoint2, Timepoint1),
    Infant_timepoint_days = ifelse(Type1 == "Infant", Timepoint_days1, Timepoint_days2),
    Mother_timepoint_days = ifelse(Type1 == "Infant", Timepoint_days2, Timepoint_days1),
    Infant_timepoint_months = ifelse(Type1 == "Infant", Timepoint_months1, Timepoint_months2),
    Mother_timepoint_months = ifelse(Type1 == "Infant", Timepoint_months2, Timepoint_months1),
    Infant_NEXT_ID = ifelse(Type1 == "Infant", NEXT_ID1, NEXT_ID2),
    Mother_NEXT_ID = ifelse(Type1 == "Infant", NEXT_ID2, NEXT_ID1),
    Infant_FAM_ID = ifelse(Type1 == "Infant", FAM_ID1, FAM_ID2),
    Mother_FAM_ID = ifelse(Type1 == "Infant", FAM_ID2, FAM_ID1),
    Individual_ID_Infant = paste(Infant_FAM_ID, Infant_NEXT_ID, "infant", sep = "_"),
    Individual_ID_Mother = paste(Mother_FAM_ID, Mother_NEXT_ID, "mother", sep = "_")
  ) %>%
  select(-Sample1, -Sample2, -Type1, -Type2,
         -NEXT_ID1, -NEXT_ID2, -FAM_ID1, -FAM_ID2,
         -Individual_ID1, -Individual_ID2,
         -Timepoint1, -Timepoint2, -Timepoint_days1, -Timepoint_days2)

# Add information about potential transmission (based on inStrain viral cut-off)
inStrain_results$Mother_infant_sharing <- ifelse(
  inStrain_results$Comparison == "Mother-Infant" & inStrain_results$popANI > 0.99999 , "Yes", "No")

# Finally, add infant-specific phenotypes
inStrain_results <- inStrain_results %>%
  mutate(
    Feeding_mode = Sample_metadata_infants$infant_ffq_ever_never_breastfed[match(Infant_sample, Sample_metadata_infants$NG_ID)],
    Delivery_mode = Sample_metadata_infants$delivery_mode[match(Infant_sample, Sample_metadata_infants$NG_ID)],
    Delivery_place = Sample_metadata_infants$birth_deliverybirthcard_place_delivery_simple[match(Infant_sample, Sample_metadata_infants$NG_ID)],
    Infant_eczema = Sample_metadata_infants$infant_health_eczema_diagnosis_strict[match(Infant_sample, Sample_metadata_infants$NG_ID)],
    Infant_food_allergy = Sample_metadata_infants$infant_health_food_allergy[match(Infant_sample, Sample_metadata_infants$NG_ID)]
  )

# Select vOTUs shared between mother-infant pairs (715)
shared_vOTUs <- unique(inStrain_results$Virus_ID[inStrain_results$Mother_infant_sharing == "Yes" &
                                                   inStrain_results$Mother_Infant_pair == "Pair"])

# Count number of total pairwise comparisons, number of vOTUs involved and number of mother-infant pairs (related) involved
nrow(inStrain_results) #193682
length(unique(inStrain_results$Virus_ID)) #1997
length(unique(inStrain_results$Pair_ID[inStrain_results$Mother_Infant_pair == "Pair"])) #432
nrow(inStrain_results[inStrain_results$Mother_Infant_pair == "Pair"]) #4590 (comparisons across related pairs)

inStrain_results_pair <- inStrain_results[inStrain_results$Mother_Infant_pair == "Pair",]
inStrain_results_pair <- inStrain_results_pair[,c("Virus_ID", "popANI", "Host", "Infant_FAM_ID", "Infant_timepoint",
                                                  "Mother_timepoint","Mother_infant_sharing")]

##************************************************************************
# 3. Exploring Mother-Infant viral strain distances and sharing
#*************************************************************************

######################################################################
# 3.1. Check if related mother-infant pairs have more similar strains
######################################################################

inStrain_results$Mother_Infant_pair <- factor(inStrain_results$Mother_Infant_pair,
                                                            levels = c("Pair", "Not pair"))

# Perform statistical test (use permutation-based test due to lack of independence (repeated measures within infants)
num_permutations <- 10000
permuted_stats <- numeric(num_permutations)
set.seed(123) 

# Extract the observed test statistic: difference in medians between pairs and non-pairs
observed_stat <- median(inStrain_results$popANI[inStrain_results$Mother_Infant_pair == "Pair"]) - 
  median(inStrain_results$popANI[inStrain_results$Mother_Infant_pair == "Not pair"])

# Perform the permutation test
for (i in 1:num_permutations) {
  permuted_labels <- sample(inStrain_results$Mother_Infant_pair)
  permuted_stat <- median(inStrain_results$popANI[permuted_labels == "Pair"]) - 
    median(inStrain_results$popANI[permuted_labels == "Not pair"])
  permuted_stats[i] <- permuted_stat
}

# Calculate the p-value: proportion of permuted stats greater than or equal to observed (one-sided hypothesis)
p_value <- mean(permuted_stats >= observed_stat) # p = 0

# Subset Non-pair to 20,000 rows (for plotting)
pair_rows <- inStrain_results[inStrain_results$Mother_Infant_pair == "Pair", ]
nonpair_rows <- inStrain_results[inStrain_results$Mother_Infant_pair == "Not pair", ]
nonpair_sample <- nonpair_rows[sample(nrow(nonpair_rows), min(20000, nrow(nonpair_rows))), ]
inStrain_subset <- rbind(pair_rows, nonpair_sample)
inStrain_subset$Mother_Infant_pair <- factor(inStrain_subset$Mother_Infant_pair,
                                              levels = c("Pair", "Not pair"))

# Generate boxplot of strain similarity in pairs vs non-pairs
pdf("10_STRAIN_TRANSMISSION/Plots/Mother_infant_distances_relatedness.pdf", width = 2.8, height = 3.2)
ggplot(inStrain_subset, aes(x = Mother_Infant_pair, y = popANI, fill = Mother_Infant_pair)) +
  geom_jitter(alpha = 0.2, aes(color = Mother_Infant_pair), size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = "Relatedness", y = "Strain similarity (popANI)") +
  scale_fill_manual(values = c("Pair" = "#800020","Not pair" = "lightgrey")) +
  scale_color_manual(values = c("Pair" = "#800020","Not pair" = "lightgrey")) + 
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


###########################################################################
# 3.2. Check if related mother-infant pairs have more strain-sharing events
###########################################################################

# Report the overall number of pairs in which we detect some strain sharing (365 pairs)
length(unique(inStrain_results$Pair_ID[inStrain_results$Mother_Infant_pair == "Pair" &
                                 inStrain_results$Mother_infant_sharing == "Yes"]))

# Count potential sharing events in related/unrelated pairs
# Each mother-infant transmission event is counted only once 
mother_infant_sharing <- inStrain_results %>%
  group_by(Mother_infant_comparison_ID, Mother_Infant_pair) %>%
  summarise(
    Mother_infant_sharing = ifelse(any(Mother_infant_sharing == "Yes"), "Yes", "No"),
    .groups = "drop"
  ) %>%
  na.omit()

# # Build contingency table and test significance
contingency_table_sharing <- table(mother_infant_sharing$Mother_Infant_pair, mother_infant_sharing$Mother_infant_sharing)

# Perform chi-square test
chi2_test <- chisq.test(contingency_table_sharing)
chi2_test
chi2_test$p.value # p=0

# Transform contingency table into data frame for ggplot
contingency_table_sharing_plot <- as.data.frame(as.table(contingency_table_sharing))
colnames(contingency_table_sharing_plot) <- c("Pairedness", "Mother_infant_sharing", "Count")

# Calculate proportions for stacked bar plot
contingency_table_sharing_plot <- contingency_table_sharing_plot %>%
  group_by(Pairedness) %>%
  mutate(Proportion = Count / sum(Count),
         Total = sum(Count))

# Generate stacked bar plot with proportions and counts
contingency_table_sharing_plot$Pairedness <- factor(contingency_table_sharing_plot$Pairedness,
                                                    levels = c("Pair", "Not pair"))

pdf('10_STRAIN_TRANSMISSION/Plots/Mother_infant_sharing_relatedness.pdf', width = 2.5, height = 3.6)
ggplot(contingency_table_sharing_plot, aes(x = Pairedness, y = Proportion, fill = Mother_infant_sharing)) +
  geom_bar(stat = "identity", width = 0.9) +
  geom_text(aes(label = paste0("n=", Count)), position = position_stack(vjust = 0.5), color = "black",size = 4.55) +
  scale_fill_manual(values = c("No" = "#E0E0E0", "Yes" = "#A0A0A0")) +
  labs(x = "Relatedness", y = "Proportion of pairs (%)", fill = "Sharing") +
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
# 4. Save results
#******************************************************************
# Save inStrain results table with metadata
write.table(inStrain_results,"10_STRAIN_TRANSMISSION/inStrain_results_processed.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(inStrain_results_pair,"10_STRAIN_TRANSMISSION/inStrain_results_pairs_processed.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(shared_vOTUs,"10_STRAIN_TRANSMISSION/shared_vOTUs.txt", sep = "\t", row.names = F, quote = FALSE) 

