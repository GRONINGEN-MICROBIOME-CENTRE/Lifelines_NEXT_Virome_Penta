################################################################################
##### LL-NEXT: Viral strain transmission in breastmilk - inStrain
### Author(s): Asier Fernández-Pato
### Last updated: 8th December, 2025
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
Sample_metadata_BM <- read.delim("BREASTMILK/Metadata_BM_VG_Gut_20_05_2025.txt")
Sample_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_07052025.txt")
Sample_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_phenos_07052025.txt")

# Read abundance tables
Abundance_table <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_07052025.txt")
Abundance_table_BM <- read.delim("BREASTMILK/BM_LLNEXT_Viral_Abundance_Table_29072025.txt")
Abundance_table_infants <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_INFANTS_07052025.txt")
Abundance_table_mothers <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_MOTHERS_07052025.txt")

# Read viral metadata tables
Viral_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")
Viral_metadata_BM <- read.delim("BREASTMILK/BM_LLNEXT_Viral_Metadata_29072025.txt")
Viral_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_INFANTS_07052025.txt")
Viral_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_MOTHERS_07052025.txt")

##************************************************************************
# 2. Process in inStrain results: Add metadata 
#*************************************************************************

# After mapping the vOTUs_inStrain genomes to the LL-NEXT samples and running inStrain, we load the results
inStrain_results_BM <- fread("10_STRAIN_TRANSMISSION/all_inStrain_comparisons_BM.tsv")

# Filter those comparisons with a minimum of 75% of the genome compared
inStrain_results_BM <- inStrain_results_BM[inStrain_results_BM$percent_genome_compared >=0.75,]

# Format the results table
inStrain_results_BM$name1 <- gsub("\\.sorted\\.bam$", "", inStrain_results_BM$name1)
inStrain_results_BM$name2 <- gsub("\\.sorted\\.bam$", "", inStrain_results_BM$name2)
inStrain_results_BM <- inStrain_results_BM[, c("scaffold","name1","name2","popANI")]
colnames(inStrain_results_BM) <- c("Virus_ID", "Sample1", "Sample2", "popANI")

# Add sample metadata
inStrain_results_BM <- inStrain_results_BM %>%
  mutate(
    NEXT_ID1 = Sample_metadata_BM$NEXT_ID[match(Sample1, Sample_metadata_BM$NG_ID)],
    NEXT_ID2 = Sample_metadata_BM$NEXT_ID[match(Sample2, Sample_metadata_BM$NG_ID)],
    FAM_ID1 = Sample_metadata_BM$FAMILY[match(Sample1, Sample_metadata_BM$NG_ID)],
    FAM_ID2 = Sample_metadata_BM$FAMILY[match(Sample2, Sample_metadata_BM$NG_ID)],
    Type1 = Sample_metadata_BM$mother_infant[match(Sample1, Sample_metadata$NG_ID)],
    Type2 = Sample_metadata_BM$mother_infant[match(Sample2, Sample_metadata$NG_ID)],
    Origin1 = ifelse(Sample1 %in% Sample_metadata$NG_ID, "Gut", "Milk"),
    Origin2 = ifelse(Sample2 %in% Sample_metadata$NG_ID, "Gut", "Milk"),
    Timepoint1 = Sample_metadata_BM$Timepoint_categorical[match(Sample1, Sample_metadata_BM$NG_ID)],
    Timepoint2 = Sample_metadata_BM$Timepoint_categorical[match(Sample2, Sample_metadata_BM$NG_ID)],
    Timepoint_days1 = Sample_metadata$exact_age_days_at_collection[match(Sample1, Sample_metadata$NG_ID)],
    Timepoint_months1 = Sample_metadata$exact_age_months_at_collection[match(Sample1, Sample_metadata$NG_ID)],
  ) %>%
  mutate(
    Type1 = ifelse(is.na(Type1), "mother", Type1),
    Type2 = ifelse(is.na(Type2), "mother", Type2)
  )


# Determine comparison type and pairs
inStrain_results_BM <- inStrain_results_BM %>%
  mutate(
    Comparison = case_when(
      Origin1 == "Gut" & Origin2 == "Gut" ~ "Gut-Gut",
      Origin1 == "Milk" & Origin2 == "Milk" ~ "Milk-Milk",
      TRUE ~ "Milk-Gut"
    ),
    Mother_Infant_pair = case_when(
      Comparison == "Milk-Gut" & FAM_ID1 == FAM_ID2 ~ "Pair",
      Comparison == "Milk-Gut" & FAM_ID1 != FAM_ID2 ~ "Not pair",
      TRUE ~ NA_character_
    )
  )

# Subset only comparisons involving breastmilk samples
inStrain_results_BM <- inStrain_results_BM[inStrain_results_BM$Comparison == "Milk-Gut",]

# Create individual IDs
inStrain_results_BM <- inStrain_results_BM %>%
  mutate(
    Individual_ID1 = paste(FAM_ID1, NEXT_ID1, tolower(Type1), sep = "_"),
    Individual_ID2 = paste(FAM_ID2, NEXT_ID2, tolower(Type2), sep = "_"),
    Pair_ID = paste(pmin(Individual_ID1, Individual_ID2),
                    pmax(Individual_ID1, Individual_ID2),
                    sep = "_")
  )

# Mother-infant comparison ID
inStrain_results_BM <- inStrain_results_BM %>%
  mutate(
    Mother_infant_comparison_ID = ifelse(
      Comparison == "Milk-Gut",
      ifelse(grepl("mother", Individual_ID1),
             paste(Virus_ID, Individual_ID1, Individual_ID2, sep = "_"),
             paste(Virus_ID, Individual_ID2, Individual_ID1, sep = "_")),
      NA
    )
  )

# Add mother/infant-specific metadata
inStrain_results_BM <- inStrain_results_BM %>%
  mutate(
    infant_sample = ifelse(Type1 == "infant", Sample1, Sample2),
    Mother_sample = ifelse(Type1 == "infant", Sample2, Sample1),
    infant_timepoint = ifelse(Type1 == "infant", Timepoint1, Timepoint2),
    Mother_timepoint = ifelse(Type1 == "infant", Timepoint2, Timepoint1),
    infant_timepoint_days = Timepoint_days1,
    infant_timepoint_months = Timepoint_months1,
    infant_NEXT_ID = ifelse(Type1 == "infant", NEXT_ID1, NEXT_ID2),
    Mother_NEXT_ID = ifelse(Type1 == "infant", NEXT_ID2, NEXT_ID1),
    infant_FAM_ID = ifelse(Type1 == "infant", FAM_ID1, FAM_ID2),
    Mother_FAM_ID = ifelse(Type1 == "infant", FAM_ID2, FAM_ID1),
    Individual_ID_infant = paste(infant_FAM_ID, infant_NEXT_ID, "infant", sep = "_"),
    Individual_ID_Mother = paste(Mother_FAM_ID, Mother_NEXT_ID, "mother", sep = "_")
  ) %>%
  select(-Sample1, -Sample2, -Type1, -Type2,
         -NEXT_ID1, -NEXT_ID2, -FAM_ID1, -FAM_ID2,
         -Individual_ID1, -Individual_ID2,
         -Timepoint1, -Timepoint2, -Timepoint_days1)

# Add information about potential transmission (based on inStrain viral cut-off)
inStrain_results_BM$Mother_infant_sharing <- ifelse(inStrain_results_BM$popANI > 0.99999 , "Yes", "No")

# Add bacterial host information
inStrain_results_BM$Host <- ifelse(
  is.na(Viral_metadata_BM$Bacterial_genus_host[match(inStrain_results_BM$Virus_ID, Viral_metadata_BM$Virus_ID)]),
  Viral_metadata_BM$Bacterial_family_host[match(inStrain_results_BM$Virus_ID, Viral_metadata_BM$Virus_ID)],
  Viral_metadata_BM$Bacterial_genus_host[match(inStrain_results_BM$Virus_ID, Viral_metadata_BM$Virus_ID)]
)
  
# Finally, add infant-specific phenotypes
inStrain_results_BM <- inStrain_results_BM %>%
  mutate(
    Feeding_mode = Sample_metadata_infants$infant_birthcard_feeding_mode_after_birth[match(infant_sample, Sample_metadata_infants$NG_ID)], #after birth
    Delivery_mode = Sample_metadata_infants$delivery_mode[match(infant_sample, Sample_metadata_infants$NG_ID)],
    Infant_eczema = Sample_metadata_infants$infant_health_eczema_diagnosis_strict[match(infant_sample, Sample_metadata_infants$NG_ID)],
    Infant_food_allergy = Sample_metadata_infants$infant_health_food_allergy[match(infant_sample, Sample_metadata_infants$NG_ID)]
  )

# Filter results only for mother-infant pairs
inStrain_results_BM_pair <- inStrain_results_BM[inStrain_results_BM$Mother_Infant_pair == "Pair",]
inStrain_results_BM_pair <- inStrain_results_BM_pair[,c("Virus_ID", "popANI", "Host", "infant_FAM_ID", "infant_timepoint",
                                                  "Mother_timepoint","Mother_infant_sharing")]

##************************************************************************
# 3. General stats. Comparison with feces-feces sharing
#*************************************************************************

#####################
# General stats
#####################

# Count number of vOTUs with strain profiles compared between mother-infants
length(unique(inStrain_results_BM$Virus_ID)) #26

# Count the number of families with transmission
transmission_pairs_BM <- inStrain_results_BM %>%
  filter(Comparison == "Milk-Gut", Mother_infant_sharing == "Yes", Mother_Infant_pair == "Pair")

# Count unique families with transmission 
num_families_transmission_BM <- length(unique(transmission_pairs_BM$Mother_FAM_ID))
num_families_transmission_BM #7

# Select vOTUs shared between mother-infant pairs #13
shared_vOTUs_BM <- unique(inStrain_results_BM$Virus_ID[inStrain_results_BM$Mother_infant_sharing == "Yes" &
                                                      inStrain_results_BM$Mother_Infant_pair == "Pair"])


#####################
# FC-BM Comparison
#####################

# Compare sharing results with the fecal-fecal route
# Load results and subset comparisons within the 91 families of interest
inStrain_results_FC <- read.delim("10_STRAIN_TRANSMISSION/inStrain_results_processed.txt")
milk_families <- Sample_metadata_BM %>%
  filter(Type == "Milk") %>%
  pull(FAMILY)
inStrain_results_FC <- inStrain_results_FC %>%
  filter(Mother_FAM_ID %in% milk_families)

# Count unique families with transmission 
transmission_pairs_FC <- inStrain_results_FC %>%
  filter(Mother_infant_sharing == "Yes", Mother_Infant_pair == "Pair")

num_families_transmission_FC <- length(unique(transmission_pairs_FC$Mother_FAM_ID))
num_families_transmission_FC #51

# Select vOTUs shared between mother-infant pairs #142
shared_vOTUs_FC <- unique(inStrain_results_FC$Virus_ID[inStrain_results_FC$Mother_infant_sharing == "Yes" &
                                                         inStrain_results_FC$Mother_Infant_pair == "Pair"])

# Statistical test: comparison of sharing frequency BM vs FC
# As paired data is used, alternative to Fisher test --> McNemar’s test
pairs_sharing_BM <- unique(transmission_pairs_BM$Mother_FAM_ID)
pairs_sharing_FC <- unique(transmission_pairs_FC$Mother_FAM_ID)
all_pairs <- Sample_metadata_BM %>%
  filter(Type == "Milk") %>%
  pull(FAMILY) %>%
  unique()

# Generate the contingency table
both       <- length(intersect(pairs_sharing_BM, pairs_sharing_FC))           # feces=yes, milk=yes
feces_only <- length(setdiff(pairs_sharing_FC, pairs_sharing_BM))             # feces=yes, milk=no
milk_only  <- length(setdiff(pairs_sharing_BM, pairs_sharing_FC))             # feces=no,   milk=yes
neither    <- length(setdiff(all_pairs, union(pairs_sharing_BM, pairs_sharing_FC))) # feces=no, milk=no

BM_FC_sharing_table <- matrix(
  c(both, feces_only,
    milk_only, neither),
  nrow = 2, byrow = TRUE,
  dimnames = list(Feces = c("Yes","No"), Milk = c("Yes","No"))
)

# Perform the test
mcnemar.test(BM_FC_sharing_table, correct = TRUE) #p-value = 5.42e-10

# Check overlap of shared vOTUs between FC and BM
intersect(shared_vOTUs_BM, shared_vOTUs_FC) #0

######################################################################################
# PLOT: Proportion of pairs (n=91) showing fecal-fecal or milk-feces sharing
######################################################################################

BM_FC_plot <- data.frame(
  SharingType = rep(c("Feces sharing","Milk sharing"), each = 2),
  Sharing = rep(c("Yes","No"), times = 2),
  Count = c( both + feces_only, milk_only + neither, both + milk_only, feces_only + neither)
)

# Add Proportion and reorder the levels
BM_FC_plot$Proportion <- ave(BM_FC_plot$Count, BM_FC_plot$SharingType,
                             FUN = function(x) x / sum(x))
BM_FC_plot$SharingType <- factor(BM_FC_plot$SharingType,
                                 levels = c("Feces sharing","Milk sharing"))
BM_FC_plot$Sharing <- factor(BM_FC_plot$Sharing,
                             levels = c("No","Yes"))  

# Plot
pdf('BREASTMILK/Mother_infant_sharing_proportion_BM_FC.pdf', width = 2.5, height = 3.6)
ggplot(BM_FC_plot, aes(x = SharingType, y = Proportion, fill = Sharing)) +
  geom_bar(stat = "identity", width = 0.9 ) +
  geom_text(aes(label = paste0("n=", Count)),
            position = position_stack(vjust = 0.5),
            color = "black", size = 4.2) +
  scale_fill_manual(values = c("Yes" = "#A0A0A0", "No" = "#E0E0E0"),
                    name = "Sharing") +
  labs(x = "Comparison", y = "Proportion of families (%)") +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(0, 1.1),
                     breaks=c(0, 0.25, 0.50, 0.75,1),
                     expand = expansion(add = c(0.05, 0.08))) +
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

# Finally, compare strain sharing rate in BM vs FC samples
# For this, load abundance table (estimation of the 31,209 vOTUs - including BM viruses - on BM and FC samples from the 91 pairs)
Abundance_table_91 <- read.delim("BREASTMILK/91_BM_FC_Abundance_Table_RPKM.txt")

# Generate inStrain results df only for pairs
inStrain_results_BM_pairs <-inStrain_results_BM[inStrain_results_BM$Mother_Infant_pair == "Pair",]
inStrain_results_FC_pairs <-inStrain_results_FC[inStrain_results_FC$Mother_Infant_pair == "Pair",]

# Count the number of pairs in which strains could be profiled and compared
length(unique(inStrain_results_BM_pairs$Pair_ID)) # 8
length(unique(inStrain_results_FC_pairs$Pair_ID)) # 65

# Generate a dataframe with the number of vOTUs shared per pair (Pair_ID) and per timepoint (related pairs)
vOTU_sharing_per_family_BM <- inStrain_results_BM_pairs %>%
  mutate(FAM_ID = str_extract(Pair_ID, "^[^_]+")) %>% 
  group_by(infant_timepoint, Pair_ID) %>%  
  summarize(
    Num_vOTUs = n_distinct(Virus_ID),
    Num_vOTUs_shared = n_distinct(Virus_ID[Mother_infant_sharing == "Yes"]),  # Unique Virus_IDs where sharing is "Yes"
    .groups = "drop" 
  ) 
vOTU_sharing_per_family_FC <- inStrain_results_FC_pairs %>%
  mutate(FAM_ID = str_extract(Pair_ID, "^[^_]+")) %>% 
  group_by(Infant_timepoint, Pair_ID) %>%  
  summarize(
    Num_vOTUs = n_distinct(Virus_ID),
    Num_vOTUs_shared = n_distinct(Virus_ID[Mother_infant_sharing == "Yes"]),  # Unique Virus_IDs where sharing is "Yes"
    .groups = "drop" 
  ) 

# Identify the total number of vOTUs present in the infant sample for each pair and timepoint
# Note: For some infants, multiple samples of the same timepoint are available (we consider them together)
vOTU_sharing_per_family_BM <- vOTU_sharing_per_family_BM %>%
  rowwise() %>%
  mutate(
    Num_vOTUs_present_infant = {
      # Extract the timepoint, Pair ID and FAM_ID for the current row
      timepoint <- infant_timepoint
      pair <- Pair_ID
      family <- sub("_.*", "", pair)
      infant_NEXT_ID <- inStrain_results_BM_pairs %>%
        filter(Pair_ID == pair) %>%
        pull(infant_NEXT_ID) %>%
        unique()
      
      # Find the SAMPLE_IDs corresponding to Infant NEXT_ID in Sample_metadata_infants
      # Same infants have more than 1 sample at the same timepoint
      ng_ids <- Sample_metadata_infants %>%
        filter(Timepoint_categorical == timepoint & NEXT_ID == infant_NEXT_ID) %>%
        pull(NG_ID)
      
      # Calculate the total number of vOTUs present for each infant and each timepoint
      total_vOTUs <- length(unique(unlist(sapply(ng_ids, function(ng_id) {
        which(Abundance_table_91[[ng_id]] > 0)
      }))))
      
      # Return the total number of vOTUs present
      total_vOTUs
    }
  )

vOTU_sharing_per_family_FC <- vOTU_sharing_per_family_FC %>%
  rowwise() %>%
  mutate(
    Num_vOTUs_present_infant = {
      # Extract the timepoint, Pair ID and FAM_ID for the current row
      timepoint <- Infant_timepoint
      pair <- Pair_ID
      family <- sub("_.*", "", pair)
      infant_NEXT_ID <- inStrain_results_FC_pairs %>%
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
        which(Abundance_table_91[[ng_id]] > 0)
      }))))
      
      # Return the total number of vOTUs present
      total_vOTUs
    }
  )

# Estimate the proportion of infant vOTUs shared per family
colnames(vOTU_sharing_per_family_BM)[colnames(vOTU_sharing_per_family_BM) == "infant_timepoint"] <- "Infant_timepoint"
vOTU_sharing_per_family_BM$infant_timepoint <- NULL

timepoint_order <- c("W2", "M1", "M2", "M3", "M6", "M9", "M12")
vOTU_sharing_per_family_BM$Proportion_sharing <- vOTU_sharing_per_family_BM$Num_vOTUs_shared /vOTU_sharing_per_family_BM$Num_vOTUs_present_infant
vOTU_sharing_per_family_FC$Proportion_sharing <- vOTU_sharing_per_family_FC$Num_vOTUs_shared /vOTU_sharing_per_family_FC$Num_vOTUs_present_infant

vOTU_sharing_per_family_BM$infant_timepoint <- factor(vOTU_sharing_per_family_BM$Infant_timepoint, levels = timepoint_order)
vOTU_sharing_per_family_FC$Infant_timepoint <- factor(vOTU_sharing_per_family_FC$Infant_timepoint, levels = timepoint_order)


# Get average sharing rate at W2
W2_sharing_rate_BM <- as.numeric(unlist(vOTU_sharing_per_family_BM[vOTU_sharing_per_family_BM$Infant_timepoint == "W2",
                                                             "Proportion_sharing"]))
W2_sharing_rate_FC <- as.numeric(unlist(vOTU_sharing_per_family_FC[vOTU_sharing_per_family_FC$Infant_timepoint == "W2",
                                                                "Proportion_sharing"]))
mean(W2_sharing_rate_BM) * 100 
mean(W2_sharing_rate_FC) * 100 

# Combine BM and FC data
colnames(vOTU_sharing_per_family_BM)[colnames(vOTU_sharing_per_family_BM) == "Infant_timepoint"] <- "Timepoint"
colnames(vOTU_sharing_per_family_FC)[colnames(vOTU_sharing_per_family_FC) == "Infant_timepoint"] <- "Timepoint"
vOTU_sharing_per_family_BM$Sample_Type <- "BM"
vOTU_sharing_per_family_FC$Sample_Type <- "FC"
vOTU_sharing_per_family_BM$infant_timepoint <- NULL
vOTU_sharing_combined <- rbind(vOTU_sharing_per_family_BM, vOTU_sharing_per_family_FC)

# Ensure timepoint order
vOTU_sharing_combined$Timepoint <- factor(vOTU_sharing_combined$Timepoint, levels = timepoint_order)

# Estimate number of pairs per timepoint
n_df <- vOTU_sharing_combined %>%
  group_by(Timepoint, Sample_Type) %>%
  summarise(n_pairs = n(), .groups = "drop")

# Generate the plot
pdf('BREASTMILK/Family_sharing_rate_BM_vs_FC_boxplot_n.pdf', width = 3.2, height = 3.5)
ggplot(vOTU_sharing_combined, aes(x = Timepoint, y = Proportion_sharing, fill = Sample_Type)) +
  stat_summary(aes(group = Sample_Type, color = Sample_Type),fun = mean, geom = "line", size = 1.2,
               position = position_dodge(width = 0.6)) +
  stat_summary(aes(group = Sample_Type, color = Sample_Type),fun = mean, geom = "point", size = 3, shape = 18,
               position = position_dodge(width = 0.6)) +
  geom_jitter(color = "lightgrey", aes(size = Num_vOTUs_shared),alpha = 0.5, shape = 16,
              position = position_jitterdodge(jitter.width = 0.15, dodge.width = 0.6)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, size = 0.6, position = position_dodge(width = 0.6), color = "black") +
  geom_text(data = n_df, aes(x = Timepoint, y = 0, label = n_pairs, group = Sample_Type),
            position = position_dodge(width = 0.6), vjust = -9, size = 4.5, inherit.aes = FALSE) +
  labs(x = "Timepoint", y = "vOTU sharing rate (%)", color = "Sample type", fill = "Sample type",
       size = "Number of vOTUs") +
  scale_size_continuous(range = c(1, 5), breaks = c(1, 5, 10)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     limits = c(-0.005, NA)) + 
  scale_fill_manual(values = c("BM"= "#D8A7B1","FC"="#A1887F")) +
  scale_color_manual(values = c("BM"= "#D8A7B1","FC"="#A1887F")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.title = element_text(size = 16, hjust = 0.5),
    legend.position = "bottom",
    legend.box = "vertical",
    legend.title = element_text(size = 12),
    legend.text = element_text(size = 10)
  )
dev.off()


#******************************************************************
# 4. Save results
#******************************************************************
# Save inStrain results table with metadata
write.table(inStrain_results_BM,"BREASTMILK/inStrain_results_BM_processed.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(inStrain_results_BM_pair,"BREASTMILK/inStrain_results_BM_pairs_processed.txt", sep = "\t", row.names = F, quote = FALSE) 
