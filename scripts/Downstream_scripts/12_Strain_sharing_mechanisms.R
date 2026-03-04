################################################################################
##### LL-NEXT: Viral strain transmission - Mechanisms
### Author(s): Asier Fernández-Pato
### Last updated: 8th December, 2025
################################################################################

#****************
# Load libraries
#****************
library(dplyr)
library(tidyr)
library(data.table)
library(lme4)
library(lmerTest)
library(stringr)
library(Biostrings)
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


#*********************************************
# 2. Association vOTU sharing and lifestyle
#*********************************************

# Check if vOTU strain sharing occurs more often among temperate phages
# For this, add sharing to Viral_metadata
Viral_metadata$Sharing <- ifelse(
  Viral_metadata$Virus_ID %in% inStrain_results$Virus_ID[inStrain_results$Mother_infant_sharing == "Yes" &
    inStrain_results$Mother_Infant_pair == "Pair"], 1, 0)

Viral_metadata$Lifestyle <- factor(Viral_metadata$Lifestyle, levels = c("Virulent", "Temperate"))

# Test association using a logistic regression:
min_nonzero <- min(Viral_metadata$Mean_rel_abundance_mothers[Viral_metadata$Mean_rel_abundance_mothers > 0], na.rm = TRUE)
pseudo <- min_nonzero / 2

Lifestyle_sharing_GLM <- glm(
  Sharing ~ Lifestyle + log10(Mean_rel_abundance_mothers + pseudo), 
  data = Viral_metadata, 
  family = binomial
)
summary(Lifestyle_sharing_GLM) #LifestyleTemperate 0.43256    0.10815   4.000 6.35e-05  ***

# Plot results
pdf('10_STRAIN_TRANSMISSION/Plots/Lifestyle_sharing.pdf', width = 3.9, height = 4.1)
ggplot(Viral_metadata[!is.na(Viral_metadata$Lifestyle),],
       aes(x = log10(Mean_rel_abundance_mothers + pseudo),y = Sharing, color = Lifestyle,fill = Lifestyle)) +
  geom_boxplot(aes(group = interaction(Sharing, Lifestyle)), position = position_dodge(width = 0.25),
    width = 0.15, alpha = 0.4, outlier.shape = NA, linewidth = 0.6, color = "black") +
  geom_smooth( method = "glm",method.args = list(family = "binomial"), se = TRUE, linewidth = 1.2, alpha = 0.25) +
  labs(
    x = expression(log[10]("Relative abundance in mothers")),
    y = "Probability of mother-infant sharing",
    color = "Lifestyle",
    fill = "Lifestyle"
  ) +
  scale_color_manual(values = c("Virulent" = "#A3C293", "Temperate" = "#D7A67A")) +
  scale_fill_manual(values = c("Virulent" = "#A3C293", "Temperate" = "#D7A67A")) +
  scale_y_continuous(breaks = c(0, 0.25, 0.5, 0.75, 1)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text.x = element_text(size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 8)),
    axis.title.x = element_text(margin = margin(t = 8)),
    legend.position = "top",
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )
dev.off()


#*********************************************
# 3. Exploring mapping results (to MAGs)
#*********************************************

# Select temperate vOTUs shared between mother-infant pairs (435)
shared_vOTUs <- read.delim("10_STRAIN_TRANSMISSION/shared_vOTUs.txt")
temperate_vOTUs <- Viral_metadata[Viral_metadata$Lifestyle == "Temperate" & !is.na(Viral_metadata$Lifestyle),"Virus_ID"]
shared_temp_vOTUs <- shared_vOTUs$x[shared_vOTUs$x %in% temperate_vOTUs] 
inStrain_results_filt <- inStrain_results[inStrain_results$Virus_ID %in% shared_temp_vOTUs,]

# Select NG_IDs of mothers and infants where strain sharing has been detected (within-pairs)
inStrain_results_pair <- inStrain_results_filt[inStrain_results_filt$Mother_Infant_pair == "Pair",]
infant_samples <- unique(inStrain_results_pair[inStrain_results_pair$Mother_infant_sharing == "Yes", "Infant_sample"])
maternal_samples <- unique(inStrain_results_pair[inStrain_results_pair$Mother_infant_sharing == "Yes", "Mother_sample"])
all_samples <- data.frame(Sample_ID = c(infant_samples, maternal_samples))
write.table(all_samples, "10_STRAIN_TRANSMISSION/samples_with_sharing.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)
write.table(shared_temp_vOTUs, "10_STRAIN_TRANSMISSION/shared_temp_vOTUs.txt", sep = "\t", 
            row.names = FALSE, col.names = FALSE, quote = FALSE)

# We perform the mapping of the 435 shared vOTUs to MAGs reconstructed from all_samples with at least medium-quality
# Load mapping results (minimap)
mapping_results <- fread("10_STRAIN_TRANSMISSION/MECHANISM_TRANSMISSION/phages_vs_MAGs.paf", sep = "\t",header = FALSE, fill = TRUE) 
mapping_results <- mapping_results[,1:12]
colnames(mapping_results) <- c("qname","qlen","qstart","qend",
                               "strand","tname","tlen","tstart",
                               "tend","nmatch","alen","mapq")

# Add query coverage (fraction of phage covered) and percent identity
mapping_results$q_cov <- (mapping_results$qend - mapping_results$qstart) / mapping_results$qlen
mapping_results$perc_ident <- (mapping_results$nmatch / mapping_results$alen) * 100

# Select only matches with >=75% of the phage genome covered and >=95% identity
mapping_results <- subset(mapping_results, q_cov >= 0.75 & perc_ident >= 95)

# Add MAG taxonomy and quality to results table
MAG_taxonomy <-  read.delim("10_STRAIN_TRANSMISSION/GTDB_taxonomy.tsv")
MAG_quality <- read.delim("10_STRAIN_TRANSMISSION/CheckM_results.txt")
mapping_results$MAG <- sub("\\|.*", "", mapping_results$tname)

idx_quality <- match(mapping_results$MAG, MAG_quality$MAG_name)
idx_tax <- match(mapping_results$MAG, MAG_taxonomy$user_genome)

mapping_results$Completeness <- MAG_quality$Completeness[idx_quality]
mapping_results$Contamination <- MAG_quality$Contamination[idx_quality]
mapping_results$Taxonomy <- MAG_taxonomy$classification[idx_tax]
mapping_results$Family  <- sub("f__", "", sapply(strsplit(mapping_results$Taxonomy, ";"), function(x) x[5]))
mapping_results$Genus   <- sub("g__", "", sapply(strsplit(mapping_results$Taxonomy, ";"), function(x) x[6]))
mapping_results$Species <- sub("s__", "", sapply(strsplit(mapping_results$Taxonomy, ";"), function(x) x[7]))

# Add phage predicted host (genus-level), lifestyle and DGR presence
idx_host <- match(mapping_results$qname, Viral_metadata$Virus_ID)
mapping_results$Predicted_host_genus <- Viral_metadata$Bacterial_genus_host[idx_host]
mapping_results$DGR <- Viral_metadata$DGR[idx_host]
mapping_results$Lifestyle <- Viral_metadata$Lifestyle[idx_host] #278/435 vOTUs mapped

# Add also sample type (mother/infant) and infant sample timepoint
mapping_results$NG_ID <- sub("_.*", "", mapping_results$tname)
idx_samples <- match(mapping_results$NG_ID, Sample_metadata$NG_ID)
mapping_results$Type <- Sample_metadata$Type[idx_samples]
mapping_results$Timepoint_categorical <- Sample_metadata$Timepoint_categorical[idx_samples]
mapping_results$FAM <- Sample_metadata$FAMILY[idx_samples]
mapping_results$NEXT_ID <- Sample_metadata$NEXT_ID[idx_samples]

# Among the 278 vOTUs mapping to at least 1 MAG, check MAG genus-level tax concordance with host prediction
# Exclude those without host prediction at genus level
vOTU_genus_summary <- mapping_results %>%
  filter(!is.na(Predicted_host_genus)) %>%   
  group_by(qname, Predicted_host_genus) %>%
  summarise(MAG_genera = paste(unique(Genus), collapse = "; "), .groups = "drop") %>%
  # Check if predicted genus is among the MAG genera
  mutate(genus_concordance = sapply(1:n(), function(i) {
    predicted <- Predicted_host_genus[i]
    mag_genera <- strsplit(MAG_genera[i], "; ")[[1]]
    predicted %in% mag_genera
  }))

table(vOTU_genus_summary$genus_concordance)

# Check how many of the wrong predictions are Bacteroides-Phocaiecola changes
false_cases <- vOTU_genus_summary %>% 
  filter(genus_concordance == FALSE)

false_bacteroides_phocaeicola <- false_cases %>%
  filter(Predicted_host_genus == "Bacteroides" & grepl("Phocaeicola", MAG_genera))

nrow(false_bacteroides_phocaeicola)

# Check how many vOTUs are mapping to multiple different species and genera
vOTU_mapping_summary <- mapping_results %>%
  group_by(qname) %>%
  summarise(
    n_species = n_distinct(Species),
    n_genera  = n_distinct(Genus),
    species_list = paste(unique(Species), collapse = "; "),
    genera_list  = paste(unique(Genus), collapse = "; "),
    .groups = "drop"
  )

table(vOTU_mapping_summary$n_species > 1) #44
table(vOTU_mapping_summary$n_genera > 1) #33


#*********************************************
# 4. Exploring sharing events mediated by bacteria
#*********************************************

########################################
# 4A. Explore number of sharing events
########################################

# Get total number of potential sharing events (inStrain) across the 435 temperate vOTUs
# To avoid counting as "2 different events" those ocurring across pregnancies of the same mother with 1 infant we define 3 new variables
inStrain_results_filt <- inStrain_results_filt %>%
  mutate(
    Pregnancy = ifelse(grepl("_2$", Mother_NEXT_ID), "Second", "First"),
    Mother_NEXT_ID_simple = sub("_2$", "", Mother_NEXT_ID),
    Mother_infant_comparison_ID_simple = sub("(_[A-Za-z0-9]+)_2(_mother)", "\\1\\2", Mother_infant_comparison_ID)
  )

inStrain_results_filt_pair_sharing <- inStrain_results_filt[inStrain_results_filt$Mother_Infant_pair == "Pair" &
                                                              inStrain_results_filt$Mother_infant_sharing == "Yes",]


length(unique(inStrain_results_filt_pair_sharing$Mother_infant_comparison_ID_simple)) #670 sharing events
inStrain_potential_sharing_events <- unique(inStrain_results_filt_pair_sharing$Mother_infant_comparison_ID_simple)

# Generate a final table of the mapping results with one row per vOTU-Mother-infant comparison
# Again, we avoid counting twice samples from the same mother-infant pair but different pregnancy
mother_infant_pairs_mapping <- mapping_results %>%
  filter(Type %in% c("Mother", "Infant")) %>%
  group_by(qname, FAM) %>%
  summarise(
    mothers  = list(unique(NEXT_ID[Type == "Mother"])),
    infants  = list(unique(NEXT_ID[Type == "Infant"])),
    .groups = "drop"
  ) %>%
  rowwise() %>%
  mutate(combo = list(expand.grid(mother_next_id = mothers,
                                  infant_next_id = infants,
                                  stringsAsFactors = FALSE))) %>%
  select(-mothers, -infants) %>%
  unnest(combo) %>%
  mutate(
    # Define a simplified Mother ID (remove "_2")
    Mother_NEXT_ID_simple = sub("_2$", "", mother_next_id),
    Mother_infant_comparison_ID_simple = paste0(
      qname, "_", FAM, "_",
      Mother_NEXT_ID_simple, "_mother_", FAM, "_",
      infant_next_id, "_infant"
    )
  ) %>%
  # Remove duplicate rows based on Mother_infant_comparison_ID
  distinct(Mother_infant_comparison_ID_simple, .keep_all = TRUE)

MAG_mediated_potential_sharing_events <- unique(mother_infant_pairs_mapping$Mother_infant_comparison_ID_simple) # 194 in total

# Select which of the potential sharing events (670) were detected to map MAGs of mother-infant pairs (n=136)
inStrain_MAG_mediated_potential_sharing_event_IDs <- intersect(inStrain_potential_sharing_events,
                                                            MAG_mediated_potential_sharing_events)
inStrain_MAG_mediated_potential_sharing_events <- mother_infant_pairs_mapping[
  mother_infant_pairs_mapping$Mother_infant_comparison_ID_simple %in% inStrain_MAG_mediated_potential_sharing_event_IDs,]

# Finally, get vOTU-MAG matches from the potential sharing events and check distribution across timepoints
mother_infant_pairs_mapping_all <- mapping_results %>%
  group_by(qname, FAM) %>%
  filter(any(Type == "Mother") & any(Type == "Infant")) %>%
  ungroup()

mother_infant_pairs_mapping_infants <- mother_infant_pairs_mapping_all[mother_infant_pairs_mapping_all$Type == "Infant",]

mother_infant_pairs_mapping_infants_filtered <- mother_infant_pairs_mapping_infants %>% #238
  semi_join(
    inStrain_MAG_mediated_potential_sharing_events %>% 
      select(qname, infant_next_id),
    by = c("qname" = "qname", "NEXT_ID" = "infant_next_id")
  )

###########
# Plot
###########

# Generate a barplot with the distribution of vOTU matches to MAGs in mother-infant pairs by infant timepoint
sharing_counts <- mother_infant_pairs_mapping_infants_filtered %>%
  group_by(Timepoint_categorical) %>%
  summarise(count = n())

# Get the number of MAGs included per infant timepoint
list_MAGs <- read.delim("10_STRAIN_TRANSMISSION/MECHANISM_TRANSMISSION/list_MAGs_analysis.txt", header = F,
                        col.names = "NG_ID")
list_MAGs <- list_MAGs %>%
  left_join(Sample_metadata %>% select(NG_ID, Timepoint_categorical, Type), by = "NG_ID")

list_MAGs_infants <- list_MAGs %>%
  filter(Type == "Infant") %>%
  group_by(Timepoint_categorical) %>%
  summarise(n_MAGs = n()) %>%
  ungroup()

# Generate the plot
timepoint_colors <- c("#4B4FC5", "#3B81B3", "#4CA6B1", "#88CFA4", "#C1D97F", "#E8C26D", "#F4A219")
sharing_counts$Timepoint_categorical <- factor(sharing_counts$Timepoint_categorical, 
  levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12"))

pdf('10_STRAIN_TRANSMISSION/Plots/Number_MAG_vOTU_matches_in_pairs.pdf', width = 4.1, height = 3.9)
ggplot(sharing_counts, aes(x = Timepoint_categorical, y = count, fill = Timepoint_categorical)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8,
              position = position_jitter(width = 0.2, height = 0)) +
  geom_bar(stat = "identity", width = 0.9, alpha = 0.7) +
  scale_fill_manual(values = timepoint_colors) +
  labs(x = "Infant timepoint",y = "Number of vOTU-MAG matches in pairs (across sharing events)") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "none")
dev.off()


########################################
# 4B. Enrichment of mapping vs random pairs
########################################

# Goal: Test if vOTU detection in MAGs from related mother-infant pairs are more common than in unrelated pairs
mothers <- unique(mapping_results$NEXT_ID[mapping_results$Type == "Mother"])
infants <- unique(mapping_results$NEXT_ID[mapping_results$Type == "Infant"])

# Observed number of mappings (n=194) (in total)
observed_shared <-  nrow(mother_infant_pairs_mapping)

#  Run permutation test
n_permutations <- 10000
random_counts <- numeric(n_permutations)

for (i in seq_len(n_permutations)) {
  
  # Suffle the infant IDs (randomly)
  mapping_results_random <- mapping_results
  shuffled_infants <- sample(infants, length(infants), replace = FALSE)
  mapping_results_random$FAM[mapping_results_random$Type == "Infant"] <- 
    sample(unique(mapping_results_random$FAM[mapping_results_random$Type == "Infant"]))  # random FAMs
  
  # Estimate "mother–infant pair" mappings
  random_pairs <- mapping_results_random %>%
    group_by(qname, FAM) %>%
    filter(any(Type == "Mother") & any(Type == "Infant")) %>%
    ungroup()
  
  random_counts[i] <- nrow(random_pairs)
}

# Estimate permutation p-value
p_value <- mean(random_counts >= observed_shared)
observed_shared
mean(random_counts)

# Generate the plot to visualize results of permutation test
random_df <- data.frame(Random_counts = random_counts)

pdf('10_STRAIN_TRANSMISSION/Plots/vOTU_MAG_mapping_permutation.pdf', width = 4.1, height = 3.9)
ggplot(random_df, aes(x = Random_counts)) +
  geom_histogram(binwidth = 10, fill = "#C1D97F", color = "black", alpha = 0.7) +
  geom_vline(xintercept = 194, color = "#4B4FC5", size = 1.2, linetype = "dashed") +
  labs(x = "Number of vOTU–MAG mapping events",
       y = "Frequency") +
  annotate("text", x = 194, y = max(table(cut(random_counts, breaks = 20))) * 0.9,
           label = "Observed = 194", hjust = -0.1, color = "#4B4FC5", size = 4.5) +
  coord_cartesian(ylim = c(0, 4000), xlim= c(0, 220)) +
  theme_bw(base_size = 14) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14))
dev.off()

########################################
# 4C. MAG selection for ANI comparison
########################################

# To confirm the strain-level similarity of MAGs (>=99.9% ANI), we use SKANI
# We compare only those MAGs from pairs that are mapping to the same vOTU (from the 136 potential sharing events) (n=383)
# First, we select those MAGs to be compared
sharing_events <- inStrain_MAG_mediated_potential_sharing_event_IDs

mag_results_list <- list()
for(event in sharing_events) {
  
  # Extract vOTU and FAMILY
  vOTU <- str_extract(event, ".*?(?=_FAM)")
  fam <- str_extract(event, "FAM[0-9]+")
  
  # Extract MAGs for the mother
  mother_mags <- mother_infant_pairs_mapping_all %>%
    filter(qname == vOTU, FAM == fam, Type == "Mother") %>%
    pull(MAG) %>% unique()
  
  # Extract MAGs for the infant
  infant_mags <- mother_infant_pairs_mapping_all %>%
    filter(qname == vOTU, FAM == fam, Type == "Infant") %>%
    pull(MAG) %>% unique()
  
  # Combine into a single tibble for this event
  mag_results_list[[event]] <- tibble(
    sharing_event = event,
    mother_MAGs = paste(mother_mags, collapse = "; "),
    infant_MAGs = paste(infant_mags, collapse = "; ")
  )
}

# Combine all results into one dataframe
MAGs_per_sharing_event <- bind_rows(mag_results_list)

# Save all MAGs potentially involved in sharing events into a file
all_MAGs_sharing_events <- MAGs_per_sharing_event %>%
  select(mother_MAGs, infant_MAGs) %>%
  pivot_longer(cols = everything(), values_to = "MAG") %>%
  separate_rows(MAG, sep = "; ") %>%
  distinct(MAG) 

# Save to a file
write_lines(all_MAGs_sharing_events$MAG, "10_STRAIN_TRANSMISSION/MECHANISM_TRANSMISSION/mother_infant_MAGs_for_SKANI.txt")

########################################
# 4D. Confirm bacterial strain-level similarity
########################################

# Read SKANI results
SKANI_results <- read.delim("10_STRAIN_TRANSMISSION/MECHANISM_TRANSMISSION/skani_matrix.txt", header = F)

# Add MAG names as row and col names
mag_names <- basename(SKANI_results[,1])               
mag_names <- sub("\\.fa$", "", mag_names)            
mag_names <- sub("^MAGs_SKANI/", "", mag_names)
SKANI_results_clean <- SKANI_results[, -1]
SKANI_results_clean[] <- lapply(SKANI_results_clean, as.numeric)
rownames(SKANI_results_clean) <- mag_names
colnames(SKANI_results_clean) <- mag_names

# Add to MAGs_per_sharing_event SKANI results + taxonomy of transmitted strains
MAGs_per_sharing_event$Strain_level <- FALSE
MAGs_per_sharing_event$Strain_species <- NA_character_

for (i in seq_len(nrow(MAGs_per_sharing_event))) {
  
  mothers <- unlist(strsplit(MAGs_per_sharing_event$mother_MAGs[i], ";\\s*"))
  infants <- unlist(strsplit(MAGs_per_sharing_event$infant_MAGs[i], ";\\s*"))
  
  # Skip if no valid MAGs
  if (length(mothers) == 0 || length(infants) == 0 ||
      all(mothers == "") || all(infants == "")) next
  
  valid_mothers <- mothers[mothers %in% rownames(SKANI_results_clean)]
  valid_infants <- infants[infants %in% colnames(SKANI_results_clean)]
  
  if (length(valid_mothers) == 0 || length(valid_infants) == 0) next
  
  # Extract all relevant ANI values
  ani_vals <- SKANI_results_clean[valid_mothers, valid_infants, drop = FALSE]
  
  # Identify infant MAGs with any mother ANI >= 99.9
  infant_strain_hits <- colnames(ani_vals)[
    apply(ani_vals, 2, function(x) any(x >= 99.9, na.rm = TRUE))
  ]
  
  if (length(infant_strain_hits) > 0) {
    MAGs_per_sharing_event$Strain_level[i] <- TRUE
    
    # Retrieve species taxonomy for those infant MAGs
    shared_species <- mother_infant_pairs_mapping_all %>%
      filter(MAG %in% infant_strain_hits) %>%
      pull(Species) %>%
      unique()
    
    MAGs_per_sharing_event$Strain_species[i] <- paste(shared_species, collapse = "; ")
  }
}

# Add vOTU and FAM info
MAGs_per_sharing_event <- MAGs_per_sharing_event %>%
  mutate(
    vOTU = sub("_FAM.*", "", sharing_event),
    FAM = str_extract(sharing_event, "FAM\\d+")
  )

table(MAGs_per_sharing_event$Strain_level) # 111 events mediated my bacterial-sharing 

############
# Plot
###########

# List of timepoints
timepoints <- c("M12","M9","M6","M3","M2","M1","W2")

# Build presence/absence matrix for infant MAG detection
presence_matrix <- mother_infant_pairs_mapping_all %>%
  filter(Type == "Infant") %>%
  semi_join(MAGs_per_sharing_event, by = c("qname" = "vOTU", "FAM")) %>%
  mutate(present = 1) %>%
  distinct(qname, FAM, Timepoint_categorical, .keep_all = TRUE) %>%
  select(qname, FAM, Timepoint_categorical, present) %>%
  pivot_wider(
    names_from = Timepoint_categorical,
    values_from = present,
    values_fill = 0
  )

# Merge with main table and assign labels
heatmap_data <- MAGs_per_sharing_event %>%
  left_join(presence_matrix, by = c("vOTU" = "qname", "FAM")) %>%
  mutate(
    Label = ifelse(!is.na(Strain_species) & Strain_species != "",
                   Strain_species,
                   sharing_event)
  )

# Sort sharing events: strain-level first, non-strain last
heatmap_data <- heatmap_data %>%
  arrange(desc(Strain_level)) %>%
  mutate(sharing_event = factor(sharing_event, levels = sharing_event))

# Pivot longer for plotting
heatmap_long <- heatmap_data %>%
  pivot_longer(
    cols = all_of(timepoints),
    names_to = "Timepoint",
    values_to = "Detected"
  ) %>%
  mutate(Detected = factor(Detected, levels = c(0, 1)),
         Type = "Detected")

# Add strain-level row
strain_row <- heatmap_data %>%
  select(sharing_event, Label, Strain_level) %>%
  mutate(
    Timepoint = "Strain",
    Detected = factor(ifelse(Strain_level, 1, 0), levels = c(0, 1)),
    Type = "Strain"
  )

# Combine detection + strain rows
heatmap_plot_data <- bind_rows(
  heatmap_long %>% select(sharing_event, Label, Timepoint, Detected, Type),
  strain_row %>% select(sharing_event, Label, Timepoint, Detected, Type)
)

# Plot
pdf("10_STRAIN_TRANSMISSION/Plots/Co-transmission.pdf", width = 12, height = 3)
ggplot(heatmap_plot_data, aes(x = sharing_event, y = Timepoint)) +
  geom_tile(aes(fill = interaction(Type, Detected)), color = "grey80") +
  scale_fill_manual(values = c("Detected.0" = "white","Detected.1" = "#88CFA4",
               "Strain.0"   = "white","Strain.1"   = "#4B4FC5"),name = NULL) +
  geom_text( data = distinct(heatmap_plot_data, sharing_event, Label),
    aes(x = sharing_event, y = 8, label = Label), size = 2.8, angle = 90, hjust = 0) +
  scale_y_discrete(limits = c(timepoints, "Strain")) +
  labs(x = NULL, y = "Infant timepoint") +
  theme_minimal(base_size = 13) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.text.y = element_text(size = 11),
    axis.title.y = element_text(size = 12),
    legend.position = "top",
    panel.grid = element_blank(),
    panel.border = element_blank(),
    plot.margin = margin(10, 20, 10, 10)
  )
dev.off()

########################################
# 4E. Confirm common integration sites (flanking regions)
########################################

# First, filter only sharing events with MAGs shared at strain-level
strain_sharing <- MAGs_per_sharing_event %>%
  filter(Strain_level) %>%
  select(vOTU, FAM, mother_MAGs, infant_MAGs)

# Filter also the mapping results only involving those MAGs
mapping_strain_pairs <- mapping_results %>%
  semi_join(strain_sharing, by = "FAM")

#******************************************************************
# 4. Save results
#******************************************************************
write.table(vOTU_mapping_summary,"10_STRAIN_TRANSMISSION/MECHANISM_TRANSMISSION/vOTU_mapping_summary.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(MAGs_per_sharing_event,"10_STRAIN_TRANSMISSION/MECHANISM_TRANSMISSION/MAG_sharing_per_transmission_event.txt",
            sep = "\t", row.names = F, quote = FALSE) 


