################################################################################
##### LL-NEXT: Phage Anti-Defense Systems - Analysis
### Author(s):Asier Fernández-Pato
### Last updated: 19th December, 2025
################################################################################

#****************
# Load modules
#****************
library(dplyr)
library(ggplot2)
library(tidyr)
library(tidyverse)
library(ggpattern)


#****************
# Load data
#****************
# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_analysis/VIRUSES/")

# Read metadata tables 
Sample_metadata <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt")
Sample_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_07052025.txt")
Sample_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_07052025.txt")

# Read abundance tables
Abundance_table <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_07052025.txt")
Abundance_table_infants <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_INFANTS_07052025.txt")
Abundance_table_mothers <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_MOTHERS_07052025.txt")

# Read viral metadata tables
Viral_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")
Viral_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_INFANTS_07052025.txt")
Viral_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_MOTHERS_07052025.txt")


##************************************************************************
# 1. Anti-Defense Systems: General analysis
##************************************************************************
#Read AntiDefenseFinder results
antidefense_systems <- read.delim("6_ANTIDEFENSE_ANALYSIS/Present_vOTUs_proteins_defense_finder_systems.tsv", 
                             header = T, check.names = F)

antidefense_systems <- distinct(antidefense_systems) #Get unique systems

# Add Virus_ID to each system 
antidefense_systems$Virus_ID <- sub("_(?!.*_).*", "", antidefense_systems$sys_beg, perl = TRUE)


# Add variable with the presence/absence of antidefense system in phages to the Virus Metadata
Viral_metadata$Antidefense_system <- ifelse(Viral_metadata$Virus_ID %in% 
                                              antidefense_systems$Virus_ID, "Yes", "No")

#-------------------------------------------------------
# Add variables for the most common ADS systems (top 5)
#-------------------------------------------------------
system_ids_list <- antidefense_systems %>%
  filter(type %in% c("Anti_CBASS", "Anti_CRISPR", "Anti_RM", 
                     "Anti_RecBCD", "Anti_Thoeris")) %>%
  distinct(Virus_ID, type) %>%
  group_by(type) %>%
  summarise(ids = list(unique(Virus_ID)), .groups = "drop") %>%
  deframe()

for (system in names(system_ids_list)) {
  Viral_metadata[[paste0(system, "_system")]] <- ifelse(
    Viral_metadata$Virus_ID %in% system_ids_list[[system]], "Yes", "No"
  )
}

#-------------------------------------------------------
# Add variables for Anti-RM subtypes with at least 50 occurrences
#-------------------------------------------------------
# For hia5_hin1523_nma1821 anti_RM subtype:
# --keep only those with name_of_profiles_in_sys == hin1523 (vast majority)
filtered_antidefense <- antidefense_systems %>%
  filter(type == "Anti_RM") %>%
  group_by(subtype) %>%
  filter(n() >= 50) %>%
  ungroup()

filtered_antidefense <- filtered_antidefense %>%
  filter(!(subtype == "hia5_hin1523_nma1821")) %>%  # exclude all from this subtype first
  bind_rows(
    antidefense_systems %>%
      filter(type == "Anti_RM", subtype == "hia5_hin1523_nma1821", name_of_profiles_in_sys == "hin1523")
  )

rm_subtype_ids_list <- filtered_antidefense %>%
  distinct(Virus_ID, subtype) %>%
  group_by(subtype) %>%
  summarise(ids = list(unique(Virus_ID)), .groups = "drop") %>%
  deframe()

for (subtype in names(rm_subtype_ids_list)) {
  Viral_metadata[[paste0("Anti_RM_", subtype)]] <- ifelse(
    Viral_metadata$Virus_ID %in% rm_subtype_ids_list[[subtype]], "Yes", "No"
  )
}

Viral_metadata_infants <- Viral_metadata[Viral_metadata$Virus_ID %in% rownames(Abundance_table_infants),]
Viral_metadata_mothers <- Viral_metadata[Viral_metadata$Virus_ID %in% rownames(Abundance_table_mothers),]

#######################################
# Anti-defense systems: Overall
#######################################

# Total count of systems identified in total
length(table(antidefense_systems$sys_id)) # 3,400 total systems identified
length(table(antidefense_systems$Virus_ID)) # 3,090 phages with systems 

# Count the occurrences of each type and subtype of antidefense system (overall)
antidefense_systems_counts <- antidefense_systems %>%
  count(type) %>%
  na.omit()

antidefense_system_subtype_counts <- antidefense_systems %>%
  count(subtype) %>%
  na.omit()

# Select the top 8 more common for the plot
antidefense_systems_counts_top8 <- antidefense_systems_counts %>%
  arrange(desc(n)) %>%
  slice_head(n = 8) 

common_ADS_types <- antidefense_systems_counts_top8$type

pdf('6_ANTIDEFENSE_ANALYSIS/Plots/Antidefense_systems_barplot.pdf', width = 4.5, height = 3.5)
ggplot(antidefense_systems_counts_top8 , aes(x = reorder(type, -n), y = n)) +
  geom_bar_pattern(stat = "identity", width = 0.7, pattern = "stripe", pattern_color = "#707070", 
                   pattern_fill = "#707070", pattern_density = 0.15,
                   pattern_angle = 45, color = "#707070", fill = "#D3D3D3" ) +
  coord_flip() +
  ylab("Number of Antidefense Systems") +
  scale_fill_manual(values = rep("#D3D3D3", length(unique(antidefense_systems_counts$type)))) + 
  theme_minimal(base_size = 16) +  
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 17),
    axis.text.y = element_text(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),                       
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    plot.background = element_blank()
  )
dev.off()


#########################################
# Anti-defense systems: Infants
#########################################

# Get results for infant phages
antidefense_systems_infants <- antidefense_systems[antidefense_systems$Virus_ID %in%
                                                     rownames(Abundance_table_infants),]
# Total count of systems identified
length(table(antidefense_systems_infants$sys_id)) # 2,126 total systems identified
table(Viral_metadata_infants$Antidefense_system) # 1,873 phages with systems 

# Count the ocurrence of anti-defense systems in infant phages
antidefense_systems_counts_infants <- antidefense_systems %>%
  filter(Virus_ID %in% Viral_metadata_infants$Virus_ID) %>%  
  count(type) %>%                                           
  na.omit() 

antidefense_systems_subtytpe_counts_infants <- antidefense_systems %>%
  filter(Virus_ID %in% Viral_metadata_infants$Virus_ID) %>%  
  count(subtype) %>%                                           
  na.omit() 

# Select the top 8 more common for the plot
antidefense_systems_counts_infants_top8 <- antidefense_systems_counts_infants %>%
  arrange(desc(n)) %>%
  slice_head(n = 8) 

# Generate the barplot with the number of anti-defense systems identified (type)
pdf('6_ANTIDEFENSE_ANALYSIS/Plots/Antidefense_systems_barplot_INFANTS.pdf', width = 4.5, height = 3.5)
ggplot(antidefense_systems_counts_infants_top8 , aes(x = reorder(type, -n), y = n)) +
  geom_bar_pattern(stat = "identity", width = 0.7, pattern = "stripe", pattern_color = "#707070", 
                   pattern_fill = "#707070", pattern_density = 0.15, pattern_angle = 45, color = "#707070", fill = "#D3D3D3") +
  coord_flip() +
  ylab("Number of Antidefense Systems") +
  scale_fill_manual(values = rep("#D3D3D3", length(unique(antidefense_systems_counts$type)))) + 
  theme_minimal(base_size = 16) +  
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 17),
    axis.text.y = element_text(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(), 
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    plot.background = element_blank()
  )
dev.off()

# Get a barplot with the proportion of anti-defense systems per timepoint. For this:
# Get the anti-defense systems per Virus_ID
antidefense_systems_per_phage_infants <- antidefense_systems_infants %>%
  group_by(Virus_ID, type) %>%
  summarise(count = n(), .groups = "drop")

# Reshape Abundance_table_infants into a long format
Abundance_table_infants_long <- Abundance_table_infants %>%
  rownames_to_column("Virus_ID") %>%
  filter(Virus_ID %in% antidefense_systems_per_phage_infants$Virus_ID) %>% # Filter by Virus_ID
  pivot_longer(-Virus_ID, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Present = ifelse(Abundance > 0, 1, 0)) 

# Merge abundance data with defense system counts
antidefense_systems_sample_counts_infants <- Abundance_table_infants_long %>%
  left_join(antidefense_systems_per_phage_infants, by = "Virus_ID") 

# Select only those viruses that are present
antidefense_systems_sample_counts_infants <- antidefense_systems_sample_counts_infants[antidefense_systems_sample_counts_infants$Present !=0,]

# Add Sample timepoint information
antidefense_systems_sample_counts_infants <- antidefense_systems_sample_counts_infants %>%
  left_join(Sample_metadata_infants %>%
              select(NG_ID, Timepoint_categorical), 
            by = c("Sample" = "NG_ID"))

# Calculate the counts per system (type) and normalize to get proportions per timepoint
# Filter only ADS across the 8 mot common types
antidefense_systems_sample_counts_infants_proportions <- antidefense_systems_sample_counts_infants %>%
  group_by(Timepoint_categorical, type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Timepoint_categorical) %>%
  mutate(proportion = count / sum(count)) %>%
  ungroup() %>%
  filter(type %in% common_ADS_types) 

# Reorder Timepoint categorical variable
antidefense_systems_sample_counts_infants_proportions$Timepoint_categorical <- factor(
  antidefense_systems_sample_counts_infants_proportions$Timepoint_categorical,
  levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12")
)

# Generate the barplot with the number of anti-defense systems identified (type)
type_colors <- c("#6C7A89", "#95A5A6", "#D2B48C", "#A9A9A9", "#556B2F",
                 "#8FBC8F", "#4682B4", "#D8BFD8", "#C2A5CF", "#B5D8EB",
                 "#8C564B", "#C49C94")
type_colors <- adjustcolor(type_colors, alpha.f = 0.7)
type_names <- unique(antidefense_systems_counts$type)
color_mapping <- setNames(type_colors[1:length(type_names)], type_names)

pdf('6_ANTIDEFENSE_ANALYSIS/Plots/Infant_Defense_System_Proportions_barplot.pdf', width = 6, height = 3.1)
ggplot(antidefense_systems_sample_counts_infants_proportions, aes(x = Timepoint_categorical, y = proportion, fill = type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +  
  scale_fill_manual(values = color_mapping) +
  theme_minimal(base_size = 14) +  
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    panel.grid = element_blank(),          
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),  
    plot.background = element_blank()) +
  ylab("Proportion of Anti-defense Systems (%)") +
  xlab("Timepoint") +
  guides(fill = guide_legend(title = "Anti-defense System"))
dev.off()

#---------------------------
# Analysis per bacterial host
#---------------------------

# Add information about the bacterial host
antidefense_systems_per_phage_infants <- antidefense_systems_per_phage_infants %>%
  left_join(Viral_metadata_infants %>%
              select(Virus_ID, Bacterial_genus_host), 
            by = "Virus_ID")

# Calculate the counts of anti-defense systems per bacterial host
antidefense_systems_per_phage_infants <- antidefense_systems_per_phage_infants %>%
  filter(!is.na(Bacterial_genus_host)) %>%
  group_by(Bacterial_genus_host, type) %>%
  summarise(count = n(), .groups = 'drop')

# Filter to keep only genera with at least 20 defense systems
filtered_bacterial_genus <- antidefense_systems_per_phage_infants %>%
  group_by(Bacterial_genus_host) %>%         
  summarise(total_count = sum(count, na.rm = TRUE)) %>%  
  filter(total_count >= 20) %>%               
  pull(Bacterial_genus_host)      

# Filter the data for defense systems types, keeping only relevant genera
antidefense_systems_per_phage_infants <- antidefense_systems_per_phage_infants %>%
  filter(Bacterial_genus_host %in% filtered_bacterial_genus) 

# Add all combinations (including those with count = 0)
# Select only the top 8 most common ADS types
antidefense_systems_per_phage_infants_complete <- expand_grid(
  Bacterial_genus_host = unique(antidefense_systems_per_phage_infants$Bacterial_genus_host),
  type = unique(antidefense_systems_per_phage_infants$type)
) %>%
  left_join(antidefense_systems_per_phage_infants, by = c("Bacterial_genus_host", "type")) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  filter(type %in% common_ADS_types) 

# Calculate total counts per genus for ordering x-axis
genus_order <- antidefense_systems_per_phage_infants_complete %>%
  group_by(Bacterial_genus_host) %>%
  summarise(total_count = sum(count, na.rm = TRUE)) %>%
  arrange(desc(total_count)) %>%
  pull(Bacterial_genus_host)

# Calculate total counts per type for ordering y-axis
type_order <- antidefense_systems_per_phage_infants_complete %>%
  group_by(type) %>%
  summarise(total_count = sum(count, na.rm = TRUE)) %>%
  arrange(desc(total_count)) %>%
  pull(type)

# Convert columns to factors with the desired order
antidefense_systems_per_phage_infants_complete <- antidefense_systems_per_phage_infants_complete %>%
  mutate(
    Bacterial_genus_host = factor(Bacterial_genus_host, levels = genus_order),
    type = factor(type, levels = type_order)
  )

# Create the heatmap
pdf('6_ANTIDEFENSE_ANALYSIS/Plots/Infant_Antidefense_System_Counts_per_Type_Bacterial_Host.pdf', width = 7, height = 3.5)
ggplot(antidefense_systems_per_phage_infants_complete, aes(x = Bacterial_genus_host, y = type, fill = log10(count + 1))) +  
  geom_tile(color = "lightgrey", linewidth = 0.5, linetype = 1) +
  scale_fill_gradientn(colors = c("white", "lightblue", "maroon"), 
    values = scales::rescale(c(0, 0.5, 1)), 
    limits = c(0, log10(max(antidefense_systems_per_phage_infants_complete$count + 1))),
    name = "Log(Count)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  ) +
  ylab("Anti-defense System") + 
  xlab("Bacterial Genus Host") + 
  guides(fill = guide_colorbar(title = "Log10(Count)"))
dev.off()

#########################################
# Anti-defense systems: Mothers
#########################################

# Get results for maternal phages
antidefense_systems_mothers <- antidefense_systems[antidefense_systems$Virus_ID %in%
                                                     rownames(Abundance_table_mothers),]
# Total count of systems identified
length(table(antidefense_systems_mothers$sys_id)) # 1,967 total systems identified
table(Viral_metadata_mothers$Antidefense_system) # 1,868 phages with systems 

# Count the ocurrence of anti-defense systems in infant phages
antidefense_systems_counts_mothers <- antidefense_systems %>%
  filter(Virus_ID %in% Viral_metadata_mothers$Virus_ID) %>%  
  count(type) %>%                                           
  na.omit() 

antidefense_systems_subtytpe_counts_mothers <- antidefense_systems %>%
  filter(Virus_ID %in% Viral_metadata_mothers$Virus_ID) %>%  
  count(subtype) %>%                                           
  na.omit() 

# Select the top 8 more common for the plot
antidefense_systems_counts_mothers_top8 <- antidefense_systems_counts_mothers %>%
  arrange(desc(n)) %>%
  slice_head(n = 8) 

# Generate the barplot with the number of anti-defense systems identified (type)
pdf('6_ANTIDEFENSE_ANALYSIS/Plots/Antidefense_systems_barplot_MOTHERS.pdf', width = 4.5, height = 3.5)
ggplot(antidefense_systems_counts_mothers_top8, aes(x = reorder(type, -n), y = n)) +
  geom_bar_pattern(stat = "identity", width = 0.7, pattern = "stripe", pattern_color = "#707070", 
                   pattern_fill = "#707070", pattern_density = 0.15,
                   pattern_angle = 45, color = "#707070", fill = "#D3D3D3" ) +
  coord_flip() +
  ylab("Number of Antidefense Systems") +
  scale_fill_manual(values = rep("#D3D3D3", length(unique(antidefense_systems_counts_mothers$type)))) + 
  theme_minimal(base_size = 16) +  
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 17),
    axis.text.y = element_text(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),                        
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),  
    plot.background = element_blank()
  )
dev.off()

# Get a barplot with the proportion of anti-defense systems per timepoint. For this:
# Get the anti-defense systems per Virus_ID
antidefense_systems_per_phage_mothers <- antidefense_systems_mothers %>%
  group_by(Virus_ID, type) %>%
  summarise(count = n(), .groups = "drop")

# Reshape Abundance_table_mothers into a long format
Abundance_table_mothers_long <- Abundance_table_mothers %>%
  rownames_to_column("Virus_ID") %>%
  filter(Virus_ID %in% antidefense_systems_per_phage_mothers$Virus_ID) %>% # Filter by Virus_ID
  pivot_longer(-Virus_ID, names_to = "Sample", values_to = "Abundance") %>%
  mutate(Present = ifelse(Abundance > 0, 1, 0)) 

# Merge abundance data with antidefense system counts
antidefense_systems_sample_counts_mothers <- Abundance_table_mothers_long %>%
  left_join(antidefense_systems_per_phage_mothers, by = "Virus_ID") 

# Select only those viruses that are present
antidefense_systems_sample_counts_mothers <- antidefense_systems_sample_counts_mothers[antidefense_systems_sample_counts_mothers$Present !=0,]

# Add Sample timepoint information
antidefense_systems_sample_counts_mothers <- antidefense_systems_sample_counts_mothers %>%
  left_join(Sample_metadata_mothers %>%
              select(NG_ID, Timepoint_categorical), 
            by = c("Sample" = "NG_ID"))

# Calculate the counts per system (type) and normalize to get proportions per timepoint
# Filter also only the 8 most common ADS
antidefense_systems_sample_counts_mothers_proportions <- antidefense_systems_sample_counts_mothers %>%
  group_by(Timepoint_categorical, type) %>%
  summarise(count = n(), .groups = 'drop') %>%
  group_by(Timepoint_categorical) %>%
  mutate(proportion = count / sum(count)) %>%
  filter(type %in% common_ADS_types) 

# Generate the plot

# Reorder Timepoint categorical variable
antidefense_systems_sample_counts_mothers_proportions$Timepoint_categorical <- factor(
  antidefense_systems_sample_counts_mothers_proportions$Timepoint_categorical,
  levels = c("P12", "P28", "B", "M3")
)

pdf('6_ANTIDEFENSE_ANALYSIS/Plots/Mother_Anti-defense_System_Proportions_barplot.pdf', width = 5.5, height = 3.1)
ggplot(antidefense_systems_sample_counts_mothers_proportions, aes(x = Timepoint_categorical, y = proportion, fill = type)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +  
  scale_fill_manual(values = color_mapping) +
  theme_minimal(base_size = 14) + 
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    panel.grid = element_blank(),                        
    panel.border = element_rect(color = "black", fill = NA, size = 0.75), 
    plot.background = element_blank()) +
  ylab("Proportion of Anti-defense Systems (%)") +
  xlab("Timepoint") +
  guides(fill = guide_legend(title = "Anti-defense System"))
dev.off()

#-----------------------------
# Analysis per bacterial host
#-----------------------------

# Add information about the bacterial host
antidefense_systems_per_phage_mothers <- antidefense_systems_per_phage_mothers %>%
  left_join(Viral_metadata_mothers %>%
              select(Virus_ID, Bacterial_genus_host), 
            by = "Virus_ID")

# Calculate the counts of anti-defense systems per bacterial host
antidefense_systems_per_phage_mothers <- antidefense_systems_per_phage_mothers %>%
  filter(!is.na(Bacterial_genus_host)) %>%
  group_by(Bacterial_genus_host, type) %>%
  summarise(count = n(), .groups = 'drop')

# Filter to keep only genera with at least 20 defense systems
filtered_bacterial_genus <- antidefense_systems_per_phage_mothers %>%
  group_by(Bacterial_genus_host) %>%         
  summarise(total_count = sum(count, na.rm = TRUE)) %>%  
  filter(total_count >= 20) %>%               
  pull(Bacterial_genus_host)      

# Filter the data for defense systems types, keeping only relevant genera
antidefense_systems_per_phage_mothers <- antidefense_systems_per_phage_mothers %>%
  filter(Bacterial_genus_host %in% filtered_bacterial_genus) 

# Add Anti-CBASS results (for the plot - even though is absent from mothers)
anti_cbass <- data.frame(Bacterial_genus_host = "Agathobacter", type = "Anti_CBASS",count = 0)
antidefense_systems_per_phage_mothers <- bind_rows(antidefense_systems_per_phage_mothers, anti_cbass)

# Add all combinations (including those with count = 0)
# Again, select only the top 8 most common ADS
antidefense_systems_per_phage_mothers_complete <- expand_grid(
  Bacterial_genus_host = unique(antidefense_systems_per_phage_mothers$Bacterial_genus_host),
  type = unique(antidefense_systems_per_phage_mothers$type)
) %>%
  left_join(antidefense_systems_per_phage_mothers, by = c("Bacterial_genus_host", "type")) %>%
  mutate(count = ifelse(is.na(count), 0, count)) %>%
  filter(type %in% common_ADS_types)  

# Calculate total counts per genus for ordering x-axis
genus_order <- antidefense_systems_per_phage_mothers_complete %>%
  group_by(Bacterial_genus_host) %>%
  summarise(total_count = sum(count, na.rm = TRUE)) %>%
  arrange(desc(total_count)) %>%
  pull(Bacterial_genus_host)

# Set same order as used for infants
type_order <- c("Anti_RM", "Anti_CRISPR", "Anti_RecBCD",  "Anti_Thoeris",
                "Anti_CBASS", "Anti_Retron", "Anti_Dnd", "NADP")

# Convert columns to factors with the desired order
antidefense_systems_per_phage_mothers_complete <- antidefense_systems_per_phage_mothers_complete %>%
  mutate(
    Bacterial_genus_host = factor(Bacterial_genus_host, levels = genus_order),
    type = factor(type, levels = type_order)
  )

# Create the heatmap
pdf('6_ANTIDEFENSE_ANALYSIS/Plots/Mother_Antidefense_System_Counts_per_Type_Bacterial_Host.pdf', width = 7.2, height = 3.5)
ggplot(antidefense_systems_per_phage_mothers_complete, aes(x = Bacterial_genus_host, y = type, fill = log10(count + 1))) +  
  geom_tile(color = "lightgrey", linewidth = 0.5, linetype = 1) +
  scale_fill_gradientn(colors = c("white", "lightblue", "maroon"), 
                       values = scales::rescale(c(0, 0.5, 1)), 
                       limits = c(0, log10(max(antidefense_systems_per_phage_mothers_complete$count + 1))),
                       name = "Log(Count)") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.title.x = element_text(size = 14),
    axis.text.y = element_text(size = 12),
    axis.title.y = element_text(size = 14, margin = margin(r = 10)),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    strip.text = element_text(size = 12, face = "bold"),
    panel.border = element_rect(color = "black", fill = NA, size = 1)  
  ) +
  ylab("Anti-defense System") + 
  xlab("Bacterial Genus Host") + 
  guides(fill = guide_colorbar(title = "Log10(Count)"))
dev.off()


##************************************************************************
# 2. Enrichment tests
##************************************************************************

#######################################
# 2.A. ADS systems in Mothers vs Infants
#######################################

# Check if any anti-defense system is enriched (more common) in maternal/infant vOTUs

# Merge infant and maternal ADS system counts
merged_ADS_counts <- merge(antidefense_systems_counts_infants, 
                antidefense_systems_counts_mothers, 
                by = "type", all = TRUE)
colnames(merged_ADS_counts) <- c("type", "infant_n", "mother_n")

# Add total counts
merged_ADS_counts$Total <- rowSums(merged_ADS_counts[, c("infant_n", "mother_n")], na.rm = TRUE)
merged_ADS_counts_table <- merged_ADS_counts
  
# Get total number for mothers and infants
total_infant <- sum(merged_ADS_counts$infant_n, na.rm = TRUE)
total_mother <- sum(merged_ADS_counts$mother_n, na.rm = TRUE)

# Exclude those with less than 20 total counts for the enrichment
merged_ADS_counts <- subset(merged_ADS_counts, Total >= 20)

# Initialize result dataframe
ADS_enrichment_results <- data.frame(
  type = merged_ADS_counts$type,
  p_value = NA,
  test_used = NA,
  enrichment = NA  # "Infants", "Mothers", or "None"
)

# Perform the statistical test for each system
for (i in 1:nrow(merged_ADS_counts)) {
  # Get counts of system in each group
  infant_with_system <- merged_ADS_counts$infant_n[i]
  mother_with_system <- merged_ADS_counts$mother_n[i]
  
  # Get counts of systems NOT of this type
  infant_without_system <- total_infant - infant_with_system
  mother_without_system <- total_mother - mother_with_system
  
  # Construct the contingency table
  contingency <- matrix(c(infant_with_system, mother_with_system,
                          infant_without_system, mother_without_system),
                        nrow = 2, byrow = TRUE)
  
  # Use Fisher's exact test
  test <- fisher.test(contingency)
  ADS_enrichment_results$test_used[i] <- "Fisher"
  
  # Store p-value
  ADS_enrichment_results$p_value[i] <- test$p.value
  
  # Determine which group is enriched
  prop_infant <- infant_with_system / total_infant
  prop_mother <- mother_with_system / total_mother
  if (test$p.value < 0.05) {
    ADS_enrichment_results$enrichment[i] <- ifelse(prop_infant > prop_mother, "Infants", "Mothers")
  } else {
    ADS_enrichment_results$enrichment[i] <- "None"
  }
}

# Estimate FDR and sort results
ADS_enrichment_results$FDR <- p.adjust(ADS_enrichment_results$p_value, method = "BH")
ADS_enrichment_results <- ADS_enrichment_results[order(ADS_enrichment_results$FDR), ]

# Generate a table with all the ADS counts together with enrichment results
ADS_final_table <- merged_ADS_counts_table %>%
  left_join(
    ADS_enrichment_results,
    by = "type"
  )

#######################################
# 2.B. "Anti_RecBCD" and "Anti_CBASS" enrichment per host
#######################################

# Check if "Anti_RecBCD" and "Anti_CBASS" are enriched in any bacterial host
# For this analysis, use hosts with at least 50 ADS systems detected across phages predicted to infect them
systems_to_test <- c("Anti_RecBCD", "Anti_CBASS")
infant_ADS_by_host <- antidefense_systems_per_phage_infants

# Get bacterial genera with total counts >= 50 (across all types)
filtered_bacterial_genus <- infant_ADS_by_host %>%
  group_by(Bacterial_genus_host) %>%
  summarise(total_count = sum(count)) %>%
  filter(total_count >= 50) %>%
  pull(Bacterial_genus_host)

infant_ADS_enrichment_by_host <- expand.grid(system = systems_to_test,
                                             genus = filtered_bacterial_genus,
                                             stringsAsFactors = FALSE) %>%
  rowwise() %>%
  mutate(test_result = list({
    a <- infant_ADS_by_host %>% filter(type == system, Bacterial_genus_host == genus) %>% pull(count) %>% sum()
    b <- infant_ADS_by_host %>% filter(type != system, Bacterial_genus_host == genus) %>% pull(count) %>% sum()
    c <- infant_ADS_by_host %>% filter(type == system, Bacterial_genus_host != genus) %>% pull(count) %>% sum()
    d <- infant_ADS_by_host %>% filter(type != system, Bacterial_genus_host != genus) %>% pull(count) %>% sum()
    
    # Generate the # 2x2 contingency table
    contingency <- matrix(c(a, c, b, d), nrow = 2, byrow = TRUE)
    # Perform the statistical test
    fisher.test(contingency, simulate.p.value = TRUE, B = 1e6)
  })) %>%
  ungroup() %>%
  mutate(
    p_value = sapply(test_result, function(x) x$p.value),
    odds_ratio = sapply(test_result, function(x) x$estimate)
  ) %>%
  select(-test_result) %>%
  mutate(FDR = p.adjust(p_value, method = "fdr")) %>%
  arrange(FDR)

#****************
# Write results
#****************
# Anti-defense systems results
write.table(antidefense_systems_counts_infants,
            "6_ANTIDEFENSE_ANALYSIS/Antidefense_systems_counts_infants.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(antidefense_systems_counts_mothers,
            "6_ANTIDEFENSE_ANALYSIS/Antidefense_systems_counts_mothers.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(ADS_enrichment_results,
            "6_ANTIDEFENSE_ANALYSIS/Antidefense_systems_type_enrichment.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(ADS_enrichment_results,
            "6_ANTIDEFENSE_ANALYSIS/Antidefense_systems_type_enrichment.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(ADS_final_table,
            "6_ANTIDEFENSE_ANALYSIS/Antidefense_systems_counts_and_enrichment.txt",
            sep = "\t", row.names = F, quote = FALSE) 
write.table(infant_ADS_enrichment_by_host,
            "6_ANTIDEFENSE_ANALYSIS/Antidefense_systems_enrichment_by_host_INFANTS.txt",
            sep = "\t", row.names = F, quote = FALSE) 

# Viral metadata
write.table(Viral_metadata,"Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Viral_metadata_infants,"Metadata_NEXT/LLNEXT_Viral_Metadata_INFANTS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Viral_metadata_mothers,"Metadata_NEXT/LLNEXT_Viral_Metadata_MOTHERS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 

