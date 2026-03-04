################################################################################
##### LL-NEXT: Host Assignment Analysis
### Author(s):Asier Fernández-Pato
### Last updated: 1st July, 2025
################################################################################

#****************
# Load modules
#****************
library(dplyr)
library(tidyverse)
library(phyloseq)
library(microbiome)
library(ggplot2)
library(ggpattern)

# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_analysis/VIRUSES/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_07052025.txt")
Sample_metadata <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt")
Sample_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_07052025.txt")
Sample_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_07052025.txt")
Viral_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")

# Reorder Timepoint categorical variable
Sample_metadata$Timepoint_categorical <- factor(
  Sample_metadata$Timepoint_categorical,
  levels = c("P12", "P28", "B", "W2", "M1", "M2", "M3", "M6", "M9", "M12")
)

##************************************************************************
# 1. Host Assignment Analysis: iPHOP results preprocessing
##************************************************************************

# For the processing of the host prediction results, we follow the next steps:
### - 1) Default DB: Get prediction with max score per virus (could be > 1 different predictions with same score)
### - 2) Enriched DB: Get prediction with max score per virus (could be > 1 different predictions with same score)
### - 3) Combine results: Combine (1) and (2), retaining predictions with max score (could be > 1 different predictions with same score)
### - 4) Extract only 1 prediction per phage by:
##        - Extracting genus-level and family-level prediction if only 1 prediction with max score is present
##        - Extracting the prediction of the common ancestor (family-level) if > 1 different prediction with max score is present
##        - Set to NA cases with different family-level predictions with the same score

#Note: iPHOP was run on a set of 30,824 phages (the final selection included 30,847 after manual curation of prokaryotic origin)

#################
# A) Default DB
#################

#A) Read iPHOP prediction output at genus level with default DB
viral_hosts_df <- read.delim("8_HOST_PREDICTION/Host_prediction_to_genus_m90_default_db.csv", 
                                sep = ",", header = T, check.names = F)

# Filter iPHOP predictions to get only the predictions with max score per virus (could be > 1 with same score)
viral_hosts_df <- viral_hosts_df %>%
  group_by(Virus) %>%
  filter(`Confidence score` == max(`Confidence score`)) %>%
  ungroup() %>%
  mutate(Source = "df")

viral_hosts_df <- viral_hosts_df[,c(1,3,4)]
colnames(viral_hosts_df) <- c("Virus_ID", "Bacterial_host", "Score") 

# With default DB:
# - Predictions for 27890/30824 phages (90.48%)
# - 824 phages have more than 1 genus-level prediction with the same score

#################
# B) Enriched DB
#################

#B) Read iPHOP prediction output at genus level with enriched DB
viral_hosts_enr <- read.delim("8_HOST_PREDICTION/Host_prediction_to_genus_m90_enriched_db.csv", 
                          sep = ",", header = T, check.names = F)

# Filter iPHOP predictions to get only the predictions with max score per virus (could be > 1 with same score)
viral_hosts_enr <- viral_hosts_enr %>%
  group_by(Virus) %>%
  filter(`Confidence score` == max(`Confidence score`)) %>%
  ungroup() %>%
  mutate(Source = "df")

viral_hosts_enr <- viral_hosts_enr[,c(1,3,4)]
colnames(viral_hosts_enr) <- c("Virus_ID", "Bacterial_host", "Score") 

# With the enriched DB:
# - Predictions for 27826/30824 phages (90.27%)
# - 812 phages have more than 1 genus-level prediction with the same score

#################
# C) Combined results
#################

# Select prediction with highest scores for the previous 2 prediction approaches
viral_hosts_combined <- rbind(viral_hosts_df, viral_hosts_enr)

viral_hosts_combined <- viral_hosts_combined %>%
  group_by(Virus_ID) %>%
  filter(Score == max(Score)) %>%
  ungroup() %>%
  distinct()

# Check which phages have >1 prediction with the same score
duplicated_status <- viral_hosts_combined$Virus_ID %in% viral_hosts_combined$Virus_ID[duplicated(viral_hosts_combined$Virus_ID)]

# Add a new column indicating prediction status
viral_hosts_combined <- viral_hosts_combined %>%
  mutate(Multiple_prediction = ifelse(duplicated_status, "Yes", "No"))

# Extract family and genus-level hosts
viral_hosts_combined$Bacterial_family_host <- str_extract(viral_hosts_combined$Bacterial_host, "(?<=f__)[^;]+")
viral_hosts_combined$Bacterial_genus_host <- str_extract(viral_hosts_combined$Bacterial_host, "(?<=g__)[^;]+")

# Generate final prediction table (with 1 prediction per virus):
## - Extract genus-level and family-level prediction if only 1 prediction with max score is present
## - Extract the prediction of the common ancestor (family-level) if >1 prediction with max score is present
viral_hosts_combined <- viral_hosts_combined %>%
  group_by(Virus_ID) %>%
  mutate(
    # If Bacterial_genus_host values are different for a Virus_ID, set them to NA 
    Bacterial_genus_host = ifelse(length(unique(Bacterial_genus_host)) > 1, NA, Bacterial_genus_host),
    
    # If Bacterial_family_host values are different for a Virus_ID, set them to NA
    Bacterial_family_host = ifelse(length(unique(Bacterial_family_host)) > 1, NA, Bacterial_family_host)
  ) %>%
  ungroup() %>%
  distinct(Virus_ID, Multiple_prediction, Bacterial_family_host, Bacterial_genus_host, Score, .keep_all = TRUE)  # Remove duplicates

# Process Bacterial_host values
viral_hosts_combined <- viral_hosts_combined %>%
  mutate(
    # If Bacterial_family_host is NA, set Bacterial_host to NA
    Bacterial_host = ifelse(is.na(Bacterial_family_host), NA, Bacterial_host),
    
    # If Bacterial_genus_host is NA, but Bacterial_family_host is not, keep only the family
    Bacterial_host = ifelse(is.na(Bacterial_genus_host) & !is.na(Bacterial_family_host),
                            gsub("(.*;f__[A-Za-z0-9_-]+);.*", "\\1", Bacterial_host), Bacterial_host),
  )

# Merge with Virus metadata
viral_hosts_combined_to_merge <- viral_hosts_combined[c("Virus_ID", "Bacterial_genus_host", "Bacterial_family_host")]
Viral_metadata <- left_join(Viral_metadata, viral_hosts_combined_to_merge, by = "Virus_ID")

# Set host columns as a factor (ordering them by occurence)
Viral_metadata$Bacterial_family_host <- factor(Viral_metadata$Bacterial_family_host)
Viral_metadata$Bacterial_genus_host <- factor(Viral_metadata$Bacterial_genus_host)

# With the combined results:
# - Genus-level predictions for 27632/30824 phages (89.65%)  (27631/30847 from the final selection: 89.57%)
# - Family-level predictions for 28434/30824 phages (92.26%)  (28433/30847 from the final selection: 92.17%)
# - 2387/30824 phages (7.74%) with no genus/family-level prediction (2414/30847 from the final selection: 7.82%)

# Generate a bar plot with the distribution of hosts (top 5)
top5_genera <- Viral_metadata %>%
  filter(!is.na(Bacterial_genus_host)) %>%
  count(Bacterial_genus_host, sort = TRUE) %>%
  slice_head(n = 5)

palette <- c("#A9C6DF", "#D3A9A2","#B7C4A1","#CAB7D4", "#E0CDA8")
names(palette) <- top5_genera$Bacterial_genus_host

pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/Host_top5_barplot.pdf', width = 4, height = 2.7)
ggplot(top5_genera, aes(x = reorder(Bacterial_genus_host, n), y = n, fill = Bacterial_genus_host)) +
  geom_bar(stat = "identity", width = 0.4, color = "#707070") + 
  scale_fill_manual(values = palette) +
  coord_flip() +
  ylab("Number of vOTUs") +
  theme_minimal(base_size = 16) +  
  theme(
    axis.text = element_text(size = 16),
    axis.title = element_text(size = 17),
    axis.text.y = element_text(),
    axis.title.y = element_blank(),
    legend.position = "none",
    panel.grid = element_blank(),                         
    panel.border = element_rect(color = "black", fill = NA, size = 0.75), 
    plot.background = element_blank(),
    axis.ticks.x = element_line(color = "black"))
dev.off()

##################################
# Generation of Phyloseq object
##################################

# Split host taxonomy
tax_levels <- strsplit(as.character(viral_hosts_combined$Bacterial_host), ";")
tax_levels <- lapply(tax_levels, function(x) substring(x, 4))

# Determine the maximum number of levels in the taxonomy column and fill empty values with NA
max_levels <- max(sapply(tax_levels, length))
tax_levels <- lapply(tax_levels, function(x) c(x, rep(NA, max_levels - length(x))))

# Bind the taxonomy vectors to the original data frame
taxonomy_hosts  <- cbind(viral_hosts_combined, do.call("rbind", tax_levels))
taxonomy_hosts[,c(2:6)] <- NULL

# Rename the new columns
colnames(taxonomy_hosts)[2:7] <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus")

# As some viruses do not have host taxonomy, add them manually to the taxonomy table
# Create a vector of Virus_IDs with no match in taxonomy_viruses
taxonomy_hosts <- taxonomy_hosts[taxonomy_hosts$Virus_ID %in% rownames(Abundance_table), ] 
missing_IDs <- setdiff(rownames(Abundance_table), taxonomy_hosts$Virus_ID)

# Create a new dataframe with NA values for missing Virus_IDs
missing_hosts_viruses <- data.frame(Virus_ID = missing_IDs)
missing_hosts_viruses[, names(taxonomy_hosts)[2:7]] <- NA

# Combine taxonomy_viruses and missing_viruses
taxonomy_hosts <- rbind(taxonomy_hosts, missing_hosts_viruses)

# Reorder the rows based on Virus_ID matching row_names
taxonomy_hosts <- data.frame(taxonomy_hosts[match(rownames(Abundance_table), taxonomy_hosts$Virus_ID), ])

# Set virus IDs as rownames 
rownames(taxonomy_hosts) <- taxonomy_hosts$Virus_ID
taxonomy_hosts$Virus_ID <- NULL

#**Due to differences in taxonomy assignment of MAGs in enriched DB some viruses are assigned = "Genus" but != "Phylum" (changed in GTDB tax)
# I manually adapt taxonomy_hosts table it to make it uniform (and avoid problems for visualization with phyloseq)
#taxonomy_hosts <- read.delim("6_HOST_ASSIGNMENT/taxonomy_hosts.tsv", row.names = 1)
#colnames(taxonomy_hosts) <- c("Domain","Phylum","Class","Order","Family","Genus")
#taxonomy_hosts <- taxonomy_hosts[match(rownames(Abundance_table), rownames(taxonomy_hosts)), ]

# Transform dataframes into tables (for Phyloseq)
Abundance_table_matrix <- as.matrix(Abundance_table)
taxonomy_hosts_matrix <- as.matrix(taxonomy_hosts)

# Generate phyloseq object
vOTU <- otu_table(Abundance_table_matrix, taxa_are_rows = TRUE)
TAX <- tax_table(taxonomy_hosts_matrix)
samples <- sample_data(Sample_metadata)
sample_names(samples) <- Sample_metadata$NG_ID

Phyloseq_hosts <- phyloseq(vOTU, TAX, samples)


##************************************************************************
# 2. Host Taxonomic Analysis in Infants
##************************************************************************

# Generate Virus metadata for those present in infant samples
Abundance_table_infants <- Abundance_table[,Sample_metadata$Type == "infant"]
viruses_present_infants <- rownames(Abundance_table_infants)[rowSums(Abundance_table_infants)>0]
Viral_metadata_infants <- Viral_metadata[Viral_metadata$Virus_ID %in% viruses_present_infants, ]

# Generate a barplot of the host assignment of vOTUs at family level
Viral_family_host_counts_infants <- Viral_metadata_infants %>%
  count(Bacterial_family_host) %>%
  na.omit()

Viral_family_host_counts_infants <- head(arrange(Viral_family_host_counts_infants, desc(n)), 8) # select top 8

pdf('8_HOST_PREDICTION/Plots/Family_host_barplot_infants.pdf', width=4.5, height=3.1)
ggplot(Viral_family_host_counts_infants, aes(x = reorder(Bacterial_family_host, -n), y = n)) +
  geom_bar_pattern(stat = "identity", width = 0.7, pattern = "stripe", pattern_color = "#707070", 
                   pattern_fill = "#707070", pattern_density = 0.15, pattern_angle = 45, color = "#707070", fill = "#D3D3D3" 
  ) +
  coord_flip() +
  ylab("Number of vOTUs") +
  scale_fill_manual(values = rep("#D3D3D3", length(unique(Viral_family_host_counts_infants$Bacterial_family_host)))) +  
  theme_classic() +
  theme(
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 15),
    axis.text.y = element_text(face = "italic"),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none"
  ) 
dev.off()

# Generate a barplot of the host assignment of vOTUs at genus level
Viral_genus_host_counts_infants <- Viral_metadata_infants %>%
  count(Bacterial_genus_host) %>%
  na.omit()

Viral_genus_host_counts_infants <- head(arrange(Viral_genus_host_counts_infants, desc(n)), 8) # select top 8

pdf('8_HOST_PREDICTION/Plots/Genus_host_barplot_infants.pdf', width=4.5, height=3.1)
ggplot(Viral_genus_host_counts_infants, aes(x = reorder(Bacterial_genus_host, -n), y = n)) +
  geom_bar_pattern(stat = "identity", width = 0.7, pattern = "stripe", pattern_color = "#707070", 
    pattern_fill = "#707070", pattern_density = 0.15, pattern_angle = 45, color = "#707070", fill = "#D3D3D3" 
  ) +
  coord_flip() +
  ylab("Number of vOTUs") +
  scale_fill_manual(values = rep("#D3D3D3", length(unique(Viral_genus_host_counts_infants$Bacterial_genus_host)))) +  
  theme_classic() +
  theme(
    axis.text = element_text(size = 14, colour = "black"),
    axis.title = element_text(size = 15),
    axis.text.y = element_text(face = "italic"),
    axis.title.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "none"
  ) 
dev.off()


###################################
# Host abundance analysis
##################################

# Select only infant samples
Phyloseq_hosts_infants <- subset_samples(Phyloseq_hosts, Type =="Infant")

#********
# Family-level
#********
# Generate family-level abundances
Phyloseq_hosts_infants_family <- Phyloseq_hosts_infants %>% 
  subset_taxa(!is.na(Family)) %>% 
  aggregate_taxa(level = "Family")  

# Identify the top 10 genera by total abundance
top_family <- names(sort(taxa_sums(Phyloseq_hosts_infants_family), decreasing = TRUE)[1:10])

# Transform data to compositional for relative abundance
Phyloseq_hosts_infants_family_rel <- microbiome::transform(Phyloseq_hosts_infants_family, "compositional")

# Convert the phyloseq object to a data frame 
family_infant_df <- psmelt(Phyloseq_hosts_infants_family_rel)

# Create an "Other" category for genera outside of the top 10
family_infant_df  <- family_infant_df %>%
  mutate(Family = ifelse(Family %in% top_family, Family, "Other")) %>%
  group_by(Sample, Family) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup()

# Calculate the genera with the highest relative abundance
total_abundance <- family_infant_df %>%
  group_by(Family) %>%
  summarize(Total_Abundance = sum(Abundance)) %>%
  arrange(desc(Total_Abundance))

# Set the order of Genus factor based on the relative abundance (with "Other" on top)
family_infant_df$Family <- factor(family_infant_df$Family, 
                                  levels = c("Other", setdiff(unique(family_infant_df$Family), "Other")))

# Set colors for each family 
family_colors <- c("#BF9000", "#5B9279", "#CC7722", "#80B3AD", "#91A7B3",
                   "#D68A9E", "#E3B96F", "#A3C97B", "#9C8AA5", "#B3AF8F", "#D3D3D3")
family_colors <- adjustcolor(family_colors, alpha.f = 0.6)
names(family_colors) <- c(setdiff(unique(family_infant_df$Family), "Other"), "Other") 


# Include Timepoint_categorical in the data frame
Timepoint <- Sample_metadata_infants[,c("NG_ID", "Timepoint_categorical")]
colnames(Timepoint)[1] <- "Sample"
family_infant_df <- left_join(family_infant_df, Timepoint, by = "Sample")

# Ensure Timepoint_categorical is a factor and reorder levels if needed
family_infant_df$Timepoint_categorical <- factor(family_infant_df$Timepoint_categorical, 
                                                 levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12"))



# Plot the abundance data at the family level including the "Other" category
pdf('8_HOST_PREDICTION/Plots/Infant_Family_Abundance_barplot.pdf', width = 6, height = 3.1)
ggplot(family_infant_df, aes(x = Timepoint_categorical, y = Abundance)) +
  geom_bar(stat = "identity", position = "fill", color = NA, aes(fill = Family)) +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +  
  scale_fill_manual(values = family_colors) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13, face = "italic"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        panel.background = element_blank(),
        plot.background = element_blank()) +
  ylab("Relative Abundance (%)") +
  xlab("Timepoint") +
  guides(fill = guide_legend(title = "Family"))
dev.off()

#********
# Genus-level
#********

# Generate genus-level abundances
Phyloseq_hosts_infants_genus <- Phyloseq_hosts_infants %>% 
  subset_taxa(!is.na(Genus)) %>% 
  aggregate_taxa(level = "Genus")  

# Identify the top 10 genera by total abundance
top_genus <- names(sort(taxa_sums(Phyloseq_hosts_infants_genus), decreasing = TRUE)[1:10])

# Transform data to compositional for relative abundance
Phyloseq_hosts_infants_genus_rel <- microbiome::transform(Phyloseq_hosts_infants_genus, "compositional")

# Convert the phyloseq object to a data frame 
genus_infant_df <- psmelt(Phyloseq_hosts_infants_genus_rel)

# Create an "Other" category for genera outside of the top 10
genus_infant_df  <- genus_infant_df %>%
  mutate(Genus = ifelse(Genus %in% top_genus, Genus, "Other")) %>%
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup()

# Calculate the genera with the highest relative abundance
total_abundance <- genus_infant_df %>%
  group_by(Genus) %>%
  summarize(Total_Abundance = sum(Abundance)) %>%
  arrange(desc(Total_Abundance))

# Set the order of Genus factor based on the relative abundance (with Other on top)
genus_infant_df$Genus <- factor(genus_infant_df$Genus, levels = c("Other", setdiff(total_abundance$Genus, "Other")))

# Set colors for each genus 
genus_colors <- c("#74B2A6", "#95AFD1", "#A5A0C2", "#E57164", 
                     "#E5A356", "#C4A484", "#F6CDE0", "#A86FAE",
                     "#B8DDB5", "#FFE57F","#D3D3D3")
genus_colors <- adjustcolor(genus_colors, alpha.f = 0.7)
names(genus_colors) <- c(setdiff(unique(genus_infant_df$Genus), "Other"), "Other") 


# Include Timepoint_categorical in the data frame
Timepoint <- Sample_metadata_infants[,c("NG_ID", "Timepoint_categorical")]
colnames(Timepoint)[1] <- "Sample"
genus_infant_df <- left_join(genus_infant_df, Timepoint, by = "Sample")

# Ensure Timepoint_categorical is a factor and reorder levels if needed
genus_infant_df$Timepoint_categorical <- factor(genus_infant_df$Timepoint_categorical, 
                                                 levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12"))


# Plot the abundance data at the genus level including the "Other" category
pdf('8_HOST_PREDICTION/Plots/Infant_Genus_Abundance_barplot.pdf', width = 6, height = 3.1)
ggplot(genus_infant_df, aes(x = Timepoint_categorical, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +  
  scale_fill_manual(values = genus_colors) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13, face = "italic"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        panel.background = element_blank(),
        plot.background = element_blank()) +
  ylab("Relative Abundance (%)") +
  xlab("Timepoint") +
  guides(fill = guide_legend(title = "Genus"))
dev.off()


##************************************************************************
# 3. Host Taxonomic Analysis in Mothers
##************************************************************************

# Generate Virus metadata for those present in mother samples
Abundance_table_mothers <- Abundance_table[,Sample_metadata$Type == "Mother"]
viruses_present_mothers <- rownames(Abundance_table_mothers)[rowSums(Abundance_table_mothers)>0]
Viral_metadata_mothers <- Viral_metadata[Viral_metadata$Virus_ID %in% viruses_present_mothers, ]

# Generate a barplot of the host assignment of vOTUs at family level
Viral_family_host_counts_mothers <- Viral_metadata_mothers %>%
  count(Bacterial_family_host) %>%
  na.omit()

Viral_family_host_counts_mothers <- head(arrange(Viral_family_host_counts_mothers, desc(n)), 8) # select top 8

pdf('8_HOST_PREDICTION/Plots/Family_host_barplot_mothers.pdf', width=4.5, height=3.1)
ggplot(Viral_family_host_counts_mothers, aes(x = reorder(Bacterial_family_host, -n), y = n)) +
  geom_bar_pattern(stat = "identity", width = 0.7, pattern = "stripe", pattern_color = "#707070", 
                   pattern_fill = "#707070", pattern_density = 0.15, pattern_angle = 45, color = "#707070", fill = "#D3D3D3" 
  ) +
  coord_flip() +
  ylab("Number of vOTUs") +
  scale_fill_manual(values = rep("#D3D3D3", length(unique(Viral_family_host_counts_mothers$Bacterial_family_host)))) + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(face = "italic"),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none")
dev.off()


# Generate a barplot of the host assignment of vOTUs at genus level
Viral_genus_host_counts_mothers <- Viral_metadata_mothers %>%
  count(Bacterial_genus_host) %>%
  na.omit()

Viral_genus_host_counts_mothers <- head(arrange(Viral_genus_host_counts_mothers, desc(n)), 8) # select top 8

pdf('8_HOST_PREDICTION/Plots/Genus_host_barplot_mothers.pdf', width=4.5, height=3.1)
ggplot(Viral_genus_host_counts_mothers, aes(x = reorder(Bacterial_genus_host, -n), y = n)) +
  geom_bar_pattern(stat = "identity", width = 0.7, pattern = "stripe", pattern_color = "#707070", 
                   pattern_fill = "#707070", pattern_density = 0.15, pattern_angle = 45, color = "#707070", fill = "#D3D3D3") +
  coord_flip() +
  ylab("Number of vOTUs") +
  scale_fill_manual(values = rep("#D3D3D3", length(unique(Viral_genus_host_counts_mothers$Bacterial_genus_host)))) + 
  theme_classic() +
  theme(axis.text = element_text(size = 14, colour = "black"),
        axis.title = element_text(size = 15),
        axis.text.y = element_text(face = "italic"),
        axis.title.y = element_blank(),
        axis.line.y = element_blank(),
        legend.position = "none") 
dev.off()

##################################
# Host abundance analysis
##################################

# Select only mother samples
Phyloseq_hosts_mothers <- subset_samples(Phyloseq_hosts, Type =="Mother")

#********
# Family-level
#********
# Generate family-level abundances
Phyloseq_hosts_mothers_family <- Phyloseq_hosts_mothers %>% 
  subset_taxa(!is.na(Family)) %>% 
  aggregate_taxa(level = "Family")  

# Identify the top 10 genera by total abundance
top_family <- names(sort(taxa_sums(Phyloseq_hosts_mothers_family), decreasing = TRUE)[1:10])

# Transform data to compositional for relative abundance
Phyloseq_hosts_mothers_family_rel <- microbiome::transform(Phyloseq_hosts_mothers_family, "compositional")

# Convert the phyloseq object to a data frame 
family_mother_df <- psmelt(Phyloseq_hosts_mothers_family_rel)

# Create an "Other" category for genera outside of the top 10
family_mother_df  <- family_mother_df %>%
  mutate(Family = ifelse(Family %in% top_family, Family, "Other")) %>%
  group_by(Sample, Family) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup()

# Calculate the genera with the highest relative abundance
total_abundance <- family_mother_df %>%
  group_by(Family) %>%
  summarize(Total_Abundance = sum(Abundance)) %>%
  arrange(desc(Total_Abundance))

# Set the order of Genus factor based on the relative abundance (with Other on top)
family_mother_df$Family <- factor(family_mother_df$Family, levels = c("Other", setdiff(total_abundance$Family, "Other")))

# Set colors for each family 
mother_families <- unique(family_mother_df$Family)
new_families <- setdiff(mother_families, names(family_colors))
new_family_colors <- adjustcolor(c("#8C6A4A","#4A7A8C","#8C4A72", "#4A8C6A","#6A4A8C"),
                                 alpha.f = 0.7)  
names(new_family_colors) <- new_families
family_colors_mothers <- c(family_colors, new_family_colors)
family_colors_mothers <- family_colors_mothers[names(family_colors_mothers) %in% mother_families]

# Include Timepoint_categorical in the data frame
Timepoint <- Sample_metadata_mothers[,c("NG_ID", "Timepoint_categorical")]
colnames(Timepoint)[1] <- "Sample"
family_mother_df <- left_join(family_mother_df, Timepoint, by = "Sample")

# Ensure Timepoint_categorical is a factor and reorder levels if needed
family_mother_df$Timepoint_categorical <- factor(family_mother_df$Timepoint_categorical, 
                                                 levels = c("P12", "P28", "B", "M3"))


# Plot the abundance data at the family level including the "Other" category
pdf('8_HOST_PREDICTION/Plots/Mother_Family_Abundance_barplot.pdf', width = 6, height = 3.1)
ggplot(family_mother_df, aes(x = Timepoint_categorical, y = Abundance, fill = Family)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +  
  scale_fill_manual(values = family_colors_mothers) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13, face = "italic"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        panel.background = element_blank(),
        plot.background = element_blank()) +
  ylab("Relative Abundance (%)") +
  xlab("Timepoint") +
  guides(fill = guide_legend(title = "Family"))
dev.off()


#********
# Genus-level
#********

# Generate genus-level abundances
Phyloseq_hosts_mothers_genus <- Phyloseq_hosts_mothers %>% 
  subset_taxa(!is.na(Genus)) %>% 
  aggregate_taxa(level = "Genus")  

# Identify the top 10 genera by total abundance
top_genus <- names(sort(taxa_sums(Phyloseq_hosts_mothers_genus), decreasing = TRUE)[1:10])

# Transform data to compositional for relative abundance
Phyloseq_hosts_mothers_genus_rel <- microbiome::transform(Phyloseq_hosts_mothers_genus, "compositional")

# Convert the phyloseq object to a data frame 
genus_mother_df <- psmelt(Phyloseq_hosts_mothers_genus_rel)

# Create an "Other" category for genera outside of the top 10
genus_mother_df  <- genus_mother_df %>%
  mutate(Genus = ifelse(Genus %in% top_genus, Genus, "Other")) %>%
  group_by(Sample, Genus) %>%
  summarise(Abundance = sum(Abundance)) %>%
  ungroup()

# Calculate the genera with the highest relative abundance
total_abundance <- genus_mother_df %>%
  group_by(Genus) %>%
  summarize(Total_Abundance = sum(Abundance)) %>%
  arrange(desc(Total_Abundance))

# Set the order of Genus factor based on the relative abundance (with Other on top)
genus_mother_df$Genus <- factor(genus_mother_df$Genus, levels = c("Other", setdiff(total_abundance$Genus, "Other")))

# Set colors for each family 
mother_genera <- unique(genus_mother_df$Genus)
new_genera <- setdiff(mother_genera, names(genus_colors))
new_genera_colors <- adjustcolor(c("#8C6A4A","#4A7A8C","#8C4A72", "#4A8C6A","#6A4A8C"),
                                 alpha.f = 0.7)  
names(new_genera_colors) <- new_genera
genus_colors_mothers <- c(genus_colors, new_genera_colors)
genus_colors_mothers <- genus_colors_mothers[names(genus_colors_mothers) %in% mother_genera]


# Include Timepoint_categorical in the data frame
Timepoint <- Sample_metadata_mothers[,c("NG_ID", "Timepoint_categorical")]
colnames(Timepoint)[1] <- "Sample"
genus_mother_df <- left_join(genus_mother_df, Timepoint, by = "Sample")

# Ensure Timepoint_categorical is a factor and reorder levels if needed
genus_mother_df$Timepoint_categorical <- factor(genus_mother_df$Timepoint_categorical, 
                                                levels = c("P12", "P28", "B", "M3"))


# Plot the abundance data at the genus level including the "Other" category
pdf('8_HOST_PREDICTION/Plots/Mother_Genus_Abundance_barplot.pdf', width = 6, height = 3.1)
ggplot(genus_mother_df, aes(x = Timepoint_categorical, y = Abundance, fill = Genus)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_y_continuous(labels = function(x) paste0(x * 100, "%")) +  
  scale_fill_manual(values = genus_colors_mothers) +
  theme_classic() +
  theme(axis.text = element_text(size = 14),
        axis.title = element_text(size = 15),
        legend.title = element_text(size = 15),
        legend.text = element_text(size = 13, face = "italic"),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        panel.background = element_blank(),
        plot.background = element_blank()) +
  ylab("Relative Abundance (%)") +
  xlab("Timepoint") +
  guides(fill = guide_legend(title = "Genus"))
dev.off()

##*************
# Save output
#**************
# Save virus metadata
write.table(Viral_metadata,"Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Viral_metadata_infants,"Metadata_NEXT/LLNEXT_Viral_Metadata_INFANTS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Viral_metadata_mothers,"Metadata_NEXT/LLNEXT_Viral_Metadata_MOTHERS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 

# Save genus-level abundances
Bacterial_genus_host_vOTU_abundance_infants <- data.frame(t(otu_table(Phyloseq_hosts_infants_genus)))
write.table(Bacterial_genus_host_vOTU_abundance_infants, "Abundance_table/Bacterial_genus_host_vOTU_abundance_INFANTS_17092025.txt",
            sep = "\t", row.names = T, quote = FALSE)

Bacterial_genus_host_vOTU_abundance_mothers <- data.frame(t(otu_table(Phyloseq_hosts_mothers_genus)))
write.table(Bacterial_genus_host_vOTU_abundance_mothers, "Abundance_table/Bacterial_genus_host_vOTU_abundance_MOTHERS_17092025.txt",
            sep = "\t", row.names = T, quote = FALSE)

write.table(Abundance_table_mothers_host, "Abundance_table/Bacterial_genus_host_vOTU_abundance_MOTHERS_17092025.txt",
            sep = "\t", row.names = T, quote = FALSE)
