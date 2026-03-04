################################################################################
##### LL-NEXT: Alpha and Beta diversity analysis
### Author(s):Asier Fernández-Pato
### Last updated: 19th December, 2025
################################################################################


#****************
# Load modules
#****************
library(vegan)       
library(lme4)        
library(lmerTest)    
library(dplyr)       
library(ggplot2)     
library(ggExtra)    
library(reshape2)

#****************
# Define functions
#****************

# Function to perform a permutation test to compare distances 
perm_test <- function(group1, group2, data, group_col = "Type", 
                              timepoint_pair = NULL, n_perm = 10000, pairwise = TRUE) {
  # If pairwise comparison (specific timepoint pairs), use Timepoint_pair
  if (pairwise) {
    if (is.null(timepoint_pair)) {
      stop("Please specify 'timepoint_pair' for pairwise comparisons.")
    }
    
    # Extract the distances for the specified pair (group_col can be "Type", "Relatedness", etc.)
    g1 <- data$Distance[data$Timepoint_pair == timepoint_pair[1] & data[[group_col]] == group1]
    g2 <- data$Distance[data$Timepoint_pair == timepoint_pair[2] & data[[group_col]] == group2]
    
  } else {  # If comparing overall distances between two groups
    # Extract distances for the two groups
    g1 <- data$Distance[data[[group_col]] == group1]
    g2 <- data$Distance[data[[group_col]] == group2]
  }
  
  # Calculate observed difference in medians
  obs_diff <- median(g1) - median(g2)
  
  # Combine the two groups together
  combined <- c(g1, g2)
  labels <- c(rep("g1", length(g1)), rep("g2", length(g2)))
  
  # Permutation loop
  perm_diffs <- replicate(n_perm, {
    perm_labels <- sample(labels)
    median(combined[perm_labels == "g1"]) - median(combined[perm_labels == "g2"])
  })
  
  # Calculate p-value based on permutation distribution
  p_value <- mean(abs(perm_diffs) >= abs(obs_diff))
  
  # Return results
  if (pairwise) {
    return(data.frame(group_col = group_col,pair1 = timepoint_pair[1],pair2 = timepoint_pair[2],
                      group1 = group1,group2 = group2, obs_diff = round(obs_diff, 3),p_value))
  } else {
    return(data.frame(group_col = group_col,group1 = group1,group2 = group2,obs_diff = round(obs_diff, 3),p_value))
  }
}

#****************
# Load data
#****************
# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_analysis/VIRUSES/")

# Read abundance table and metadata tables
Abundance_table <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_07052025.txt")
Abundance_table_infants <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_INFANTS_07052025.txt")
Abundance_table_mothers <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_MOTHERS_07052025.txt")
Abundance_table_infants_host <- read.delim("Abundance_table/Bacterial_genus_host_vOTU_abundance_INFANTS_17092025.txt")
Abundance_table_mothers_host <- read.delim("Abundance_table/Bacterial_genus_host_vOTU_abundance_MOTHERS_17092025.txt")

Sample_metadata <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt")
Viral_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")

##************************************************************************
# 0. Preprocessing of Sample Metadata
#*************************************************************************

# Set Timepoint and Type variable from metadata as factors
Sample_metadata$Timepoint_categorical <- factor(Sample_metadata$Timepoint_categorical, 
                                                levels=c("P12", "P28", "B", "W2", "M1", "M2", "M3","M6",
                                                         "M9", "M12"))
Sample_metadata$Type <- factor(Sample_metadata$Type)
levels(Sample_metadata$Type) <- c("Infant", "Mother")

# Rename variable with DNA concentration
colnames(Sample_metadata)[colnames(Sample_metadata) == "dna_conc"] <- "DNA_concentration_ng_ul"

# Generate a variable with the proportion of clean reads mapping to viral genomes per samples
Bowtie_mapping_rate <- read.delim("Metadata_NEXT/Bowtie_alignment_rate.txt", sep=" ", header = F)
colnames(Bowtie_mapping_rate) <- c("NG_ID", "Mapping_rate")
Sample_metadata <- left_join(Sample_metadata, Bowtie_mapping_rate, by = "NG_ID")


##************************************************************************
# 1. Alpha Diversity analyses: Richness and Shannon Index vs Time
#*************************************************************************

##################################
# 1.1. Estimation of alpha diversity and richness 
##################################

# Estimate shannon diversity and richness
Sample_metadata$richness  <-specnumber(t(Abundance_table))                 
Sample_metadata$shannon <- vegan::diversity(t(Abundance_table), index = "shannon") 

# Generate the metadata tables for mothers and infants
Sample_metadata_infants <- Sample_metadata[Sample_metadata$Type == "Infant",]
Sample_metadata_mothers <- Sample_metadata[Sample_metadata$Type == "Mother",]

# Plot the richness and diversity in mothers and babies
pdf('2_DIVERSITY/Plots/Shannon_Index_mother_infant.pdf', width=2.7, height=3.2)
ggplot(Sample_metadata[!is.na(Sample_metadata$Timepoint_categorical),], 
       aes(x = Type, y = shannon, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Sample type', y = 'Shannon Index') +
  scale_fill_manual(values = c("#66A6AD", "#8E7CA6")) +
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


pdf('2_DIVERSITY/Plots/Richness_mother_infant.pdf', width=2.9, height=3.2)
ggplot(Sample_metadata[!is.na(Sample_metadata$Timepoint_categorical),], 
       aes(x = Type, y = richness, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Sample type', y = 'Viral Richness') +
  scale_fill_manual(values = c("#66A6AD", "#8E7CA6")) +
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


##################################
# 1.2. Analysis in INFANTS 
##################################

# Set W2 as reference
Sample_metadata_infants$Timepoint_categorical <-
  relevel(as.factor(Sample_metadata_infants$Timepoint_categorical),
          ref = "W2")

# Check statistical significance over time (time as continuous variable)
# A) Richness
MM_time_richness_infant <- lmer(richness ~ exact_age_months_at_collection + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID),
                                REML = F, data = Sample_metadata_infants)
summary(MM_time_richness_infant)   # exact_age_months_at_collection 1.990e+00  8.107e-02 2.721e+03  24.547  < 2e-16 ***

# B) Shannon Index
MM_time_shannon_infant <- lmer(shannon ~ exact_age_months_at_collection + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID),
                               REML = F, data = Sample_metadata_infants)
summary(MM_time_shannon_infant)  # exact_age_months_at_collection 6.970e-02  3.658e-03 2.868e+03  19.056   <2e-16 ***

# Check statistical significance over time (time as categorical)
# A) Richness
MM_time_richness_infant_cat <- lmer(richness ~ Timepoint_categorical + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID),
                                REML = F, data = Sample_metadata_infants)
summary(MM_time_richness_infant_cat)   

# B) Shannon Index
MM_time_shannon_infant_cat <- lmer(shannon ~ Timepoint_categorical + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID),
                               REML = F, data = Sample_metadata_infants)
summary(MM_time_shannon_infant_cat)  


# Generate plots for richness and shannon index 
pdf('2_DIVERSITY/Plots/Shannon_Index_Infants.pdf', width=3.6, height=3.1)
ggplot(Sample_metadata[Sample_metadata$Type == "Infant" & 
                         !is.na(Sample_metadata$Timepoint_categorical) &
                         Sample_metadata$Timepoint_categorical %in% c("W2", "M1", "M2", "M3", "M6",
                                                                      "M9", "M12"), ],
       aes(x = Timepoint_categorical, y = shannon, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Timepoint', y = 'Shannon Index') +
  scale_fill_manual(values = c("Infant" = "#66A6AD")) +
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

pdf('2_DIVERSITY/Plots/Richness_Infants.pdf', width=3.8, height=3.1)
ggplot(Sample_metadata[Sample_metadata$Type == "Infant" & 
                         !is.na(Sample_metadata$Timepoint_categorical) &
                         Sample_metadata$Timepoint_categorical %in% c("W2", "M1", "M2", "M3", "M6",
                                                                      "M9", "M12"), ],
       aes(x = Timepoint_categorical, y = richness, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Timepoint', y = 'Viral Richness') +
  scale_fill_manual(values = c("Infant" = "#66A6AD")) +
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

##################################
# 1.3. Analysis in MOTHERS
##################################

# Check statistical significance over time
Sample_metadata_mothers$Timepoint_categorical <-
  relevel(as.factor(Sample_metadata_mothers$Timepoint_categorical),
          ref = "P12")

# A) Richness
MM_time_richness_mother <- lmer(richness ~ Timepoint_categorical + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID),
                                REML = F, data = Sample_metadata_mothers)
summary(MM_time_richness_mother) 

# B) Shannon Index
MM_time_shannon_mother <- lmer(shannon ~ Timepoint_categorical + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID),
                               REML = F, data = Sample_metadata_mothers)
summary(MM_time_shannon_mother)


# Generate plots for richness and shannon index 
pdf('2_DIVERSITY/Plots/Shannon_Index_Mothers.pdf', width=3.6, height=3.1)
ggplot(Sample_metadata[Sample_metadata$Type == "Mother" & 
                         !is.na(Sample_metadata$Timepoint_categorical) &
                         Sample_metadata$Timepoint_categorical %in% c("P12", "P28", "B", "M3"), ],
       aes(x = Timepoint_categorical, y = shannon, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Timepoint', y = 'Shannon Index') +
  scale_fill_manual(values = c("Mother" = "#7A6395")) +
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

pdf('2_DIVERSITY/Plots/Richness_Mothers.pdf', width=3.8, height=3.1)
ggplot(Sample_metadata[Sample_metadata$Type == "Mother" & 
                         !is.na(Sample_metadata$Timepoint_categorical) &
                         Sample_metadata$Timepoint_categorical %in% c("P12", "P28", "B", "M3"), ],
       aes(x = Timepoint_categorical, y = richness, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Timepoint', y = 'Viral Richness') +
  scale_fill_manual(values = c("Mother" = "#7A6395")) +
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


##************************************************************************
# 2. Beta Diversity analyses: NMDS vs Time
#*************************************************************************

##################################
# 2.1. Overall analysis 
##################################

# Subset abundance table to select only prevalent (>5%) vOTUs
prevalence <- rowSums(Abundance_table > 0) / ncol(Abundance_table) * 100
Abundance_table_5prev <- Abundance_table[prevalence >=5, ]
Abundance_table_5prev_nonzero <- Abundance_table_5prev[, colSums(Abundance_table_5prev) > 0 ]

# Subset metadata to select phenotypes of interest 
NMDS_phenos <- Sample_metadata[Sample_metadata$NG_ID %in% colnames(Abundance_table_5prev_nonzero),]
row.names(NMDS_phenos) <- NMDS_phenos$NG_ID
NMDS_phenos <- NMDS_phenos[,c("Type","DNA_concentration_ng_ul", 
                              "shannon", "Timepoint_categorical")]

# Run NMDS in vegan (metaMDS)
NMS_host_5prev <- metaMDS(t(Abundance_table_5prev_nonzero), distance = "bray", k=2)
en_all<- envfit(NMS_host_5prev, NMDS_phenos, permutations = 1000, na.rm = TRUE)
en_all$factors 

centroids_all <- as.data.frame(scores(en_all, "factors"))
centroids_all$Type <- c(gsub('Type', '', row.names(centroids_all)))
centroids_all$Type[-c(1,2)] <- NA
centroids_all$Timepoint <- c(gsub('Timepoint_categorical', '', row.names(centroids_all)))
centroids_all$Timepoint[c(1,2)] <- NA

data.scores_all = as.data.frame(scores(NMS_host_5prev, "sites"))
data.scores_all <- merge(data.scores_all, NMDS_phenos, by='row.names')
row.names(data.scores_all) <- data.scores_all$Row.names
data.scores_all$Row.names <- NULL

# Calculate centroid and distances for each timepoint
data.scores_all <- data.scores_all %>%
  group_by(Timepoint_categorical) %>%
  mutate(centroid_NMDS1 = mean(NMDS1),
         centroid_NMDS2 = mean(NMDS2),
         distance_from_centroid = sqrt((NMDS1 - centroid_NMDS1)^2 + (NMDS2 - centroid_NMDS2)^2))

# Use IQR method to identify and remove outliers within each timepoint
data.scores_all <- data.scores_all %>%
  group_by(Timepoint_categorical) %>%
  mutate(Q1 = quantile(distance_from_centroid, 0.25),
         Q3 = quantile(distance_from_centroid, 0.75),
         IQR = Q3 - Q1,
         lower_bound = Q1 - 0.05 * IQR,
         upper_bound = Q3 + 0.05 * IQR) %>%
  filter(distance_from_centroid >= lower_bound & distance_from_centroid <= upper_bound) %>%
  ungroup()

# Outliers identified based on beta diversity analysis(based on IQR method):
outliers <- Sample_metadata$bioSampleId[which(!Sample_metadata$clean_reads_FQ_2 %in% data.scores_all$read_depth)]

# Generate NMDS plot with maternal (blue) and infant (green) samples
pdf('2_DIVERSITY/Plots/Viral_vOTUs_ALL_Bray_NMDS_5prev.pdf', width=4.5, height=5.5)
NMDS_plot_all <- ggplot(data = data.scores_all, aes(x = NMDS1, y = NMDS2, color=Type)) + 
  geom_point(size = 1, alpha=0.7) + 
  geom_point(data=centroids_all, aes(fill=Type),shape=NA, size=4, color='black') + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Type, color=Type), linetype = 2, linewidth = 1)+
  theme_bw()+
  scale_color_manual(values = c("Mother" = "#3B4992", "Infant" = "#EE0000")) +
  theme_bw() +
  theme(axis.text=element_text(size=13),
        axis.title = element_text(size = 14), 
        panel.background = element_blank(), 
        panel.grid = element_blank(),  
        panel.border = element_rect(fill = NA, linewidth =1.2, colour = "grey30"), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14, colour = "grey30"),
        legend.position = "bottom") +
  guides(fill = "none") 
ggMarginal(NMDS_plot_all, type="densigram", groupFill=T)
dev.off()

##################################
# 2.2. Analysis in INFANTS 
##################################

# Subset abundance table to select only prevalent (>5%) vOTUs in infants
prevalence <- rowSums(Abundance_table_infants > 0) / ncol(Abundance_table_infants) * 100
Abundance_table_infants_5prev <- Abundance_table_infants[prevalence >= 5, ]
Abundance_table_infants_5prev <- Abundance_table_infants_5prev[, colSums(Abundance_table_infants_5prev) > 0 ]

## NDMS analysis
# Preparing phenotypes
NMDS_phenos_infants <- Sample_metadata_infants [Sample_metadata_infants$NG_ID %in% colnames(Abundance_table_infants_5prev),]
NMDS_phenos_infants <- NMDS_phenos_infants[,c("NG_ID", "read_depth", "DNA_concentration_ng_ul", "Timepoint_categorical")]
row.names(NMDS_phenos_infants) <- NMDS_phenos_infants$NG_ID
NMDS_phenos_infants$NG_ID <- NULL

# Run NMDS in vegan (metaMDS)
NMS_infants <- metaMDS(t(Abundance_table_infants_5prev), distance = "bray", k=2)
en_infants = envfit(NMS_infants, NMDS_phenos_infants, permutations = 1000, na.rm = TRUE)
en_infants$factors 

centroids_infants <- as.data.frame(scores(en_infants, "factors"))
centroids_infants$Timepoint_categorical <- c(gsub('Timepoint_categorical', '', row.names(centroids_infants)))

data.scores_infants <- as.data.frame(scores(NMS_infants, "sites"))
data.scores_infants <- merge(data.scores_infants, NMDS_phenos_infants, by='row.names')
row.names(data.scores_infants) <- data.scores_infants$Row.names
data.scores_infants$Row.names <- NULL

# Calculate centroid and distances for each timepoint
data.scores_infants <- data.scores_infants %>%
  tibble::rownames_to_column(var = "SampleID") %>%
  group_by(Timepoint_categorical) %>%
  mutate(
    centroid_NMDS1 = mean(NMDS1),
    centroid_NMDS2 = mean(NMDS2),
    distance_from_centroid = sqrt((NMDS1 - centroid_NMDS1)^2 + (NMDS2 - centroid_NMDS2)^2)
  ) %>%
  ungroup()

data.scores_infants <- data.scores_infants %>%
  group_by(Timepoint_categorical) %>%
  mutate(centroid_NMDS1 = mean(NMDS1),
         centroid_NMDS2 = mean(NMDS2),
         distance_from_centroid = sqrt((NMDS1 - centroid_NMDS1)^2 + (NMDS2 - centroid_NMDS2)^2))

# Use IQR method to identify and remove outliers within each timepoint
data.scores_infants_filt <- data.scores_infants %>%
  group_by(Timepoint_categorical) %>%
  mutate(Q1 = quantile(distance_from_centroid, 0.25),
         Q3 = quantile(distance_from_centroid, 0.75),
         IQR = Q3 - Q1,
         lower_bound = Q1 - 1.5 * IQR,
         upper_bound = Q3 + 1.5 * IQR) %>%
  filter(distance_from_centroid >= lower_bound & distance_from_centroid <= upper_bound) %>%
  ungroup()

# Generate NMDS plot
pdf('2_DIVERSITY/Plots/Viral_vOTU_INFANTS_Bray_5prev.pdf', width=4.5, height=4.1)
NMDS_plot_infants <- ggplot(data = data.scores_infants_filt, aes(x = NMDS1, y = NMDS2, color=Timepoint_categorical)) + 
  geom_point(size = 1, alpha=0.6) + 
  geom_point(data=centroids_infants, aes(fill=Timepoint_categorical), shape=NA, size=4, color='black') + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Timepoint_categorical, color=Timepoint_categorical),
               linetype = 2, linewidth = 0.8) +
  scale_color_manual(values=c("#1A1F9E", "#3B81B3", "#4CA6B1", "#88CFA4", "#C1D97F", "#E8C26D", "#F4A219", "#D8172A")
  ) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.05))) +
  theme_bw() +
  theme(axis.text=element_text(size=13),
        axis.title = element_text(size = 14), 
        panel.background = element_blank(), 
        panel.grid = element_blank(),  
        panel.border = element_rect(fill = NA, linewidth =0.75, colour = "grey30"), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14, colour = "grey30"),
        legend.position = "bottom") +
  guides(fill = "none") 
ggMarginal(NMDS_plot_infants, type= "boxplot", groupFill=TRUE, outlier.shape = NA)
dev.off()

# Estimate the statistical significance using PERMANOVA
Abundance_table_infants <- Abundance_table_infants[, Sample_metadata_infants$NG_ID] #ensure same order
dist_infants <- vegdist(t(Abundance_table_infants), method = "bray")

PERMANOVA_infants <- adonis2(
  dist_infants ~  Timepoint_categorical + DNA_concentration_ng_ul + read_depth, 
  data = Sample_metadata_infants, 
  permutations = 1000,
  strata = Sample_metadata_infants$NEXT_ID,
  by = "margin"
)

# Due to the high interindividual variability (masks temporal significance), run LMM on NMDS1
data.scores_infants_NMDS <- data.scores_infants[,c("SampleID", "NMDS1")]
colnames(data.scores_infants_NMDS) <- c("NG_ID", "NMDS1")

Sample_metadata_infants <- left_join(Sample_metadata_infants, data.scores_infants_NMDS , by = "NG_ID")

MM_type_composition_infants <- lmer(NMDS1 ~ exact_age_months_at_collection + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID),
                                    REML = F, data = Sample_metadata_infants)
summary(MM_type_composition_infants) 
summary(MM_type_composition_infants)$coefficients["exact_age_months_at_collection", "Pr(>|t|)"]

##################################
# 2.3. Analysis in MOTHERS
##################################

# Filter from the abundance table phages with less than 5% prevalence in mothers
prevalence <- rowSums(Abundance_table_mothers > 0) / ncol(Abundance_table_mothers) * 100
Abundance_table_mothers_5prev <- Abundance_table_mothers[prevalence >= 5, ]
Abundance_table_mothers_5prev_nonzero <- Abundance_table_mothers_5prev[, colSums(Abundance_table_mothers_5prev) > 0 ]

# Calculate dissimilarities in microbiome composition between infant samples
dist_vOTUs_mothers <- vegdist(t(Abundance_table_mothers_5prev_nonzero), index = "bray")

## NDMS analysis
# Preparing phenotypes
NMDS_phenos_mothers <- Sample_metadata_mothers [Sample_metadata_mothers$NG_ID %in% colnames(Abundance_table_mothers_5prev_nonzero),]
NMDS_phenos_mothers <- NMDS_phenos_mothers[,c("NG_ID", "read_depth", "DNA_concentration_ng_ul", "Timepoint_categorical")]
row.names(NMDS_phenos_mothers) <- NMDS_phenos_mothers$NG_ID
NMDS_phenos_mothers$NG_ID <- NULL

# Run NMDS in vegan (metaMDS)
NMS_mothers <- metaMDS(t(Abundance_table_mothers_5prev_nonzero), distance = "bray", k=2)
en_mothers = envfit(NMS_mothers, NMDS_phenos_mothers, permutations = 1000, na.rm = TRUE)
en_mothers$factors 

centroids_mothers <- as.data.frame(scores(en_mothers, "factors"))
centroids_mothers$Timepoint_categorical <- c(gsub('Timepoint_categorical', '', row.names(centroids_mothers)))

data.scores_mothers <- as.data.frame(scores(NMS_mothers, "sites"))
data.scores_mothers <- merge(data.scores_mothers, NMDS_phenos_mothers, by='row.names')
row.names(data.scores_mothers) <- data.scores_mothers$Row.names
data.scores_mothers$Row.names <- NULL

# Calculate centroid and distances for each timepoint
data.scores_mothers <- data.scores_mothers %>%
  tibble::rownames_to_column(var = "SampleID") %>% 
  group_by(Timepoint_categorical) %>%
  mutate(
    centroid_NMDS1 = mean(NMDS1),
    centroid_NMDS2 = mean(NMDS2),
    distance_from_centroid = sqrt((NMDS1 - centroid_NMDS1)^2 + (NMDS2 - centroid_NMDS2)^2)
  ) %>%
  ungroup()

# Use IQR method to identify and remove outliers within each timepoint
data.scores_mothers_filt <- data.scores_mothers %>%
  group_by(Timepoint_categorical) %>%
  mutate(Q1 = quantile(distance_from_centroid, 0.25),
         Q3 = quantile(distance_from_centroid, 0.75),
         IQR = Q3 - Q1,
         lower_bound = Q1 - 1.5 * IQR,
         upper_bound = Q3 + 1.5 * IQR) %>%
  filter(distance_from_centroid >= lower_bound & distance_from_centroid <= upper_bound) %>%
  ungroup()

# Generate NMDS plot
data.scores_mothers_filt$Timepoint_categorical <- factor(data.scores_mothers_filt$Timepoint_categorical, 
                                                levels=c("P12", "P28", "B", "M3"))

pdf('2_DIVERSITY/Plots/Viral_vOTU_MOTHERS_Bray_5prev_NMDS.pdf', width=4.5, height=3.7)
NMDS_plot_mothers <- ggplot(data = data.scores_mothers_filt, aes(x = NMDS1, y = NMDS2, color=Timepoint_categorical)) + 
  geom_point(size = 1, alpha=0.7) + 
  geom_point(data=centroids_mothers, aes(fill=Timepoint_categorical), shape=NA, size=4, color='black') + 
  stat_ellipse(geom = "polygon", alpha = 0.0, aes(group = Timepoint_categorical, color=Timepoint_categorical), linetype = 2, linewidth = 0.8)+
  scale_color_manual(values=c("#F2D7D5", "#E6B0AA", "#D2B4DE", "#A569BD")) +
  scale_y_continuous(expand = expansion(mult = c(0.2, 0.05))) +
  theme_bw()+
  theme(axis.text=element_text(size=13),
        axis.title = element_text(size = 14), 
        panel.background = element_blank(), 
        panel.grid = element_blank(),  
        panel.border = element_rect(fill = NA, linewidth =0.75, colour = "grey30"), 
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.title = element_blank(), 
        legend.text = element_text(size = 14, colour = "grey30"),
        legend.position = "bottom") +
  guides(fill = "none")  
ggMarginal(NMDS_plot_mothers, type= "boxplot", groupFill=T)
dev.off()

# Estimate the statistical significance using PERMANOVA
dist_mothers <- vegdist(t(Abundance_table_mothers), method = "bray")
PERMANOVA_mothers <- adonis2(
  dist_mothers ~ Timepoint_categorical + DNA_concentration_ng_ul + read_depth, 
  data = Sample_metadata_mothers, 
  permutations = 1000,
  strata = Sample_metadata_mothers$NEXT_ID,
  by = "margin"
)

# Associate NMDS1 and NMDS2 with timepoint
data.scores_mothers_NMDS <- data.scores_mothers[,c("SampleID", "NMDS1", "NMDS2")]
colnames(data.scores_mothers_NMDS) <- c("NG_ID", "NMDS1", "NMDS2")

Sample_metadata_mothers <- left_join(Sample_metadata_mothers, data.scores_mothers_NMDS , by = "NG_ID")

MM_type_composition_mothers <- lmer(NMDS2 ~ Timepoint_categorical + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID),
                                    REML = F, data = Sample_metadata_mothers)
summary(MM_type_composition_mothers) 


##************************************************************************
# 3. Acquisition of novel vOTUs
#*************************************************************************

# Estimate the number of novel vOTUs colonizing the infant gut at each timepoint
# All infants are considered for this comparison
# Only comparisons between consecutive timepoints are considered (if not available, NA is added)

# Order timepoints
ordered_timepoints <- c("W2", "M1", "M2", "M3", "M6", "M9", "M12")
Sample_metadata_infants$Timepoint_categorical <- factor(Sample_metadata_infants$Timepoint_categorical,
                                                        levels = ordered_timepoints, ordered = TRUE)

# Sort by NEXT_ID and collection time
Sample_metadata_infants <- Sample_metadata_infants[order(Sample_metadata_infants$NEXT_ID,
                                                         Sample_metadata_infants$exact_age_days_at_collection), ]

# Estimate days since previous sample
# For the first sample of the infant, take exact_age_days_at_collection values
Sample_metadata_infants$Days_since_prev <- ave(
  Sample_metadata_infants$exact_age_days_at_collection,
  Sample_metadata_infants$NEXT_ID,
  FUN = function(x) {
    out <- c(NA, diff(x))
    out[1] <- x[1]  
    out
  }
)
# Convert abundance table to binary
Presence_matrix <- (Abundance_table > 0) * 1

Sample_metadata_infants$Novel_vOTU_count <- NA

# Estimate number of novel vOTUs per sample
for (ind in unique(Sample_metadata_infants$NEXT_ID)) {
  ind_df <- Sample_metadata_infants[Sample_metadata_infants$NEXT_ID == ind, ]
  ind_samples <- ind_df$NG_ID
  ind_timepoints <- as.character(ind_df$Timepoint_categorical)
  
  seen_vOTUs <- c()
  
  for (i in seq_along(ind_samples)) {
    curr_sample <- ind_samples[i]
    curr_time <- ind_timepoints[i]
    
    if (i == 1) {
      # For the 1st sample of the individual, all detected vOTUs are novel
      current_vOTUs <- rownames(Presence_matrix)[Presence_matrix[, curr_sample] == 1]
      Sample_metadata_infants$Novel_vOTU_count[Sample_metadata_infants$NG_ID == curr_sample] <- length(current_vOTUs)
      seen_vOTUs <- dplyr::union(seen_vOTUs, current_vOTUs)
    } else {
      # Check for consecutive timepoints
      prev_time <- ind_timepoints[i - 1]
      prev_idx <- match(prev_time, ordered_timepoints)
      curr_idx <- match(curr_time, ordered_timepoints)
      
      if (!is.na(prev_idx) && !is.na(curr_idx) && curr_idx - prev_idx == 1) {
        current_vOTUs <- rownames(Presence_matrix)[Presence_matrix[, curr_sample] == 1]
        novel_vOTUs <- setdiff(current_vOTUs, seen_vOTUs) # Estimate novel vOTUs not present in previous samples of the individual
        Sample_metadata_infants$Novel_vOTU_count[Sample_metadata_infants$NG_ID == curr_sample] <- length(novel_vOTUs)
        seen_vOTUs <- dplyr::union(seen_vOTUs, current_vOTUs)
      } else {
        # Skip non-consecutive
        Sample_metadata_infants$Novel_vOTU_count[Sample_metadata_infants$NG_ID == curr_sample] <- NA
        seen_vOTUs <- c()  
      }
    }
  }
}

# Estimate novel vOTUs per week
Sample_metadata_infants$Novel_vOTU_per_week <- with(
  Sample_metadata_infants,
  ifelse(!is.na(Days_since_prev) & Days_since_prev > 0,
         Novel_vOTU_count / (Days_since_prev / 7),  
         NA)
)

######################
# Statistical test
######################

# Mixed model for novel vOTUs per week
MM_novel_per_week <- lmer(Novel_vOTU_per_week ~ exact_age_months_at_collection +
                           DNA_concentration_ng_ul + read_depth + (1|NEXT_ID),
                         REML = FALSE, data = Sample_metadata_infants)

summary(MM_novel_per_week) # exact_age_months_at_collection -4.433e-01  2.255e-02  2.696e+03 -19.662  < 2e-16 ***

######################
# Plots
######################

# Novel vOTUs per sample
pdf('2_DIVERSITY/Plots/Novel_vOTU_acquisition_per_sample_trend.pdf', width=5, height=1.5)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$Novel_vOTU_count), ],
       aes(x = exact_age_days_at_collection, y = Novel_vOTU_count,
           color = Timepoint_categorical)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "loess", se = TRUE, color = "black", size = 1) +
  labs(x = 'Age (weeks)', y = 'Novel vOTUs per sample') +
  scale_color_manual(values=c("#5A5F9E", "#3B81B3", "#4CA6B5", "#88CFA4",
                              "#C1D97F", "#E8C26D", "#F4A259")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "right")
dev.off()


pdf('2_DIVERSITY/Plots/Novel_vOTU_acquisition_per_week_trend.pdf', width=5.7, height=3)
ggplot(Sample_metadata_infants[!is.na(Sample_metadata_infants$Novel_vOTU_per_week), ],
       aes(x = exact_age_days_at_collection, y = Novel_vOTU_per_week,
           color = Timepoint_categorical)) +
  geom_point(alpha = 0.6, size = 1.8) +
  geom_smooth(method = "loess", se = TRUE, color = "black", size = 1) +
  labs(x = 'Age (days)', y = 'Novel vOTUs per week') +
  scale_color_manual(values=c("#5A5F9E", "#3B81B3", "#4CA6B5", "#88CFA4",
                              "#C1D97F", "#E8C26D", "#F4A259")) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, size = 0.75),
        axis.text = element_text(size = 14),
        axis.title = element_text(size = 16),
        axis.title.y = element_text(margin = margin(r = 10)),
        axis.title.x = element_text(margin = margin(t = 10)),
        legend.position = "right",
        legend.text  = element_text(size = 12),   
        legend.title = element_text(size = 13))
dev.off()


##************************************************************************
# 4. Within-individual (Bray) distance comparison
#*************************************************************************

# Comparison of within-individual distances in mothers and infants over time
# First, calculate dissimilarities in microbiome composition between infant samples
dist_vOTUs_mothers <- vegdist(t(Abundance_table_mothers), index = "bray")
dist_vOTUs_infants <- vegdist(t(Abundance_table_infants), index = "bray")

# Transform distance matrix into a long format
dist_mat_infants <- as.matrix(dist_vOTUs_infants)
dist_mat_infants[lower.tri(dist_mat_infants, diag = TRUE)] <- NA  
distances_infants <- reshape2::melt(dist_mat_infants, na.rm = TRUE)
colnames(distances_infants) <- c("Sample1", "Sample2", "Distance")

# Add metadata to the long-formatted table
distances_infants$Timepoint1 <- Sample_metadata_infants$Timepoint_categorical[match(distances_infants$Sample1, Sample_metadata_infants$NG_ID)]
distances_infants$Timepoint2 <-  Sample_metadata_infants$Timepoint_categorical[match(distances_infants$Sample2, Sample_metadata_infants$NG_ID)]
distances_infants$NEXT_ID1 <- Sample_metadata_infants$NEXT_ID[match(distances_infants$Sample1, Sample_metadata_infants$NG_ID)]
distances_infants$NEXT_ID2 <- Sample_metadata_infants$NEXT_ID[match(distances_infants$Sample2, Sample_metadata_infants$NG_ID)]
distances_infants$Type <- rep("Infant", nrow(distances_infants))
distances_infants$Within_between <- ifelse(distances_infants$NEXT_ID1 == distances_infants$NEXT_ID2, "Within", "Between")

# Repeat for mothers:Transform distance matrix into a long format and add metadata
dist_mat_mothers <- as.matrix(dist_vOTUs_mothers)
dist_mat_mothers[lower.tri(dist_mat_mothers, diag = TRUE)] <- NA  # Keep upper triangle only
distances_mothers <- reshape2::melt(dist_mat_mothers, na.rm = TRUE)
colnames(distances_mothers) <- c("Sample1", "Sample2", "Distance")
distances_mothers$Timepoint1 <- Sample_metadata_mothers$Timepoint_categorical[match(distances_mothers$Sample1, Sample_metadata_mothers$NG_ID)]
distances_mothers$Timepoint2 <-  Sample_metadata_mothers$Timepoint_categorical[match(distances_mothers$Sample2, Sample_metadata_mothers$NG_ID)]
distances_mothers$NEXT_ID1 <- Sample_metadata_mothers$NEXT_ID[match(distances_mothers$Sample1, Sample_metadata_mothers$NG_ID)]
distances_mothers$NEXT_ID2 <- Sample_metadata_mothers$NEXT_ID[match(distances_mothers$Sample2, Sample_metadata_mothers$NG_ID)]
distances_mothers$Type <- rep("Mother", nrow(distances_mothers))
distances_mothers$Within_between <- ifelse(distances_mothers$NEXT_ID1 == distances_mothers$NEXT_ID2, "Within", "Between")

# Filter within-individual distances only
within_infants <- subset(distances_infants, Within_between == "Within")
within_mothers <- subset(distances_mothers, Within_between == "Within")

# Create a combined timepoint variable (where e.g., “W2-M1” and “M1-W2” become the same)
within_infants$Timepoint_pair <- apply(within_infants[, c("Timepoint1", "Timepoint2")], 1, function(x) paste(sort(x), collapse = "-"))
within_mothers$Timepoint_pair <- apply(within_mothers[, c("Timepoint1", "Timepoint2")], 1, function(x) paste(sort(x), collapse = "-"))

# Combine datasets
combined_within <- rbind(within_infants, within_mothers)

# Check median and 1st and 3rd quartiles of maternal and infants distances 
# Use similarity (1 - Distance) (as used in the plots)
combined_within %>%
  mutate(Similarity = 1 - Distance) %>%
  group_by(Type) %>%
  summarise(
    Median = median(Similarity, na.rm = TRUE),
    Q1 = quantile(Similarity, 0.25, na.rm = TRUE),
    Q3 = quantile(Similarity, 0.75, na.rm = TRUE)
  )

######################
# Statistical test
######################

# Test if maternal within-individual distances are smaller than the initial distances
overall_result <- perm_test(group1 = "Mother", group2 = "Infant", group_col = "Type",
                            data = combined_within, pairwise = FALSE)

# Test if infant M1-W2 distances are different or not to every maternal comparison distances
comparisons <- list(c("B-P12", "M3-P12"), c("P12-P28", "M3-P12"), c("P12-P28", "B-P12"))

results_perm_distances <- do.call(rbind, lapply(comparisons, function(x) {
  perm_test(group1 = "Mother", group2 = "Mother", group_col = "Type",
            data = combined_within, 
            timepoint_pair = x, pairwise = TRUE, n_perm = 10000)
}))


######################
# Plots
######################

# For the plot, represent similarity instead of Distance (1-Distance)

# A) Overall distances mothers vs infants
pdf('2_DIVERSITY/Plots/Within_individual_similarity_mother_infant.pdf', width=2.9, height=3.2)
ggplot(combined_within, aes(x = Type, y = 1 - Distance, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Sample type', y = 'Within-individual similarity') +
  scale_fill_manual(values = c("Infant" = "#66A6AD", "Mother" = "#8E7CA6")) +
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

# B) Infant distances W2-M1 vs maternal distances at each timepoint
# Select only timepoint comparisons vs 1st timepoint
ordered_pairs <- c("P12-P28", "B-P12", "M3-P12",      
                   "M1-W2", "M2-W2", "M3-W2",      
                   "M6-W2", "M9-W2", "M12-W2")      

# Filter the combined dataset
filtered_combined <- combined_within[combined_within$Timepoint_pair %in% ordered_pairs, ]

# Reorder the Timepoint_pair factor
filtered_combined$Timepoint_pair <- factor(filtered_combined$Timepoint_pair, levels = ordered_pairs)

pdf('2_DIVERSITY/Plots/Within_individual_similarity_timepoint_pairs.pdf', width=3.5, height=3.8)
ggplot(filtered_combined, aes(x = Timepoint_pair, y = 1 - Distance, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  scale_fill_manual(values = c("Infant" = "#66A6AD", "Mother" = "#8E7CA6")) +
  labs(x = "Timepoint pair", y = "Within-individual similarity") +
  scale_y_continuous(
    limits = c(0, 1),               
    breaks = seq(0, 1, 0.2),       
    expand = expansion(mult = c(0.1, 0.25))  
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 12),
    axis.text.y = element_text(size = 12),
    axis.title = element_text(size = 14),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "top",
    legend.text = element_text(size = 12)
  )
dev.off()

##************************************************************************
# 4.1. Between-individual (Bray) distance comparison (within also included)
#*************************************************************************

# Filter between-individual distances only
between_infants <- subset(distances_infants, Within_between == "Between")
between_mothers <- subset(distances_mothers, Within_between == "Between")

# Create a combined timepoint variable (where e.g., “W2-M1” and “M1-W2” become the same)
between_infants$Timepoint_pair <- apply(between_infants[, c("Timepoint1", "Timepoint2")], 1, function(x) paste(sort(x), collapse = "-"))
between_mothers$Timepoint_pair <- apply(between_mothers[, c("Timepoint1", "Timepoint2")], 1, function(x) paste(sort(x), collapse = "-"))

# Combine with within-individual distances (include them for the plot)
all_infants <- rbind(between_infants, within_infants)
all_mothers <- rbind(between_mothers, within_mothers)

# Add relatedness information (related)
all_infants$FAM_ID1 <- Sample_metadata_infants$FAMILY[match(all_infants$NEXT_ID1, Sample_metadata_infants$NEXT_ID)]
all_infants$FAM_ID2 <- Sample_metadata_infants$FAMILY[match(all_infants$NEXT_ID2, Sample_metadata_infants$NEXT_ID)]
all_mothers$FAM_ID1 <- Sample_metadata_mothers$FAMILY[match(all_mothers$NEXT_ID1, Sample_metadata_mothers$NEXT_ID)]
all_mothers$FAM_ID2 <- Sample_metadata_mothers$FAMILY[match(all_mothers$NEXT_ID2, Sample_metadata_mothers$NEXT_ID)]

all_infants$Relatedness <- ifelse(all_infants$Within_between == "Within", "Same",
  ifelse(all_infants$Within_between == "Between" & all_infants$FAM_ID1 == all_infants$FAM_ID2, "Related",
    "Unrelated"))
all_mothers$Relatedness <- ifelse(all_mothers$Within_between == "Within", "Same pregnancy",
                                  ifelse(all_mothers$Within_between == "Between" & all_mothers$FAM_ID1 == all_mothers$FAM_ID2, "Diff pregnancy",
                                         "Unrelated"))

######################
# Infant plots
######################

# Select only timepoint comparisons vs first timepoint
ordered_pairs <- c("P12-P28", "B-P12", "M3-P12",      
                   "M1-W2", "M2-W2", "M3-W2",      
                   "M6-W2", "M9-W2", "M12-W2") 

# Filter the dataset
filtered_all_infants <- all_infants[all_infants$Timepoint_pair %in% ordered_pairs, ]

# Reorder the Timepoint_pair and Relatedness factors
filtered_all_infants$Timepoint_pair <- factor(filtered_all_infants$Timepoint_pair, levels = ordered_pairs)
filtered_all_infants$Relatedness <- factor(filtered_all_infants$Relatedness,
                                           levels = c("Same", "Related", "Unrelated"))

# Subsample unrelated distances for plotting
same_diff <- dplyr::filter(filtered_all_infants, Relatedness != "Unrelated")
unrelated <- dplyr::filter(filtered_all_infants, Relatedness == "Unrelated")
unrelated_sub <- dplyr::slice_sample(unrelated, n = min(20000, nrow(unrelated)))
filtered_all_infants_sub <- dplyr::bind_rows(same_diff, unrelated_sub)

pdf('2_DIVERSITY/Plots/Between_individual_similarity_timepoint_pairs_INFANTS.pdf', width=7, height=3.8)
ggplot(filtered_all_infants_sub, aes(x = Timepoint_pair, y = 1 - Distance, fill = Relatedness)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  geom_boxplot(width = 0.8, color = "black", outlier.color = NA, size = 0.9) +
  scale_fill_manual(values = c("Same" = "#A1B5D8" , "Related" = "#3D405B", "Unrelated" ="#E07A5F"  )) +
  labs(x = "Timepoint pair", y = "Between-individual similarity") +
  scale_y_continuous(
    limits = c(0, 1),               
    breaks = seq(0, 1, 0.2),       
    expand = expansion(mult = c(0.1, 0.1))  
  ) +
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
    legend.position = "top",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
dev.off()

# Visualize also as a lineplot
filtered_all_infants_summ <- filtered_all_infants %>%
  group_by(Timepoint_pair, Relatedness) %>%
  summarize(
    median_similarity = median(1 - Distance),
    q1 = quantile(1 - Distance, 0.25),
    q3 = quantile(1 - Distance, 0.75),
    .groups = "drop"
  )

pdf('2_DIVERSITY/Plots/Between_individual_similarity_timepoint_pairs_INFANTS_lineplot.pdf', width = 5, height = 3.8)
ggplot(filtered_all_infants_summ, 
       aes(x = Timepoint_pair, y = median_similarity, group = Relatedness, color = Relatedness)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = Relatedness), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("Same" = "#A1B5D8" , "Related" = "#3D405B", "Unrelated" ="#E07A5F")) +
  scale_fill_manual(values = c("Same" = "#A1B5D8" , "Related" = "#3D405B", "Unrelated" ="#E07A5F")) +
  labs(x = "Timepoint pair", y = "Between-individual similarity") +
  scale_y_continuous(limits = c(0, 0.55), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(expand = expansion(add = c(0.4, 0.4))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text.x = element_text(angle=45, hjust=1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "top",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16))
dev.off()

######################
# Maternal plots
######################

# Filter the dataset
filtered_all_mothers <- all_mothers[all_mothers$Timepoint_pair %in% ordered_pairs, ]

# Reorder the Timepoint_pair and Relatedness factors
filtered_all_mothers$Timepoint_pair <- factor(filtered_all_mothers$Timepoint_pair, levels = ordered_pairs)
filtered_all_mothers$Relatedness <- factor(filtered_all_mothers$Relatedness,
                                           levels = c("Same pregnancy", "Diff pregnancy", "Unrelated"))

# Subsample unrelated distances for plotting
same_diff <- dplyr::filter(filtered_all_mothers, Relatedness != "Unrelated")
unrelated <- dplyr::filter(filtered_all_mothers, Relatedness == "Unrelated")
unrelated_sub <- dplyr::slice_sample(unrelated, n = min(20000, nrow(unrelated)))
filtered_all_mothers_sub <- dplyr::bind_rows(same_diff, unrelated_sub)

pdf('2_DIVERSITY/Plots/Between_individual_similarity_timepoint_pairs_mothers.pdf', width=5.1, height=3.8)
ggplot(filtered_all_mothers_sub, aes(x = Timepoint_pair, y = 1 - Distance, fill = Relatedness)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8,
              position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.8)) +
  geom_boxplot(width = 0.8, color = "black", outlier.color = NA, size = 0.9) +
  scale_fill_manual(values = c("Same pregnancy" = "#A1B5D8" , "Diff pregnancy" = "#3D405B", "Unrelated" ="#E07A5F")) +
  labs(x = "Timepoint pair", y = "Between-individual similarity") +
  scale_y_continuous(
    limits = c(0, 1),               
    breaks = seq(0, 1, 0.2),       
    expand = expansion(mult = c(0.1, 0.1))  
  ) +
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
    legend.position = "top",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16)
  )
dev.off()

# Visualize also as a lineplot
filtered_all_mothers_summ <- filtered_all_mothers %>%
  group_by(Timepoint_pair, Relatedness) %>%
  summarize(
    median_similarity = median(1 - Distance),
    q1 = quantile(1 - Distance, 0.25),
    q3 = quantile(1 - Distance, 0.75),
    .groups = "drop"
  )

pdf('2_DIVERSITY/Plots/Between_individual_similarity_timepoint_pairs_mothers_lineplot.pdf', width = 5, height = 3.8)
ggplot(filtered_all_mothers_summ, 
       aes(x = Timepoint_pair, y = median_similarity, group = Relatedness, color = Relatedness)) +
  geom_line(size = 1) +
  geom_point(size = 3) +
  geom_ribbon(aes(ymin = q1, ymax = q3, fill = Relatedness), alpha = 0.2, color = NA) +
  scale_color_manual(values = c("Same pregnancy" = "#A1B5D8" , "Diff pregnancy" = "#3D405B", "Unrelated" ="#E07A5F" )) +
  scale_fill_manual(values = c("Same pregnancy" = "#A1B5D8" , "Diff pregnancy" = "#3D405B", "Unrelated" ="#E07A5F")) +
  labs(x = "Timepoint pair", y = "Between-individual similarity") +
  scale_y_continuous(limits = c(0, 0.55), breaks = seq(0, 1, 0.2), expand = expansion(mult = c(0.05, 0.1))) +
  scale_x_discrete(expand = expansion(add = c(0.4, 0.4))) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text.x = element_text(angle=45, hjust=1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(size = 16),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.title.x = element_text(margin = margin(t = 10)),
    legend.position = "top",
    legend.text = element_text(size = 14),
    legend.title = element_text(size = 16))
dev.off()

######################
# Statistical test
######################

# Perform statistical test to check if:
# A) Related distances in infants are lower than unrelated distances
result_infant_comparison <- perm_test("Same", "Unrelated", data = filtered_all_infants,
                  group_col = "Relatedness", pairwise = FALSE, n_perm = 10000)
result_infant_comparison2 <- perm_test("Related", "Unrelated", data = filtered_all_infants,
                                      group_col = "Relatedness", pairwise = FALSE, n_perm = 10000)

# B) Related distances in mothers are lower than unrelated distances
result_mother_comparison <- perm_test("Same pregnancy", "Unrelated", data = filtered_all_mothers,
                                      group_col = "Relatedness", pairwise = FALSE, n_perm = 10000)
result_mother_comparison2 <- perm_test("Diff pregnancy", "Unrelated", data = filtered_all_mothers,
                                       group_col = "Relatedness", pairwise = FALSE, n_perm = 10000)



##************************************************************************
# 5. Virome prevalence and specificity 
##************************************************************************

##################################
# 5.1. vOTU Prevalence 
##################################

# Estimate the number of vOTUs present in mothers and infants at different prevalence thresholds
N_samples_mothers <- ncol(Abundance_table_mothers)
N_samples_infants <- ncol(Abundance_table_infants)
Viral_presence_mothers <- rowSums(Abundance_table_mothers > 0)
Viral_presence_infants <- rowSums(Abundance_table_infants > 0)
Prevalence_mothers <- 100*(Viral_presence_mothers / N_samples_mothers)
Prevalence_infants <- 100*(Viral_presence_infants / N_samples_infants)
prevalence_thresholds <- c(5)

vOTUs_above_thresholds <- list()
for (threshold in prevalence_thresholds) {
  vOTUs_mothers <- names(Prevalence_mothers[Prevalence_mothers >= threshold])
  vOTUs_infants <- names(Prevalence_infants[Prevalence_infants >= threshold])
  cat(paste0("Threshold ≥ ", threshold, "%\n"))
  cat("  Mothers: ", length(vOTUs_mothers), " viruses\n")
  cat("  Infants: ", length(vOTUs_infants), " viruses\n\n")
  vOTUs_above_thresholds[[paste0("Mothers_", threshold, "%")]] <- vOTUs_mothers
  vOTUs_above_thresholds[[paste0("Infants_", threshold, "%")]] <- vOTUs_infants
}

# Get IDs of vOTUs with ≥ prevalence in infants
prevalent_vOTUs_infants <- vOTUs_above_thresholds$`Infants_5%`
Viral_metadata_infants[Viral_metadata_infants$Virus_ID %in% prevalent_vOTUs_infants, "Class"]

# Compute a full prevalence data frame for all thresholds
thresholds <- seq(0, 100, by =0.1)   

# Calculate number of vOTUs above each threshold
prevalence_vOTUs <- expand.grid(
  Threshold = thresholds,
  Type = c("Infant", "Mother")
) %>%
  dplyr::mutate(
    N_vOTUs = ifelse(
      Type == "Infant",
      sapply(Threshold, function(t) sum(Prevalence_infants >= t)),
      sapply(Threshold, function(t) sum(Prevalence_mothers >= t))
    )
  )

max_prev_infant <- max(prevalence_vOTUs$Threshold[prevalence_vOTUs$Type == "Infant" & prevalence_vOTUs$N_vOTUs > 0])
max_prev_mother <- max(prevalence_vOTUs$Threshold[prevalence_vOTUs$Type == "Mother" & prevalence_vOTUs$N_vOTUs > 0])
max_prev <- max(max_prev_infant, max_prev_mother)


# Perform statistical test
# Wilcoxon test: treating the number of vOTUs at each threshold as paired and comparing groups
wide_data_prevalence <- prevalence_vOTUs%>%
  pivot_wider(names_from = Type, values_from = N_vOTUs)

wilcox.test(wide_data_prevalence$Infant, wide_data_prevalence$Mother, paired = TRUE) #p=8.511903e-72

# Generate the plot
pdf('2_DIVERSITY/Plots/vOTU_prevalence_curve_filled.pdf', width=3.8, height=3.8)
ggplot(prevalence_vOTUs, aes(x = Threshold, y = N_vOTUs, fill = Type)) +
  geom_area(alpha = 0.3, position = "identity") +     
  geom_line(aes(color = Type), size = 1.2) +          
  geom_point(aes(color = Type), size = 1, alpha = 0.3) +  
  scale_y_log10() +
  scale_x_continuous(limits = c(0, max_prev + 1)) +   
  scale_color_manual(values = c("Infant" = "#66A6AD", "Mother" = "#8E7CA6")) +
  scale_fill_manual(values = c("Infant" = "#66A6AD", "Mother" = "#8E7CA6")) +
  labs(x = "Prevalence threshold (%)", y = "vOTUs above threshold (log10)") +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),  
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    legend.position = "top",
    legend.title = element_blank(),
    legend.text = element_text(size = 14)
  )
dev.off()

##################################
# 5.2. Individual specificity
##################################

# Estimate number of viruses present in single/multiple individuals
Abundance_mat <- as.matrix(Abundance_table)
ng_to_next <- setNames(Sample_metadata$NEXT_ID, Sample_metadata$NG_ID)

# Function to estimate how many unique NEXT_IDs per vOTU
get_unique_next_ids <- function(votu_abundances) {
  present_samples <- names(votu_abundances[votu_abundances != 0])
  present_next_ids <- ng_to_next[present_samples]
  length(unique(present_next_ids))
}

votu_next_counts <- apply(Abundance_mat, 1, get_unique_next_ids)
vOTUs_multiple_ind <- names(votu_next_counts[votu_next_counts >= 2])
vOTUs_single_ind   <- names(votu_next_counts[votu_next_counts == 1]) #17,630


#****************
# Write results
#****************
# Sample metadata
write.table(Sample_metadata,"Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Sample_metadata_infants,"Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Sample_metadata_mothers,"Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 

