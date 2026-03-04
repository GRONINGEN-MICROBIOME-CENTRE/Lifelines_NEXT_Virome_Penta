################################################################################
##### LL-NEXT: Viral Lifestyle Analysis
### Author(s):Asier Fernández-Pato
### Last updated: 22nd December, 2025
################################################################################

#****************
# Load modules
#****************
library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)

# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_analysis/VIRUSES/")

##************************************************************************
# 0. Load metadata and abundance table for the LL-NEXT samples 
#*************************************************************************
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
                                         
# Set timepoint variable from metadata as factor
Sample_metadata$Timepoint_categorical <- factor(Sample_metadata$Timepoint_categorical, 
                                                levels=c("P12", "P28", "B", "W2", "M1", "M2", "M3","M6",
                                                         "M9", "M12"))


##************************************************************************
  # 1. Viral lifestyle: estimation of richness, diversity and abundance of temperate phages
##************************************************************************

# Estimate the richness (number) of temperate phages per sample
temperate_phages <- Viral_metadata$Virus_ID[which(Viral_metadata$Lifestyle == "Temperate")]
Sample_metadata$richness_temperate_phages <- colSums(Abundance_table[rownames(Abundance_table) %in% temperate_phages, ] > 0)

# Estimate the proportion and relative abundance of temperate phages per sample
# Note that this calculation includes phages without lifestyle assigned
Sample_metadata$relab_temperate_phages <- 100*(colSums(Abundance_table[rownames(Abundance_table) %in% temperate_phages,]) / 
                                    colSums(Abundance_table)) 
Sample_metadata$proportion_temperate_phages <- 100*(Sample_metadata$richness_temperate_phages / Sample_metadata$richness)

# Add the estimated variables to the infant metadata and maternal metadata
lifestyle_variables <- Sample_metadata[,c("NG_ID", "richness_temperate_phages", "relab_temperate_phages",
                                          "proportion_temperate_phages")]

Sample_metadata_infants <- left_join(Sample_metadata_infants, lifestyle_variables, by = "NG_ID")
Sample_metadata_mothers <- left_join(Sample_metadata_mothers, lifestyle_variables, by = "NG_ID")


##****************************************
# 2. Viral lifestyle: mothers vs infants
##****************************************

# Test significance: mothers vs infants
# A) Relative abundance temperate phages
MM_type_relab_temp_phages <- lmer(relab_temperate_phages ~ Type + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID), REML = F, data = Sample_metadata)
summary(MM_type_relab_temp_phages)

# B) Proportion temperate phages
MM_type_proportion_temp_phages <- lmer(proportion_temperate_phages ~ Type + DNA_concentration_ng_ul + read_depth + (1|NEXT_ID), REML = F, data = Sample_metadata)
summary(MM_type_proportion_temp_phages)


# Plot the proportion and relative abundances in mothers and babies
pdf('3_LIFESTYLE/Plots/Prop_temperate_phages_mother_infant.pdf', width = 2.9, height = 3.2)
ggplot(Sample_metadata[!is.na(Sample_metadata$Timepoint_categorical),], 
       aes(x = Type, y = proportion_temperate_phages, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Sample type', y = 'Proportion temperate phages') +
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

pdf('3_LIFESTYLE/Plots/Relab_temperate_phages_mother_infant.pdf', width = 2.9, height = 3.2)
ggplot(Sample_metadata[!is.na(Sample_metadata$Timepoint_categorical),], 
       aes(x = Type, y = relab_temperate_phages, fill = Type)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_boxplot(width = 0.6, color = "black", outlier.color = NA, size = 0.9) +
  labs(x = 'Sample type', y = 'Relative abundance temperate phages') +
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


##***********************************************
# 2. Viral lifestyle analysis: Infants over time
#************************************************

# Test significance: infants over time

# Set W2 as reference
Sample_metadata_infants$Timepoint_categorical <-
  relevel(as.factor(Sample_metadata_infants$Timepoint_categorical),
          ref = "W2")

# A) Relative abundance temperate phages
MM_time_relab_temp_phages_infants <- lmer(relab_temperate_phages ~ exact_age_months_at_collection + 
                                            DNA_concentration_ng_ul + read_depth + (1|NEXT_ID), REML = F, data = Sample_metadata_infants)
summary(MM_time_relab_temp_phages_infants)

MM_time_relab_temp_phages_infants_cat <- lmer(relab_temperate_phages ~ Timepoint_categorical + 
                                            DNA_concentration_ng_ul + read_depth + (1|NEXT_ID), REML = F, data = Sample_metadata_infants)
summary(MM_time_relab_temp_phages_infants_cat)

# B) Proportion temperate phages
MM_time_proportion_temp_phages_infants <- lmer(proportion_temperate_phages ~ exact_age_months_at_collection + 
                                                 DNA_concentration_ng_ul + read_depth + (1|NEXT_ID), REML = F, data = Sample_metadata_infants)
summary(MM_time_proportion_temp_phages_infants)

MM_time_proportion_temp_phages_infants_cat <- lmer(proportion_temperate_phages ~ Timepoint_categorical + 
                                                 DNA_concentration_ng_ul + read_depth + (1|NEXT_ID), REML = F, data = Sample_metadata_infants)
summary(MM_time_proportion_temp_phages_infants_cat)

# Generate plots for the proportion and relative abundance of temperate phages
pdf('3_LIFESTYLE/Plots/Proportion_temperate_phages_infants.pdf', width=3.8, height=3.1)
ggplot(Sample_metadata[Sample_metadata$Type=="Infant" & !is.na(Sample_metadata$Timepoint_categorical) &
      Sample_metadata$Timepoint_categorical %in% c("W2","M1","M2","M3","M6","M9","M12"),],
  aes(x=Timepoint_categorical, y=proportion_temperate_phages, fill=Type)) +
  geom_jitter(alpha=0.4, color="lightgrey", size=1.8, width=0.2) +
  geom_boxplot(width=0.6, color="black", outlier.color=NA, size=0.9) +
  labs(x='Timepoint', y='Proportion temperate phages') +
  scale_fill_manual(values=c("#66A6AD")) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", fill=NA, size=0.75),
    axis.text = element_text(size=14),
    axis.title = element_text(size=16),
    axis.title.y = element_text(margin = margin(r=10)),
    axis.title.x = element_text(margin = margin(t=10)),
    legend.position="none"
  )
dev.off()

pdf('3_LIFESTYLE/Plots/Relab_temperate_phages_infants.pdf', width=3.8, height=3.1)
ggplot(Sample_metadata[Sample_metadata$Type=="Infant" &!is.na(Sample_metadata$Timepoint_categorical) &
      Sample_metadata$Timepoint_categorical %in% c("W2","M1","M2","M3","M6","M9","M12"),],
  aes(x=Timepoint_categorical, y=relab_temperate_phages, fill=Type)) +
  geom_jitter(alpha=0.4, color="lightgrey", size=1.8, width=0.2) +
  geom_boxplot(width=0.6, color="black", outlier.color=NA, size=0.9) +
  labs(x='Timepoint', y='Rel.abundance temperate phages') +
  scale_fill_manual(values=c("#66A6AD")) +
  scale_y_continuous(
    limits = c(0, 100),
    expand = expansion(mult = c(0.1, 0.3)),   
    breaks = seq(0, 100, by = 20)             
  ) +
  theme_bw() +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color="black", fill=NA, size=0.75),
    axis.text = element_text(size=14),
    axis.title = element_text(size=16),
    axis.title.y = element_text(margin = margin(r=10)),
    axis.title.x = element_text(margin = margin(t=10)),
    legend.position="none"
  )
dev.off()


##***********************************************
# 2. Viral lifestyle analysis: mothers over time
#************************************************

# Test significance: mothers over time
Sample_metadata_mothers$Timepoint_categorical <-
  relevel(as.factor(Sample_metadata_mothers$Timepoint_categorical),
          ref = "P12")

# A) Relative abundance temperate phages
MM_time_relab_temp_phages_mothers <- lmer(relab_temperate_phages ~ exact_age_months_at_collection + 
                                            DNA_concentration_ng_ul + read_depth + (1|NEXT_ID), REML = F, data = Sample_metadata_mothers)
summary(MM_time_relab_temp_phages_mothers)

MM_time_relab_temp_phages_mothers_cat <- lmer(relab_temperate_phages ~ Timepoint_categorical + 
                                            DNA_concentration_ng_ul + read_depth + (1|NEXT_ID), REML = F, data = Sample_metadata_mothers)
summary(MM_time_relab_temp_phages_mothers_cat)

# B) Proportion temperate phages
MM_time_proportion_temp_phages_mothers <- lmer(proportion_temperate_phages ~ exact_age_months_at_collection + 
                                                 DNA_concentration_ng_ul + read_depth + (1|NEXT_ID), REML = F, data = Sample_metadata_mothers)
summary(MM_time_proportion_temp_phages_mothers)

MM_time_proportion_temp_phages_mothers_cat <- lmer(proportion_temperate_phages ~ Timepoint_categorical + 
                                                 DNA_concentration_ng_ul + read_depth + (1|NEXT_ID), REML = F, data = Sample_metadata_mothers)
summary(MM_time_proportion_temp_phages_mothers_cat)


##*************
# Save output
#**************
write.table(Sample_metadata,"Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Sample_metadata_infants,"Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Sample_metadata_mothers,"Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 

