################################################################################
##### LL-NEXT: Associations between vOTU abundances and time
### Author(s):Asier Fernández-Pato
### Last updated: 22nd December, 2025
################################################################################

#****************
# Load modules
#****************
library(lme4)
library(lmerTest)
library(tidyverse)
library(dplyr)
library(ggplot2)
library(patchwork)


#****************
# Define functions
#****************

# Function to test associations between vOTU abundances and phenotypes using mixed-effects models
associate_vOTUs_with_phenotypes <- function(sample_metadata, subject_id_col, abundance_data, phenotype_list) {
  
  # Prepare metadata with subject IDs as rownames
  df <- sample_metadata
  row.names(df) <- df[[subject_id_col]]
  
  # Merge metadata with abundance data
  df <- merge(df, abundance_data, by = "row.names")
  row.names(df) <- df$Row.names
  df$Row.names <- NULL
  
  # vOTUs (features) to test
  vOTUs <- colnames(abundance_data)
  results_all <- tibble()
  
  # Loop over each vOTU
  for (i in seq_along(vOTUs)) {
    vOTU <- vOTUs[i]
    cat("vOTU", i, "/", length(vOTUs), "\n")
    if (!vOTU %in% colnames(df)) next
    vOTU_safe <- paste0("`", vOTU, "`")
    
    # Test association with each phenotype
    for (phenotype in phenotype_list) {
      phenotype_safe <- paste0("`", phenotype, "`")
      
      # Keep only rows with non-missing phenotype
      valid_samples <- df[!is.na(df[[phenotype]]), subject_id_col]
      df_sub <- filter(df, !!sym(subject_id_col) %in% valid_samples)
      
      # Null model (no phenotype)
      model_null <- lmer(as.formula(paste(vOTU_safe, "~ read_depth + DNA_concentration_ng_ul + (1|NEXT_ID)")),
                         df_sub, REML = FALSE)
      
      # Full model (includes phenotype)
      model_full <- lmer(as.formula(paste(vOTU_safe, "~ read_depth + DNA_concentration_ng_ul +", phenotype_safe, "+ (1|NEXT_ID)")),
                         df_sub, REML = FALSE)
      
      # Likelihood ratio test
      lrt <- anova(model_full, model_null)
      p_val <- lrt$`Pr(>Chisq)`[2] 
      
      # Extract phenotype effect estimates
      coef_info <- summary(model_full)$coefficients
      coef_table <- as.data.frame(coef_info)[grep(phenotype, rownames(coef_info)), ]
      
      # Format and store results
      temp_result <- coef_table %>%
        rownames_to_column("Term") %>%
        as_tibble() %>%
        mutate(
          P = p_val,
          Model = "Mixed",
          vOTU = vOTU,
          Phenotype = phenotype
        )
      
      results_all <- bind_rows(results_all, temp_result)
    }
  }
  
  # Multiple testing correction
  results_all <- results_all %>%
    mutate(FDR = p.adjust(P, method = "BH"))
  
  return(results_all)
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
Sample_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_07052025.txt")
Sample_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_07052025.txt")

Virus_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")
Virus_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_INFANTS_07052025.txt")
Virus_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_MOTHERS_07052025.txt")


##************************************************************************
# 1. Abundance table processing
##************************************************************************

################
# INFANTS - vOTU level
################

# Perform a log transformation and filter by prevalence
Abundance_table_infants <- t(Abundance_table_infants)

# Calculate the prevalence of each vOTU
vOTU_nonZero  <- colSums(Abundance_table_infants>0) / nrow(Abundance_table_infants) 

# Get the list of vOTUs present in more than 5% samples
vOTU_keep <- colnames(Abundance_table_infants)[vOTU_nonZero > 0.05]

# Generate LOG transformed abundance tables
## CLR transformation not suitable here: data hardly compositional after read mapping
## Log transformation can help with data normality necessary for linear models
pseudocount_vOTU <- min(Abundance_table_infants[Abundance_table_infants > 0]) /2
Abundance_table_infants_log <- Abundance_table_infants
Abundance_table_infants_log[Abundance_table_infants_log  == 0] <- pseudocount_vOTU
Abundance_table_infants_log <- log(Abundance_table_infants_log)

# Remove the vOTUs that are not prevalent
Abundance_table_infants_log_filtered <- as.data.frame(Abundance_table_infants_log[, vOTU_keep])

# Order abundance table according to Sample metadata
Abundance_table_infants_log_filtered <- Abundance_table_infants_log_filtered[match(Sample_metadata_infants$NG_ID,
                                                                                   rownames(Abundance_table_infants_log_filtered)), ]

################
# INFANTS – Host-level
################

# Calculate prevalence
vOTU_nonZero <- colSums(Abundance_table_infants_host > 0) / nrow(Abundance_table_infants_host)

# Keep vOTUs present in >5% of samples
vOTU_keep <- colnames(Abundance_table_infants_host)[vOTU_nonZero > 0.05]

# Log transformation
pseudocount_vOTU <- min(Abundance_table_infants_host[Abundance_table_infants_host > 0]) / 2
Abundance_table_infants_host_log <- Abundance_table_infants_host
Abundance_table_infants_host_log[Abundance_table_infants_host_log == 0] <- pseudocount_vOTU
Abundance_table_infants_host_log <- log(Abundance_table_infants_host_log)

# Filter non-prevalent vOTUs
Abundance_table_infants_host_log_filtered <- as.data.frame(
  Abundance_table_infants_host_log[, vOTU_keep]
)

# Order according to metadata
Abundance_table_infants_host_log_filtered <- Abundance_table_infants_host_log_filtered[
  match(Sample_metadata_infants$NG_ID, rownames(Abundance_table_infants_host_log_filtered)), 
]

################
# MOTHERS - vOTU level
################

# Perform a log transformation and filter by prevalence
Abundance_table_mothers <- t(Abundance_table_mothers)

# Calculate the prevalence of each vOTU
vOTU_nonZero  <- colSums(Abundance_table_mothers>0) / nrow(Abundance_table_mothers) 

# Get the list of vOTUs present in more than 10% samples
vOTU_keep <- colnames(Abundance_table_mothers)[vOTU_nonZero > 0.1]

# Generate LOG transformed abundance tables
## CLR transformation not suitable here: data hardly compositional after read mapping
## Log transformation can help with data normality necessary for linear models
pseudocount_vOTU <- min(Abundance_table_mothers[Abundance_table_mothers > 0]) /2
Abundance_table_mothers_log <- Abundance_table_mothers
Abundance_table_mothers_log[Abundance_table_mothers_log  == 0] <- pseudocount_vOTU
Abundance_table_mothers_log <- log(Abundance_table_mothers_log)

# Remove the vOTUs that are not prevalent
Abundance_table_mothers_log_filtered <- as.data.frame(Abundance_table_mothers_log[, vOTU_keep])

# Order abundance table according to Sample metadata
Abundance_table_mothers_log_filtered <- Abundance_table_mothers_log_filtered[
  match(Sample_metadata_mothers$NG_ID,
        rownames(Abundance_table_mothers_log_filtered)), ]

################
# MOTHERS – Host-level
################

# Calculate prevalence
vOTU_nonZero <- colSums(Abundance_table_mothers_host > 0) / nrow(Abundance_table_mothers_host)

# Keep vOTUs present in >10% of samples
vOTU_keep <- colnames(Abundance_table_mothers_host)[vOTU_nonZero > 0.1]

# Log transformation
pseudocount_vOTU <- min(Abundance_table_mothers_host[Abundance_table_mothers_host > 0]) / 2
Abundance_table_mothers_host_log <- Abundance_table_mothers_host
Abundance_table_mothers_host_log[Abundance_table_mothers_host_log == 0] <- pseudocount_vOTU
Abundance_table_mothers_host_log <- log(Abundance_table_mothers_host_log)

# Filter non-prevalent vOTUs
Abundance_table_mothers_host_log_filtered <- as.data.frame(
  Abundance_table_mothers_host_log[, vOTU_keep]
)

# Order according to metadata
Abundance_table_mothers_host_log_filtered <- Abundance_table_mothers_host_log_filtered[
  match(Sample_metadata_mothers$NG_ID, rownames(Abundance_table_mothers_host_log_filtered)), 
]


##************************************************************************
# 2. Association analysis 
##************************************************************************

################
# INFANTS
################

phenotypes_to_test <- "Timepoint_categorical"

# We set W2 as the reference level
Sample_metadata_infants$Timepoint_categorical <- factor(
  Sample_metadata_infants$Timepoint_categorical,
  levels = c("W2", "M1", "M2", "M3", "M6", "M9", "M12")
)

# Run association with time at vOTU-level
vOTU_associations_mixed_infants <- associate_vOTUs_with_phenotypes(Sample_metadata_infants,
                                                                   "NG_ID", 
                                                                   Abundance_table_infants_log_filtered,
                                                                   phenotypes_to_test)

# Run association with time at host-level
vOTU_associations_mixed_infants_host <- associate_vOTUs_with_phenotypes(Sample_metadata_infants,
                                                                   "NG_ID", 
                                                                   Abundance_table_infants_host_log_filtered,
                                                                   phenotypes_to_test)

# Order associations results by FDR
vOTU_associations_mixed_infants <- vOTU_associations_mixed_infants %>% arrange(FDR)
vOTU_associations_mixed_infants_host <- vOTU_associations_mixed_infants_host %>% arrange(FDR)

# Save results
write.table(vOTU_associations_mixed_infants, "9_ASSOCIATION_PHENOTYPES/vOTU_Time_associations_mixed_infants_19092025.txt",
            sep="\t", row.names=F, quote = F)
write.table(vOTU_associations_mixed_infants_host, "9_ASSOCIATION_PHENOTYPES/vOTU_Time_associations_mixed_infants_host_19092025.txt",
            sep="\t", row.names=F, quote = F)


################
# MOTHERS
################

# We set P12 as the reference level
Sample_metadata_mothers$Timepoint_categorical <- factor(
  Sample_metadata_mothers$Timepoint_categorical,
  levels = c("P12", "P28", "B", "M3")
)

# Run association with time at vOTU-level
vOTU_associations_mixed_mothers <- associate_vOTUs_with_phenotypes(Sample_metadata_mothers,
                                                                   "NG_ID", 
                                                                   Abundance_table_mothers_log_filtered,
                                                                   phenotypes_to_test)

# Run association with time at host-level
vOTU_associations_mixed_mothers_host <- associate_vOTUs_with_phenotypes(Sample_metadata_mothers,
                                                                                  "NG_ID", 
                                                                                  Abundance_table_mothers_host_log_filtered,
                                                                        phenotypes_to_test)


# Subset pregnancy timepoints (P12, P28 and B) for host-level association
Sample_metadata_mothers_pregnancy <- Sample_metadata_mothers %>%
  dplyr::filter(Timepoint_categorical %in% c("P12", "P28", "B"))

Abundance_table_mothers_host_log_filtered_pregnancy <- Abundance_table_mothers_host_log_filtered[
  match(Sample_metadata_mothers_pregnancy$NG_ID,
        rownames(Abundance_table_mothers_host_log_filtered)), ]

# Run association with pregnancy timepoints at host-level
vOTU_associations_mixed_mothers_host_pregnancy <- associate_vOTUs_with_phenotypes(Sample_metadata_mothers_pregnancy,
                                                                        "NG_ID", 
                                                                        Abundance_table_mothers_host_log_filtered,
                                                                        phenotypes_to_test)


# Order associations results by FDR
vOTU_associations_mixed_mothers <- vOTU_associations_mixed_mothers %>% arrange(FDR)
vOTU_associations_mixed_mothers_host <- vOTU_associations_mixed_mothers_host %>% arrange(FDR)
vOTU_associations_mixed_mothers_host_pregnancy <- vOTU_associations_mixed_mothers_host_pregnancy %>% arrange(FDR)

# Save results
write.table(vOTU_associations_mixed_mothers, "9_ASSOCIATION_PHENOTYPES/vOTU_Time_associations_mixed_mothers_19092025.txt",
            sep="\t", row.names=F, quote = F)
write.table(vOTU_associations_mixed_mothers_host, "9_ASSOCIATION_PHENOTYPES/vOTU_Time_associations_mixed_mothers_host_19092025.txt",
            sep="\t", row.names=F, quote = F)
write.table(vOTU_associations_mixed_mothers_host_pregnancy, "9_ASSOCIATION_PHENOTYPES/vOTU_Time_associations_mixed_mothers_host_pregnancy_19092025.txt",
            sep="\t", row.names=F, quote = F)


##************************************************************************
# 3. Plots
##************************************************************************

##*****************
# 3.1 Infants
##*****************

###########################
# Settings
########################### 
association_heatmap_infant_vOTUs <- vOTU_associations_mixed_infants
association_heatmap_infant_host <- vOTU_associations_mixed_infants_host


association_heatmap_infant_host$Timepoint <- sub("^Timepoint_categorical", "", association_heatmap_infant_host$Term)
timepoint_levels <- c("M1", "M2", "M3", "M6", "M9", "M12")
association_heatmap_infant_host$Timepoint <- factor(association_heatmap_infant_host$Timepoint, levels = timepoint_levels)

timepoint_levels <- c("M1","M2","M3","M6","M9","M12")
gap <- 0.5  # gap between blocks

###########################
# Prepare vOTU-level data
###########################
vOTU_heatmap_long <- dplyr::left_join(
  association_heatmap_infant_vOTUs,
  Virus_metadata_infants[, c("Virus_ID", "Bacterial_genus_host")],
  by = c("vOTU" = "Virus_ID")
) %>%
  dplyr::mutate(
    Host = Bacterial_genus_host,
    Timepoint = sub("^Timepoint_categorical", "", Term),
    Timepoint = factor(Timepoint, levels = timepoint_levels)
  )

# Extra row for vOTU significance
vOTU_sig <- vOTU_heatmap_long %>%
  dplyr::group_by(vOTU, Host) %>%
  dplyr::summarise(sig = ifelse(any(FDR < 0.05), "*",""), .groups = "drop") %>%
  dplyr::mutate(Timepoint = "Sig", Estimate = NA)

vOTU_heatmap_long <- bind_rows(
  dplyr::select(vOTU_heatmap_long, vOTU, Timepoint, Estimate, Host),
  dplyr::select(vOTU_sig, vOTU, Timepoint, Estimate, Host, sig)
)

vOTU_heatmap_long$sig[is.na(vOTU_heatmap_long$sig)] <- ""
vOTU_heatmap_long$Timepoint <- factor(vOTU_heatmap_long$Timepoint,
                                      levels = c(timepoint_levels, "Sig"))

###########################
# Prepare host-level data
###########################
host_heatmap_long <- association_heatmap_infant_host %>%
  dplyr::mutate(
    Host = vOTU,
    Timepoint = sub("^Timepoint_categorical", "", Term),
    Timepoint = factor(Timepoint, levels = timepoint_levels)
  )

# Extra row for host significance
host_sig <- host_heatmap_long %>%
  dplyr::group_by(Host) %>%
  dplyr::summarise(sig = ifelse(any(FDR < 0.05), "*",""), .groups = "drop") %>%
  dplyr::mutate(Timepoint = "Sig", Estimate = NA)

host_heatmap_long <- bind_rows(
  dplyr::select(host_heatmap_long, Host, Timepoint, Estimate),
  dplyr::select(host_sig, Host, Timepoint, Estimate, sig)
)

host_heatmap_long$sig[is.na(host_heatmap_long$sig)] <- ""
host_heatmap_long$Timepoint <- factor(host_heatmap_long$Timepoint,
                                      levels = c(timepoint_levels, "Sig"))

###########################
# Order hosts and vOTUs
###########################

host_order <- host_heatmap_long %>%
  dplyr::filter(Timepoint=="M1") %>%
  dplyr::arrange(Estimate) %>%
  dplyr::pull(Host) %>% unique()

# Set manually for visualization
manual_host_order <- c("Staphylococcus","Lachnospira", "Streptococcus", "Dysosmobacter",
                       "Cutibacterium", "Phocaeicola", "Lactobacillus","Hungatella", 
                       "Phyllobacterium", "Haemophilus_D", "Veillonella_A", "Enterocloster",
                       "Lacticaseibacillus","Agathobacter","Faecalibacterium",
                       "Lactococcus", "Enterococcus", "Lawsonibacter", "Bifidobacterium",
                       "Sutterella","Pseudoflavonifractor_A", "Parabacteroides", "Bilophila", 
                       "Roseburia", "Blautia_A","Veillonella",
                       "Bacteroides", "Alistipes","Prevotella", "Ruthenibacterium", 
                       "Salmonella","Flavonifractor", "Escherichia", "Collinsella",
                       "Akkermansia", "Citrobacter",
                       "Ruminococcus_B", "Pauljensenia",
                       "Enterobacter","Clostridium", "Klebsiella")

# Apply to both datasets
vOTU_heatmap_long$Host  <- factor(vOTU_heatmap_long$Host,
                                  levels = manual_host_order)
host_heatmap_long$Host  <- factor(host_heatmap_long$Host,
                                  levels = manual_host_order)

# Split vOTUs by host and add gaps
vOTU_split <- vOTU_heatmap_long %>%
  dplyr::filter(Timepoint != "Sig") %>%
  dplyr::select(Host, vOTU) %>% dplyr::distinct() %>%
  dplyr::arrange(Host, vOTU)

x_pos <- c()
current_x <- 1
host_block <- vOTU_split %>% dplyr::group_by(Host) %>% dplyr::group_split()
for(block in host_block){
  n <- nrow(block)
  x_pos <- c(x_pos, current_x:(current_x + n - 1))
  current_x <- current_x + n + gap
}
vOTU_split$x <- x_pos
vOTU_heatmap_long <- dplyr::left_join(vOTU_heatmap_long, vOTU_split[, c("vOTU","x")], by="vOTU")

###########################
# Bottom heatmap: Host-level
###########################
xmax_host <- length(unique(host_heatmap_long$Host)) + 0.5
ymax_host <- length(levels(host_heatmap_long$Timepoint)) + 0.5

p_bottom <- ggplot(host_heatmap_long, aes(x=Host, y=Timepoint)) +
  geom_rect(xmin=0.5, xmax=xmax_host, ymin=0.5, ymax=ymax_host,
            fill=NA, colour="black", size=1.25) +
  geom_tile(aes(fill=Estimate), colour="lightgrey", linewidth=0.5) +
  geom_text(aes(label=sig), color="#333333", size=8) +
  scale_fill_gradient2(low="#003F45", mid="#F0F0F0", high="#F48C42", #low = "#5E4FA2", mid = "#F0F0F0", high = "#9E0142",
                       midpoint=0, name="Effect size", na.value="white") +
  scale_y_discrete(limits=c(timepoint_levels, "Sig")) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=11),
    axis.title.x=element_blank(),
    panel.grid = element_blank(),  
    axis.text.y=element_text(size=12),
    axis.title.y=element_text(size=14),
    panel.border=element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )  + 
  xlab("Bacterial host genus") 

###########################
# Top heatmap: vOTU-level
###########################

# Positions of blocks according to host
vOTU_positions <- vOTU_split %>%
  dplyr::group_by(Host) %>%
  dplyr::summarise(
    xmin = min(x) - 0.5,
    xmax = max(x) + 0.5,
    .groups = "drop"
  )

# Set height of heatmap
ymax_votu <- length(levels(vOTU_heatmap_long$Timepoint)) + 0.5

p_top <- ggplot(vOTU_heatmap_long, aes(x = x, y = Timepoint)) +
  geom_rect(data = vOTU_positions,
            aes(xmin = xmin, xmax = xmax,
                ymin = 0.5, ymax = ymax_votu),
            fill = NA, colour = "black", size = 1.25, inherit.aes = FALSE) +
  geom_tile(aes(fill = Estimate), colour = "lightgrey", linewidth = 0.5) +
  geom_text(aes(label = sig), color = "#333333", size = 8) +
  scale_fill_gradient2(low="#003F45", mid="#F0F0F0", high="#F48C42",
                       midpoint = 0, name = "Effect size", na.value = "white") +
  scale_x_continuous(breaks = vOTU_split$x,
                     labels = vOTU_split$vOTU,
                     expand = c(0, 0)) +
  scale_y_discrete(limits = c(timepoint_levels, "Sig")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=14),
    panel.grid = element_blank(),         
    panel.border = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  ylab("Timepoint")


###########################
# Combine and save plot
###########################
pdf("9_ASSOCIATION_PHENOTYPES/Plots/Linked_Host_vOTU_heatmap_simple_with_sig_topnames.pdf",
    width=15, height=7.5)
(p_top / p_bottom) + patchwork::plot_layout(heights=c(2,2))
dev.off()



##************************************************************************
# 3. Plots – Mothers
##************************************************************************

###########################
# Settings
########################### 
association_heatmap_mothers_vOTUs <- vOTU_associations_mixed_mothers
association_heatmap_mothers_host <- vOTU_associations_mixed_mothers_host

association_heatmap_mothers_host$Timepoint <- sub("^Timepoint_categorical", "", association_heatmap_mothers_host$Term)

# Maternal timepoints to show in plots (excluding reference P12)
timepoint_levels_mothers <- c("P28", "B", "M3")
association_heatmap_mothers_host$Timepoint <- factor(association_heatmap_mothers_host$Timepoint,
                                                     levels = timepoint_levels_mothers)

gap <- 0.5  # gap between host blocks

###########################
# Prepare vOTU-level data
###########################
vOTU_heatmap_long_mothers <- dplyr::left_join(
  association_heatmap_mothers_vOTUs,
  Virus_metadata_mothers[, c("Virus_ID", "Bacterial_genus_host")],
  by = c("vOTU" = "Virus_ID")
) %>%
  dplyr::mutate(
    Host = Bacterial_genus_host,
    Timepoint = sub("^Timepoint_categorical", "", Term),
    Timepoint = factor(Timepoint, levels = timepoint_levels_mothers)
  )

# Extra row for vOTU significance
vOTU_sig_mothers <- vOTU_heatmap_long_mothers %>%
  dplyr::group_by(vOTU, Host) %>%
  dplyr::summarise(sig = ifelse(any(FDR < 0.05), "*",""), .groups = "drop") %>%
  dplyr::mutate(Timepoint = "Sig", Estimate = NA)

vOTU_heatmap_long_mothers <- bind_rows(
  dplyr::select(vOTU_heatmap_long_mothers, vOTU, Timepoint, Estimate, Host),
  dplyr::select(vOTU_sig_mothers, vOTU, Timepoint, Estimate, Host, sig)
)

vOTU_heatmap_long_mothers$sig[is.na(vOTU_heatmap_long_mothers$sig)] <- ""
vOTU_heatmap_long_mothers$Timepoint <- factor(vOTU_heatmap_long_mothers$Timepoint,
                                              levels = c(timepoint_levels_mothers, "Sig"))

###########################
# Prepare host-level data
###########################
host_heatmap_long_mothers <- association_heatmap_mothers_host %>%
  dplyr::mutate(
    Host = vOTU,
    Timepoint = sub("^Timepoint_categorical", "", Term),
    Timepoint = factor(Timepoint, levels = timepoint_levels_mothers)
  )

# Extra row for host significance
host_sig_mothers <- host_heatmap_long_mothers %>%
  dplyr::group_by(Host) %>%
  dplyr::summarise(sig = ifelse(any(FDR < 0.05), "*",""), .groups = "drop") %>%
  dplyr::mutate(Timepoint = "Sig", Estimate = NA)

host_heatmap_long_mothers <- bind_rows(
  dplyr::select(host_heatmap_long_mothers, Host, Timepoint, Estimate),
  dplyr::select(host_sig_mothers, Host, Timepoint, Estimate, sig)
)

host_heatmap_long_mothers$sig[is.na(host_heatmap_long_mothers$sig)] <- ""
host_heatmap_long_mothers$Timepoint <- factor(host_heatmap_long_mothers$Timepoint,
                                              levels = c(timepoint_levels_mothers, "Sig"))

###########################
# Order hosts and vOTUs
###########################

# Set order manually (for better visualization)
host_order <- c("Lactococcus","Agathobacter","Parabacteroides","Faecalibacterium",
                "Acetatifactor","CAG.269","Barnesiella","Lachnospira",
                "UBA11524","Ruminiclostridium_E","Bacteroides","Alistipes_A",
                "Dysosmobacter","CAG.177","X51.20","Odoribacter","Alistipes",
                "Eubacterium_F","ER4","CAG.273","Faecousia","UBA1417",
                "Phascolarctobacterium_A","CAG.245","UBA7173","Sutterella",
                "Scatacola_A","Desulfovibrio","Akkermansia","CAG.302",
                "CAG.353","Eubacterium_R","CAG.115","Prevotella",
                "Ruminococcus_E","CAG.217","Eisenbergiella","Ruminococcus_D",
                "CAG.267","Victivallis","Limiplasma","Coprococcus",
                "Lawsonibacter","Phascolarctobacterium","Fusicatenibacter",
                "CAG.170","SFEL01","Vescimonas","Gemmiger","Dialister",
                "Collinsella","Dorea_A","Copromonas","CAG.103",
                "Ruminococcus_C","Phocaeicola","Anaerobutyricum","Bilophila",
                "Blautia_A","Roseburia","Bifidobacterium","Mediterraneibacter")

hosts <- c("Lactococcus","Agathobacter","Parabacteroides","Faecalibacterium",
           "Lachnospira","Bacteroides","UBA7173","Scatacola_A","Dysosmobacter",
           "Ruminococcus_E","Ruminococcus_D","Alistipes","Vescimonas","Blautia_A",
           "Roseburia","Bifidobacterium")

positions <- c(2,4,6,9,12,22,29,32,34,36,38,47,50,52,54,56)

remaining <- setdiff(host_order, hosts)
new_order <- rep(NA, length(host_order))
new_order[positions] <- hosts
new_order[is.na(new_order)] <- remaining


vOTU_heatmap_long_mothers$Host <- factor(vOTU_heatmap_long_mothers$Host,
                                         levels = new_order)

host_heatmap_long_mothers$Host <- factor(host_heatmap_long_mothers$Host,
                                         levels = new_order)

# Split vOTUs by host and add gaps
vOTU_split_mothers <- vOTU_heatmap_long_mothers %>%
  dplyr::filter(Timepoint != "Sig") %>%
  dplyr::select(Host, vOTU) %>% dplyr::distinct() %>%
  dplyr::arrange(Host, vOTU)

x_pos <- c()
current_x <- 1
host_block <- vOTU_split_mothers %>% dplyr::group_by(Host) %>% dplyr::group_split()
for(block in host_block){
  n <- nrow(block)
  x_pos <- c(x_pos, current_x:(current_x + n - 1))
  current_x <- current_x + n + gap
}
vOTU_split_mothers$x <- x_pos
vOTU_heatmap_long_mothers <- dplyr::left_join(vOTU_heatmap_long_mothers,
                                              vOTU_split_mothers[, c("vOTU","x")], by="vOTU")

###########################
# Bottom heatmap: host-level
###########################
xmax_host <- length(unique(host_heatmap_long_mothers$Host)) + 0.5
ymax_host <- length(levels(host_heatmap_long_mothers$Timepoint)) + 0.5

p_bottom <- ggplot(host_heatmap_long_mothers, aes(x=Host, y=Timepoint)) +
  geom_rect(xmin=0.5, xmax=xmax_host, ymin=0.5, ymax=ymax_host,
            fill=NA, colour="black", size=1.25) +
  geom_tile(aes(fill=Estimate), colour="lightgrey", linewidth=0.5) +
  geom_text(aes(label=sig), color="#333333", size=8) +
  scale_fill_gradient2(low="#003F45", mid="#F0F0F0", high="#F48C42",
                       midpoint=0, name="Effect size", na.value="white") +
  scale_y_discrete(limits=c(timepoint_levels_mothers, "Sig")) +
  theme_minimal() +
  theme(
    axis.text.x=element_text(angle=90, hjust=1, vjust=0.5, size=11),
    axis.title.x=element_blank(),
    panel.grid = element_blank(),  
    axis.text.y=element_text(size=12),
    axis.title.y=element_text(size=14),
    panel.border=element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  )  + 
  xlab("Bacterial host genus") 

###########################
# Top heatmap: vOTU-level
###########################
vOTU_positions <- vOTU_split_mothers %>%
  dplyr::group_by(Host) %>%
  dplyr::summarise(
    xmin = min(x) - 0.5,
    xmax = max(x) + 0.5,
    .groups = "drop"
  )

ymax_votu <- length(levels(vOTU_heatmap_long_mothers$Timepoint)) + 0.5

p_top <- ggplot(vOTU_heatmap_long_mothers, aes(x = x, y = Timepoint)) +
  geom_rect(data = vOTU_positions,
            aes(xmin = xmin, xmax = xmax,
                ymin = 0.5, ymax = ymax_votu),
            fill = NA, colour = "black", size = 1.25, inherit.aes = FALSE) +
  geom_tile(aes(fill = Estimate), colour = "lightgrey", linewidth = 0.5) +
  geom_text(aes(label = sig), color = "#333333", size = 8) +
  scale_fill_gradient2(low="#003F45", mid="#F0F0F0", high="#F48C42",
                       midpoint = 0, name = "Effect size", na.value = "white") +
  scale_x_continuous(breaks = vOTU_split_mothers$x,
                     labels = vOTU_split_mothers$vOTU,
                     expand = c(0, 0)) +
  scale_y_discrete(limits = c(timepoint_levels_mothers, "Sig")) +
  theme_minimal() +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 12),
    axis.title.x=element_blank(),
    axis.title.y=element_text(size=14),
    panel.grid = element_blank(),         
    panel.border = element_blank(),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +
  ylab("Timepoint")

###########################
# Combine and save
###########################
pdf("9_ASSOCIATION_PHENOTYPES/Plots/Linked_Host_vOTU_heatmap_mothers_with_sig_topnames.pdf",
    width=15, height=5)
(p_top / p_bottom) + patchwork::plot_layout(heights=c(2,2))
dev.off()

