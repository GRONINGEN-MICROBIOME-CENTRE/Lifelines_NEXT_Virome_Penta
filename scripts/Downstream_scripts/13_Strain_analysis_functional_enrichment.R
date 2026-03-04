################################################################################
##### LL-NEXT: Viral strain transmission - Functional enrichment
### Author(s): Asier Fernández-Pato
### Last updated: 9th December, 2025
################################################################################

#****************
# Load libraries
#****************
library(dplyr)
library(tidyr)
library(pbapply)
library(data.table)
library(ggplot2)
library(treemapify)

# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_Analysis/VIRUSES/")

#****************
# Define functions
#****************
# Function to filter foldseek hits
filter_foldseek_hits <- function(df) {
  df %>%
    # Apply thresholds
    filter(qtmscore > 0.7, prob >= 0.9, qcov >= 0.7) %>%
    # Arrange by query, then descending prob, then descending qtmscore
    arrange(query, desc(prob), desc(qtmscore)) %>%
    # Keep only the top hit per query
    group_by(query) %>%
    dplyr::slice(1) %>%
    ungroup()
}

# Function to clean GO columns
clean_GO <- function(col) {
  remove_prefixes <- c( # remove entries starting with any of these prefixes
    "IEA:", "IDA:", "ISS:", "IPI:", "IBA:", "ISO:", "TAS:", "HDA:",
    "IMP:", "NAS:", "IC:", "ISM:", "EXP:"
  )
  pattern <- paste0("^(", paste0(remove_prefixes, collapse="|"), ")")
  str_split(col, "[;,]") %>% 
    map_chr(~ {
      x <- str_trim(.x)
      x <- x[!str_detect(x, pattern)]  # remove unwanted entries
      x <- x[x != ""]
      if (length(x) == 0) x <- "-"
      paste(x, collapse="; ")
    })
}

# Function to clean Keywords ad Pthway descriptions
clean_keywords_description <- function(col) {
  str_split(col, ";") %>%
    map_chr(~ {
      x <- str_trim(.x)
      x <- x[!x %in% c("-", "Reference proteome", "Proteomics identification", "3D-structure")]
      x <- x[x != ""]
      if (length(x) == 0) x <- "-"
      paste(x, collapse="; ")
    })
}


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
Abundance_table_BM <- read.delim("BREASTMILK/BM_LLNEXT_Viral_Abundance_Table_29072025.txt")

# Read viral metadata tables
Viral_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")
Viral_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_INFANTS_07052025.txt")
Viral_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_MOTHERS_07052025.txt")

# Read processed inStrain results 
inStrain_results <- read.delim("10_STRAIN_TRANSMISSION/inStrain_results_processed.txt")
inStrain_results_BM <- read.delim("BREASTMILK/inStrain_results_BM_processed.txt")

##************************************************************************
# 2. Identification of shared and vertically transmitted vOTUs 
#*************************************************************************

##################
# FC-FC comparison
##################

# Identify shared (and not shared) vOTUs (with any maternal timepoint and infant W2-M3)
inStrain_results_pairs_shared_W2_M3 <- inStrain_results %>%
  filter(Mother_infant_sharing == "Yes" & Mother_Infant_pair == "Pair" &
           Infant_timepoint %in% c("W2", "M1", "M2", "M3"))

# Non-shared vOTUs are defined as those detected in maternal samples (present) but not in any infant sample (at vOTU level)
shared_vOTUs_FC <- unique(inStrain_results_pairs_shared_W2_M3$Virus_ID)
maternal_vOTUs_FC <- rownames(Abundance_table_mothers)
infant_vOTUs_FC <- rownames(Abundance_table_infants)
non_shared_vOTUs_FC <- maternal_vOTUs_FC[!maternal_vOTUs_FC %in% infant_vOTUs_FC]

##################
# Milk-FC comparison
##################

# Identify shared (and not shared) vOTUs (with any maternal timepoint and infant W2-M3)
inStrain_results_BM_pairs_shared_W2_M3 <- inStrain_results_BM %>%
  filter(Mother_infant_sharing == "Yes" & Mother_Infant_pair == "Pair" &
           infant_timepoint %in% c("W2", "M1", "M2", "M3"))

# Non-shared vOTUs are defined as those detected in maternal samples (present) but not in any infant sample (at vOTU level)
shared_vOTUs_BM <- unique(inStrain_results_BM_pairs_shared_W2_M3$Virus_ID)
maternal_vOTUs_BM <- rownames(Abundance_table_BM)
non_shared_vOTUs_BM <- maternal_vOTUs_BM[!maternal_vOTUs_BM %in% shared_vOTUs_BM]


##################
# Combined results
##################
shared_vOTUs <- unique(c(shared_vOTUs_FC, shared_vOTUs_BM))
non_shared_vOTUs <- unique(c(non_shared_vOTUs_FC, non_shared_vOTUs_BM))

# Keep only in shared_vOTUs those present in both groups
non_shared_vOTUs <- non_shared_vOTUs[!non_shared_vOTUs %in% shared_vOTUs]


##************************************************************************
# 3. Functional enrichment (protein families) in shared and vertically transmitted vOTUs 
#*************************************************************************

# Predicted proteins from all vOTU representative genomes (n=31,209) were clustered into protein families (2-step clustering approach using MMseqs2)

# Load clustering results
protein_clusters_long_format <- read.table("10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/vOTUs_proteins_clusters.tsv",
                                           sep = "\t", header = FALSE, stringsAsFactors = FALSE)
protein_cluster_names <- read.table("10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/vOTUs_proteins_clusters_names.tsv",
                                    sep = "\t", header = FALSE, stringsAsFactors = FALSE)
colnames(protein_clusters_long_format) <- c("Cluster_Representative", "Protein_ID")
colnames(protein_cluster_names) <- c("Cluster_ID", "Cluster_Representative")
protein_clusters_long_format <- left_join(protein_clusters_long_format, protein_cluster_names, by = "Cluster_Representative")
  
# Remove all clusters exclusively containing proteins from potential contaminant vOTUs
contaminants <- c("LLNEXT_25294", "LLNEXT_79800", "LLNEXT_BM_1", "LLNEXT_BM_103")
protein_clusters_long_format <- protein_clusters_long_format %>%
  filter(!str_detect(Protein_ID, paste0(contaminants, "_", collapse = "|")))

# Generate a mapping file with the Virus_ID corresponding to each protein_ID
Viral_protein_mapping <- data.frame(Protein_ID = protein_clusters_long_format$Protein_ID)
Viral_protein_mapping$Virus_ID <- sub("_([^_]+)$", "", protein_clusters_long_format$Protein_ID)

# Group proteins by their cluster representative
protein_clusters <- protein_clusters_long_format %>%
  group_by(Cluster_ID) %>%
  summarise(Proteins = list(Protein_ID)) # Store proteins as a list 

# Identify shared and not shared proteins
shared_proteins <- Viral_protein_mapping$Protein_ID[Viral_protein_mapping$Virus_ID %in% shared_vOTUs]
non_shared_proteins <- Viral_protein_mapping$Protein_ID[Viral_protein_mapping$Virus_ID %in% non_shared_vOTUs]

# Estimate the total number of proteins and the total number of shared proteins
n_proteins <- length(protein_clusters_long_format$Protein_ID)
n_shared_proteins <- length(shared_proteins)
n_non_shared_proteins <- length(non_shared_proteins)

# Count the number of proteins from each cluster present in shared/non-shared vOTUs
protein_clusters <- protein_clusters %>%
  mutate(
    shared = pbsapply(Proteins, function(protein_list) sum(protein_list %in% shared_proteins)),
    non_shared = pbsapply(Proteins, function(protein_list) sum(protein_list %in% non_shared_proteins))
  )

####################################################
# A. Filter Clusters: At least 10 proteins
####################################################

# Filter based on vOTU diversity and cluster size
protein_clusters_shared <- protein_clusters %>%
  filter((shared + non_shared) >= 10)

####################################################
# B. Enrichment: Fisher test
####################################################

# Enrichment for shared proteins
protein_clusters_shared$fisher_p_value <- NA
protein_clusters_shared$fisher_odds_ratio <- NA

for (i in 1:nrow(protein_clusters_shared)) {
  a <- protein_clusters_shared$shared[i]
  b <- protein_clusters_shared$non_shared[i]
  c <- n_shared_proteins - a
  d <- n_non_shared_proteins - b
  
  contingency <- matrix(c(a, b, c, d), nrow=2, byrow=TRUE)
  fisher_res <- fisher.test(contingency, alternative="two.sided")
  
  protein_clusters_shared$fisher_p_value[i] <- fisher_res$p.value
  protein_clusters_shared$fisher_odds_ratio[i] <- fisher_res$estimate
}

protein_clusters_shared$FDR <- p.adjust(protein_clusters_shared$fisher_p_value, method="BH")

# Get protein clusters enriched among shared vOTUs
enriched_clusters_shared <- protein_clusters_shared %>%
  filter(FDR < 0.05 & fisher_odds_ratio > 1)

protein_clusters_shared$Enriched <- ifelse(
  protein_clusters_shared$FDR < 0.05 & protein_clusters_shared$fisher_odds_ratio > 1,
  TRUE,FALSE)

# Reorder and add represenative protein as a column (to save the table)
enriched_clusters_shared <- enriched_clusters_shared %>%
  arrange(FDR)
enriched_clusters_shared_save <- enriched_clusters_shared
enriched_clusters_shared_save$Proteins <- NULL
enriched_clusters_shared_save <- enriched_clusters_shared_save %>%
  left_join(
    protein_clusters_long_format %>%
      select(Cluster_ID, Cluster_Representative) %>%
      distinct(Cluster_ID, .keep_all = TRUE),
    by = "Cluster_ID"
  )

# Generate a plot with enrichment results
pdf('10_STRAIN_TRANSMISSION/Plots/PF_enrichment_density.pdf', width = 4.9, height = 3)
ggplot(protein_clusters_shared, aes(x = fisher_odds_ratio, fill = Enriched)) +
  geom_density(alpha = 0.6) +
  geom_rug(aes(color = Enriched), sides = "b", alpha = 0.5, length = unit(0.05, "npc")) +  
  scale_fill_manual(values = c("grey70", "#FB8072"), labels = c("Non-enriched", "Enriched")) +
  scale_color_manual(values = c("grey50", "#FB8072")) +
  scale_x_log10() +
  scale_y_continuous(expand = expansion(mult = c(0.1, 0.05))) +
  theme_minimal(base_size = 14) +
  labs(x = "Odds ratio (log scale)", y = "Density", fill = "Enrichment", color = "Enrichment") +
  theme(
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 15),
    axis.ticks = element_line(color = "black"),  
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 12)
  ) 
dev.off()

##************************************************************************
# 4. Functional annotation enriched PFs (HMM-based vs PHROGs and KOs)
#*************************************************************************

# Extract ID of enriched PFs and their representatives
Enriched_PFs <- enriched_clusters_shared$Cluster_ID
write.table(Enriched_PFs, file = "10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/Enriched_PFs.txt", quote = FALSE, row.names = FALSE,col.names = FALSE)
rep_enriched_PFs <- unique(protein_clusters_long_format[protein_clusters_long_format$Cluster_ID %in% Enriched_PFs, "Cluster_Representative"])
write.table(Enriched_PFs, file = "10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/Enriched_PFs_representative.txt", quote = FALSE, row.names = FALSE,col.names = FALSE)

# Load annotation results of enriched clusters
PHROGs_annotation <- read.delim("10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/HMMER_processed_table.tsv")
KOfam_annotation <- read.delim("10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/ko-annotations.tsv", check.names = F)
ADF_annotation <- read.delim("10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/Enriched_PFs_rep_defense_finder_systems.tsv", check.names = F)

# Read PHROGs metadata
metadata_phrog <- read.delim("5_DGR_ANALYSIS/phrog_annot_v4.tsv", check.names = F)
metadata_phrog$phrog <- paste0("phrog_", metadata_phrog$phrog)

# Combine HMMER results and metadata to get the annotations for each protein (order by E-value)
Protein_annotations_PHROGs <- left_join(PHROGs_annotation, metadata_phrog,  by = c("Match" = "phrog"))
Protein_annotations_PHROGs <- Protein_annotations_PHROGs[order(Protein_annotations_PHROGs$E_value), ]
Protein_annotations_PHROGs <- Protein_annotations_PHROGs[, c("Protein_ID","Match", "E_value", "Description", "annot", "category")]
colnames(Protein_annotations_PHROGs) <- c("Protein_ID","Match", "E_value", "Protein_Info", "DESCRIPTION", "CATEGORY")

# Keep only top annotation (min E-value)
Protein_annotations_PHROGs <- Protein_annotations_PHROGs %>%
  group_by(Protein_ID) %>%
  slice_min(order_by = E_value, with_ties = FALSE) %>%
  ungroup()

# Filter only significant matches from KOFAM results
KOfam_annotation <- KOfam_annotation[-1, -1]
KOfam_annotation <- KOfam_annotation[!is.na(KOfam_annotation$thrshld) & KOfam_annotation$thrshld != "", ] # remove those with NA/empty in threshold
KOfam_annotation$score <- as.numeric(KOfam_annotation$score)
KOfam_annotation$thrshld <- as.numeric(KOfam_annotation$thrshld)
KOfam_annotation <- KOfam_annotation[KOfam_annotation$score >= KOfam_annotation$thrshld,]
KOfam_annotation <- KOfam_annotation %>% # If multiple annotations for a protein, select the one with lowest E-value
  group_by(`gene name`) %>%
  slice_min(order_by = `E-value`, n = 1, with_ties = FALSE) %>%
  ungroup()

# Check number of PFs with significant matches to known proteins (n=625)
PF_with_PHROGs_match <- Protein_annotations_PHROGs$Protein_ID
PF_with_KOfam_match <- KOfam_annotation$`gene name`
PF_with_ADF_match <- ADF_annotation$protein_in_syst
PF_with_match <- unique(c(PF_with_PHROGs_match, PF_with_KOfam_match, PF_with_ADF_match))
PF_KOfam_not_PHROGs_match <- setdiff(PF_with_KOfam_match, PF_with_PHROGs_match)

# Check number of PFs with informative annotations (n=302)
PF_with_PHROGs_ann <- Protein_annotations_PHROGs$Protein_ID[Protein_annotations_PHROGs$CATEGORY != "unknown function"]
PF_with_KOfam_ann <- KOfam_annotation$`gene name`[KOfam_annotation$`KO definition` != "uncharacterized protein"] #2 cases
PF_with_ADF_ann <- ADF_annotation$protein_in_syst
PF_with_ann <- unique(c(PF_with_PHROGs_ann, PF_with_KOfam_ann, PF_with_ADF_ann))
PF_KOfam_not_PHROGs <- setdiff(PF_with_KOfam_ann, PF_with_PHROGs_ann) #8
PF_ADF_not_PHROGs <- setdiff(PF_with_ADF_ann, PF_with_PHROGs_ann) #2

# Get known PF annotations
# Classify into categories: PHROGS (each category), KOFAM and ADF
Protein_annotations_PHROGs_known <- Protein_annotations_PHROGs[Protein_annotations_PHROGs$CATEGORY != "unknown function",]
Protein_annotations_PHROGs_known <- Protein_annotations_PHROGs_known[!Protein_annotations_PHROGs_known$Protein_ID %in% 
                                                                       ADF_annotation$protein_in_syst,]
Protein_annotations_KOfam_not_PHROGs <- KOfam_annotation[KOfam_annotation$`gene name` %in%PF_KOfam_not_PHROGs,]
Protein_annotations_ADF <- ADF_annotation

category_counts <- Protein_annotations_PHROGs_known %>%
  count(CATEGORY) %>%
  bind_rows(tibble(CATEGORY = "KOFAM", n = nrow(Protein_annotations_KOfam_not_PHROGs))) %>%
  bind_rows(tibble(CATEGORY = "AntiDefenseFinder", n = nrow(Protein_annotations_ADF))) %>%
  bind_rows(
    tibble(
      CATEGORY = "Unknown",
      n = length(Enriched_PFs) - nrow(Protein_annotations_PHROGs_known) -
        nrow(Protein_annotations_KOfam_not_PHROGs) - nrow(Protein_annotations_ADF)
    )
  ) %>%
  mutate(Percent = 100 * n / sum(n))

# Generate a tree map with annotated and unnanotated proteins
category_colors <- c("#264653", "#2A9D8F", "#8AB17D", "#E9C46A", "#F4A261", "#E76F51",
                     "#A68A8A", "#BFB1A8", "#8C7AA9", "#C6D57E", "#D4A5A5", "#B5C2B7")

category_order <- c("DNA, RNA and nucleotide metabolism", "connector", "head and packaging", 
                    "integration and excision", "lysis", "moron, auxiliary metabolic gene and host takeover", 
                    "other", "tail", "transcription regulation", "KOFAM", "AntiDefenseFinder", "Unknown")
category_counts$CATEGORY <- factor(category_counts$CATEGORY, levels = category_order)


pdf("10_STRAIN_TRANSMISSION/Plots/Functional_categories_map.pdf", width = 14, height = 2)
ggplot(category_counts, aes(area = n, fill = CATEGORY)) +
  geom_treemap(color = "black", size = 0.3) +
  scale_fill_manual(values = category_colors) +
  labs(fill = "Functional category",
    title = "Functional annotation of enriched protein families") +
  theme_minimal(base_size = 14) +
  theme(
    plot.title = element_blank(),
    legend.title = element_text(size = 15),
    legend.text = element_text(size = 13),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    panel.grid = element_blank(),
    plot.background = element_blank()
  )
dev.off()

################
# Save results
################

# Save file with unannotated proteins to be used in structural prediction
# Save them in order (from less to more members)
cluster_summary <- protein_clusters_long_format %>%
  group_by(Cluster_ID) %>%
  summarise(N_members = n(), .groups = "drop")

rep_enriched_PFs <- enriched_clusters_shared_save$Cluster_Representative
PFs_reps_wo_PHROGs_annotations <- rep_enriched_PFs[!rep_enriched_PFs %in% PF_with_PHROGs_ann]
PFs_wo_PHROGs_annotations <- unique(protein_clusters_long_format[protein_clusters_long_format$Cluster_Representative
                                                          %in% PFs_reps_wo_PHROGs_annotations, "Cluster_ID"])

PFs_wo_PHROGs_annotations_ordered <- cluster_summary %>%
  filter(Cluster_ID %in% PFs_wo_PHROGs_annotations) %>%
  arrange(N_members) %>%
  pull(Cluster_ID)

write.table(PFs_wo_PHROGs_annotations, file = "10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/Enriched_PFs_wo_PHROGs_annotations.txt", quote = FALSE, row.names = FALSE,col.names = FALSE)


# Add sequence-similarity based annotations and save enrichment table
# Add first PHROGs annotations, and then KOFAM and ADF
enriched_clusters_shared_save <- enriched_clusters_shared_save %>%
  left_join(
    Protein_annotations_PHROGs %>%
      select(Protein_ID, DESCRIPTION, CATEGORY),
    by = c("Cluster_Representative" = "Protein_ID")
  ) %>%
  dplyr::rename(
    Annotation = DESCRIPTION,
    PHROG_category = CATEGORY
  )

enriched_clusters_shared_save <- enriched_clusters_shared_save %>%
  left_join(
    Protein_annotations_KOfam_not_PHROGs %>%
      select(`gene name`, `KO definition`),
    by = c("Cluster_Representative" = "gene name")
  )

enriched_clusters_shared_save <- enriched_clusters_shared_save %>%
  left_join(
    Protein_annotations_ADF %>%
      select(protein_in_syst, name_of_profiles_in_sys),
    by = c("Cluster_Representative" = "protein_in_syst")
  )

# Merge all annotations
enriched_clusters_shared_save <- enriched_clusters_shared_save %>%
  mutate(
    Annotation = dplyr::coalesce(
      Annotation,               
      `KO definition`,          
      name_of_profiles_in_sys 
    )
  )

write.table(enriched_clusters_shared_save,"10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/RESULTS/PF_enrichment.txt",
            sep = "\t", row.names = F, quote = FALSE) 

##************************************************************************
# 5. Functional annotation by structural similarity
#*************************************************************************

####################################
# Process Foldseek results
####################################

# By running ColabFold-Foldseek on unnanotated clusters, we get the following results
foldseek_cols <- c("query", "target", "pident", "alnlen", "evalue", "bits",
  "qlen", "tlen", "prob", "lddt", "qtmscore", "ttmscore", "alntmscore",
  "qcov", "tcov", "qstart", "qend", "tstart", "tend","qseq", "tseq")

# Read the file and assign column names
BFVD_annotation <- read.delim(
  "10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/BFVD_results.tsv",
  header = FALSE, sep = "\t", col.names = foldseek_cols
)
AF_annotation <- read.delim(
  "10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/AFDB_results.tsv",
  header = FALSE, sep = "\t", col.names = foldseek_cols
)
PDB_annotation <- read.delim(
  "10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/PDB_results.tsv",
  header = FALSE, sep = "\t", col.names = foldseek_cols
)

# Filter results
BFVD_sig <- filter_foldseek_hits(BFVD_annotation)
AFDB_sig   <- filter_foldseek_hits(AF_annotation)
PDB_sig  <- filter_foldseek_hits(PDB_annotation)

# Modify tables with the significant results
BFVD_sig <- BFVD_sig %>%
  mutate(Cluster_ID = str_extract(query, "^[^_]+"))
AFDB_sig <- AFDB_sig %>%
  mutate(uniprot = sub("AF-([^-]+)-.*", "\\1", target),
         Cluster_ID = str_extract(query, "^[^_]+"))
PDB_sig <- PDB_sig %>%
  mutate(
    pdb_id = str_extract(target, "^[^-]+"),          # extract PDB ID
    assembly = str_extract(target, "assembly\\d+"), # extract assembly
    chain = str_extract(target, "(?<=_)[^_]+$"),
    Cluster_ID = str_extract(query, "^[^_]+")
  )

# Keep only clusters with no annotation by sequence similarity (exclude if present 10 clusters with only KOFAM, ADF annotations)
representatives_to_exclude <- c(PF_KOfam_not_PHROGs, PF_ADF_not_PHROGs)
clusters_to_exclude <- unique(protein_clusters_long_format[protein_clusters_long_format$Cluster_Representative %in%
                                                      representatives_to_exclude, "Cluster_ID"])

BFVD_sig <- BFVD_sig[!BFVD_sig$Cluster_ID %in% clusters_to_exclude,]
AFDB_sig <- AFDB_sig[!AFDB_sig$Cluster_ID %in% clusters_to_exclude,]
PDB_sig <- PDB_sig[!PDB_sig$Cluster_ID %in% clusters_to_exclude,]

# Save results tables with significant hits
write.table(BFVD_sig, "10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/RESULTS/BFVD_sig_hits.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(AFDB_sig, "10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/RESULTS/AFDB_sig_hits.tsv", sep="\t", quote=FALSE, row.names=FALSE)
write.table(PDB_sig, "10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/RESULTS/PDB_sig_hits.tsv", sep="\t", quote=FALSE, row.names=FALSE)

# Save list of UNIPROT/PDB IDs
write.table(unique(BFVD_sig$target), "10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/BFVD_uniprot_ids.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(unique(AFDB_sig$uniprot), "10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/AFDB_uniprot_ids.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)
write.table(unique(PDB_sig$pdb_id), "10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/PDB_ids.txt", quote=FALSE, row.names=FALSE, col.names=FALSE)

# Check the number of PFs with at least one significant match to known protein structures
PFs_with_sig_structural_match <- unique(c(BFVD_sig$query, AFDB_sig$query, PDB_sig$query)) # 288 (42.1% of unannotated PFs) 
PFs_with_sig_structural_match <- sub("_.*", "", PFs_with_sig_structural_match)

# Check the overlap with proteins with any matches by sequence-homology based methods
PFs_with_sequence_hom_match <- protein_cluster_names[protein_cluster_names$Cluster_Representative %in% PF_with_match,
                     "Cluster_ID"]
PFs_str_not_seq_matches <- setdiff(PFs_with_sig_structural_match, PFs_with_sequence_hom_match) #149 with no sequence-homology based matches


####################################
# Processing Uniprot annotations
####################################

# Load results with Uniprot annotations of BFVD and AFDB matches
BFVD_annotations <- read.delim("10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/BFVD_uniprot_annotation.tsv",
  header = T, sep = "\t")
AFDB_annotations <- read.delim("10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/AFDB_uniprot_annotation.tsv",
                               header = T, sep = "\t")
PDB_annotations <- read.delim("10_STRAIN_TRANSMISSION/FUNCTIONAL_ENRICHMENT/PDB_uniprot_annotation_filt.tsv",
                               header = T, sep = "\t")

#***************************
# A. Processing BFVD results
#***************************

# Remove rows with "_" (no ID) in UniProt_ID column
BFVD_annotations <- BFVD_annotations[!grepl("-", BFVD_annotations$UniProt_ID), ]

# Process annotations
BFVD_annotations <- BFVD_annotations %>%
  mutate(
    GO_CC = clean_GO(GO_CC),
    GO_MF = clean_GO(GO_MF),
    GO_BP = clean_GO(GO_BP),
    Keywords = clean_keywords_description(Keywords)
  ) %>%
  select(-Pathway_IDs, -Pathway_DBs, -Pathway_descriptions)

#***************************
# B. Processing AFDB results
#***************************

# Remove rows with "_" (no ID) in UniProt_ID column
AFDB_annotations <- AFDB_annotations[!grepl("-", AFDB_annotations$UniProt_ID), ]

# Process annotations
AFDB_annotations <- AFDB_annotations %>%
  mutate(
    GO_CC = clean_GO(GO_CC),
    GO_MF = clean_GO(GO_MF),
    GO_BP = clean_GO(GO_BP),
    Keywords = clean_keywords_description(Keywords),
    Pathway_descriptions = clean_keywords_description(Pathway_descriptions)
  ) %>%
  select(-Pathway_DBs)


#***************************
# C. Processing PDB results
#***************************

# Remove rows with "_" (no ID) in UniProt_ID column
PDB_annotations <- PDB_annotations[!grepl("-", PDB_annotations$UniProt_ID), ]

# Process annotations
PDB_annotations <- PDB_annotations %>%
  mutate(
    GO_CC = clean_GO(GO_CC),
    GO_MF = clean_GO(GO_MF),
    GO_BP = clean_GO(GO_BP),
    Keywords = clean_keywords_description(Keywords),
  )

# Add the annotations to the tables with all significant matches per DB
BFVD_ann_results <- BFVD_sig %>%
  left_join(BFVD_annotations, by = c("target" = "UniProt_ID"))
AFDB_ann_results <- AFDB_sig %>%
  left_join(AFDB_annotations, by = c("uniprot" = "UniProt_ID"))
PDB_ann_results <- PDB_sig %>%
  left_join(PDB_annotations, by = c("pdb_id" = "PDB_ID"))

####################################
# Exploring interesting proteins/potential AMGs
####################################

# Ubiquitin-like modifier-activating enzyme 5: (Ubiquitin-like proteins/domains) (some in PDB trustable)
# Endoribonuclease HigB1 (trusted)
# Thoeris anti-defense 2-like domain-containing protein: anti-Thoeris protein Tad2 (trusted)
# Dolichol phosphate-mannose biosynthesis regulatory protein (?)
# Ferritin-like proteins / (manganase catalase) (seems to have similar fold to capsid from HK97-type phages)
# Phosphocarrier protein HPr (not trusted)
# Subtilisin BPN' (No. Some subtilisin-like domains)
# Histidine decarboxylase proenzyme (not trusted)
# Sugar protein SusD (probably bacterial)
# Homologues of ESCRT-I trafficking complex (or other ESCRTs) (trusted)
# PstA (homeostasis in bacteria) (found once in phage): the first PII-like protein with enzymatic activity has been identified in a bacteriophage and shown to catalyze lysis of S-adenosylmethionine (SAM). phage-encoded SAM degrading enzyme (anti-RM/anti-crispr Svi3-3)

