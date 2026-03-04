################################################################################
##### LL-NEXT: Metadata and abundance table processing
### Author(s): Asier Fernández-Pato
### Last updated: 19th December, 2025
################################################################################

# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_Analysis/VIRUSES/")

#****************
# Load libraries
#****************
library(dplyr)
library(Biostrings)
library(UpSetR)
library(ggplot2)
library(rentrez)
library(xml2)
library(purrr)
library(tidyr)

#****************
# Define functions
#****************
# Function to count the number of viruses with each completeness group in a row
count_viruses <- function(row) {
  viruses <- row[!is.na(row)]
  complete_quality_count <- sum(viruses %in% CheckV_quality_all$virus_ID[CheckV_quality_all$quality == "Complete"])
  high_quality_count <- sum(viruses %in% CheckV_quality_all$virus_ID[CheckV_quality_all$quality == "High-quality"])
  medium_quality_count <- sum(viruses %in% CheckV_quality_all$virus_ID[CheckV_quality_all$quality == "Medium-quality"])
  return(data.frame(complete_quality_count, high_quality_count, medium_quality_count))
}

# Function to estimate summary stats of numeric variables
summary_stats <- function(vector) {
  summary <- c(
    mean = mean(vector),
    median = median(vector),
    q1 = quantile(vector, 0.25),
    q3 = quantile(vector, 0.75),
    min = min(vector),
    max = max(vector)
  )
  return(summary)
}

# Function to estimate the origin of genomes from each vOTU at a given taxonomic level (genus or family)
estimate_origin <- function(tax_level, clusters, databases) {
  results <- data.frame(matrix(ncol=1, nrow=(nrow(clusters))))
  colnames(results) <- paste0("vOTU_", toupper(tax_level), "_composition_extended")
  
  for (database in databases) {
    match_rows <- apply(clusters, 1, function(row) any(grepl(database, row, ignore.case = TRUE)))
    results[, paste0("vOTU_", toupper(tax_level), "_composition_extended")][match_rows] <-
      paste0(results[, paste0("vOTU_", toupper(tax_level), "_composition_extended")][match_rows], "&", database)
  }
  results[, paste0("vOTU_", toupper(tax_level), "_composition_extended")] <-
    substr(results[, paste0("vOTU_", toupper(tax_level), "_composition_extended")], 4, nchar(
      results[, paste0("vOTU_", toupper(tax_level), "_composition_extended")]
    ))
  # Add vOTU representative
  results$Virus_ID <- clusters$V1
  
  return(results)
}

# Function to generate combination of DBs for Upset plot
generate_DB_combinations <- function(databases, prefix = "", index = 1) {
  combinations <- c()
  if (index == length(databases)) {
    return(paste0(prefix, databases[index]))
  }
  combinations <- c(combinations, paste0(prefix, databases[index]))
  for (i in (index + 1):length(databases)) {
    sub_combinations <- generate_combinations(databases, paste0(prefix, databases[index], "&"), i)
    combinations <- c(combinations, sub_combinations)
  }
  return(combinations)
}

# Function to get individual-level presence of vOTUs
get_individual_presence <- function(abundance_table, metadata, virus_ids) {
  present_list <- lapply(virus_ids, function(v) {
    samples_with_votu <- colnames(abundance_table)[abundance_table[v, ] > 0]
    individuals <- unique(metadata$NEXT_ID[metadata$NG_ID %in% samples_with_votu])
    return(individuals)
  })
  return(present_list)
}

# Function to get NCBI taxonomy from taxid
get_taxonomy <- function(taxid) {
  raw_xml <- entrez_fetch(db = "taxonomy", id = taxid, rettype = "xml", parsed = FALSE)
  doc <- read_xml(raw_xml)
  lineage <- xml_find_all(doc, ".//LineageEx/Taxon")
  
  ranks <- xml_text(xml_find_all(lineage, "Rank"))
  names <- xml_text(xml_find_all(lineage, "ScientificName"))
  taxonomy <- setNames(names, tolower(ranks))  # Lowercase all rank names for consistency
  
  species <- xml_text(xml_find_first(doc, ".//ScientificName"))
  taxonomy["species"] <- species
  
  # Include "realm" in desired order
  desired_order <- c("realm", "superkingdom", "kingdom", "phylum", 
                     "class", "order", "family", "genus", "species")
  
  taxonomy_df <- as.data.frame(
    as.list(taxonomy[desired_order[desired_order %in% names(taxonomy)]]), 
    stringsAsFactors = FALSE
  )
  return(taxonomy_df)
}


##************************************************************************
# 1. Load metadata table for the LL-NEXT samples 
#*************************************************************************
# Load the metadata table 
# Here, we generate the final metadata table excluding samples with >75% unclassified in metaphlan 4 
# (unclassified group included relative abundances of two known contaminants)
Sample_metadata <- read.table("Metadata_NEXT/LLNEXT_metadata_03_03_2023.txt", header = T) 
Sample_metadata <- Sample_metadata [Sample_metadata$metaphlan4_unclassified_with_contaminants < 75,]

# Exclude previously identified outliers from the metadata and maternal samples from M1 and M2
Sample_metadata <- Sample_metadata[Sample_metadata$NG_ID != "AMBF025221E2",]

Sample_metadata <- Sample_metadata[!(Sample_metadata$Type == "mother"
                                     & Sample_metadata$Timepoint_categorical %in% c("M1", "M2")), ]

# Add total clean read numbers
Clean_reads <- read.delim("Metadata_NEXT/LLNEXT_nreads.txt", header = T) 
Clean_reads <- Clean_reads[Clean_reads$sample %in% Sample_metadata$NG_ID, ]
Sample_metadata$read_depth <- Clean_reads$clean_reads

# Create infant and maternal metadata
Sample_metadata_infants <- Sample_metadata[Sample_metadata$Type == "infant", ]
Sample_metadata_mothers <- Sample_metadata[Sample_metadata$Type == "mother", ]

##************************************************************************
# 2. Process abundance table
#*************************************************************************
# Read abundance table
Abundance_table <- read.delim("Abundance_table/LLNEXT_Abundance_Table_RPKM.txt") 

# Select only the samples from the metadata
Abundance_table <- Abundance_table[, colnames(Abundance_table) %in% Sample_metadata$NG_ID]

# Remove from the table those viruses not present in any sample
Absent_viruses <- rownames(Abundance_table[rowSums(Abundance_table)==0, ]) 
Abundance_table <- Abundance_table[!rownames(Abundance_table) %in% Absent_viruses,]

# Remove also potential contaminants ("LLNEXT_25294", "LLNEXT_79800")
Abundance_table <- Abundance_table[
  !rownames(Abundance_table) %in% c("LLNEXT_25294", "LLNEXT_79800"),
]

Present_viruses <- rownames(Abundance_table)
  
# Reorder abundance table to match Final metadata
Abundance_table <- Abundance_table[, Sample_metadata$NG_ID]

##************************************************************************
# 3. Generation of metadata table for each vOTU representative
#*************************************************************************

#**********************************
#A. Add DB of origin and virus IDs
#**********************************
Virus_metadata <- read.delim("Metadata_NEXT/Dereplication_input_sequences_nodup_DB_origin.txt", header=F, col.names = c("Virus_ID_original", "DB", "Virus_ID"))
Virus_metadata <- Virus_metadata %>%
  filter(Virus_ID %in% Present_viruses)
rownames(Virus_metadata) <- Virus_metadata$Virus_ID

#**********************************
#B. Add sequence length and GC content
#**********************************
# Read file and merge it with viral metadata
Length_GC <- read.delim("1_GENERAL_STATISTICS_AND_METADATA/All_viruses/Length_GC_rep_seqs.txt", header=F, col.names = c("Virus_ID", "Length", "GC")) 
Virus_metadata <- left_join(Virus_metadata, Length_GC, by = "Virus_ID")

#**********************************
#C. Add CheckV quality
#**********************************
# Estimate the CheckV quality for all viral genomes used as input for dereplication. Reason for this:
#- Quality/Completeness is missing for viruses from external DBs
# Load file and add quality information to non-LLNEXT viruses
CheckV_quality_all <- read.delim("1_GENERAL_STATISTICS_AND_METADATA/All_viruses/quality_summary.tsv", header=T)
CheckV_quality_all <- CheckV_quality_all[,c("contig_id", "viral_genes", "checkv_quality", "completeness_method")]
colnames(CheckV_quality_all) <- c("Virus_ID", "Viral_genes", "Quality", "Completeness_method")

Virus_metadata <- left_join(Virus_metadata, CheckV_quality_all, by = "Virus_ID")

#**********************************
#D. Add number of genes predicted by prodigal-gv
#**********************************\
vOTU_sequences <- readAAStringSet("1_GENERAL_STATISTICS_AND_METADATA/All_viruses/Present_vOTUs_proteins.faa")

# Extract Virus IDs
Virus_ids <- gsub(" #.*", "", names(vOTU_sequences))  
Virus_ids <- gsub("_[^_]*$", "", Virus_ids, perl = TRUE) 

# Count occurrences of each Virus_ID
Virus_proteins_counts <- table(Virus_ids)

# Convert to a data frame
Virus_proteins <- data.frame("Virus_ID" = names(Virus_proteins_counts), 
                             "n_genes" = as.numeric(Virus_proteins_counts), 
                             stringsAsFactors = FALSE)

# Join with Virus metadata
Virus_metadata <- left_join(Virus_metadata, Virus_proteins, by = "Virus_ID")

#**********************************
#E. Add the DB of origin of the genomes within each vOTUs
#**********************************
# Read file with vOTU clusters
vOTU_clusters <- read.delim("1_GENERAL_STATISTICS_AND_METADATA/All_viruses/viral_clusters.tsv", header = F)
vOTU_clusters_present <- vOTU_clusters[vOTU_clusters$V1 %in% Present_viruses,]
vOTU_clusters_present_processed <-  data.frame(vOTU_clusters_present[,-2], 
                                               tstrsplit(as.character(vOTU_clusters_present[,2]), ",", fixed = TRUE))
colnames(vOTU_clusters_present_processed) [1] <- c("rep_seq")

#############################
#E.1. All DBs from each vOTU
#############################
databases <- c("LLNEXT","GPD", "MGV", "IMG_VR", "ELGV","Shah","Benler","RefSeq", "Gulyaeva", "Guerin","Yutin")

# Estimate the origin of genomes from each vOTU
results <- data.frame(matrix(ncol=1, nrow=(nrow(vOTU_clusters_present_processed))))
colnames(results) <- "vOTU_DB_composition_extended"

for (database in databases) {
  match_rows <- apply(vOTU_clusters_present_processed, 1, function(row) any(grepl(database, row, ignore.case = TRUE)))
  results$vOTU_DB_composition_extended[match_rows] <- paste0(results$vOTU_DB_composition_extended[match_rows], "&", database)
}
results$vOTU_DB_composition_extended <- substr(results$vOTU_DB_composition_extended, 4, nchar(results$vOTU_DB_composition_extended))
results$Virus_ID <- vOTU_clusters_present_processed$rep_seq

# Merge it with viral metadata
Virus_metadata <- left_join(Virus_metadata, results, by = "Virus_ID")

#############################
#E.2. vOTUs with only LLNEXT genomes
#############################
# Estimate which vOTUs have viral genomes that do not match "LLNEXT" pattern (from other DBs)
vOTUs_with_external_genomes <- apply(vOTU_clusters_present_processed, 1, function(row) {
  any(!is.na(row) & !grepl("LLNEXT", row, ignore.case = TRUE))
})

# Summarize results
vOTU_composition <- data.frame("Virus_ID" = vOTU_clusters_present_processed$rep_seq,
                                  "vOTU_DB_composition" = ifelse(vOTUs_with_external_genomes, "OTHER_DBs", "LLNEXT"))

# Merge it with viral metadata
Virus_metadata <- left_join(Virus_metadata, vOTU_composition, by = "Virus_ID")

#**********************************
#F. Taxonomic assignment
#**********************************

###################
#F.1. RefSeq taxonomy
###################

# Identify vOTUs that contain RefSeq genomes based on metadata
vOTUs_RefSeq <- Virus_metadata[grep("RefSeq", Virus_metadata$vOTU_DB_composition_extended), "Virus_ID"]
colnames(vOTU_clusters_present_processed)[1] <- "vOTU"

# Generate a df (long-format) with vOTU representative - RefSeq genome matches
vOTU_rep_RefSeq_genome <- pivot_longer(vOTU_clusters_present_processed,
                        cols = -vOTU,
                        names_to = NULL,
                        values_to = "genome") %>%
  filter(!is.na(genome)) %>%
  filter(grepl("RefSeq", genome, ignore.case = TRUE)) %>%
  distinct()

# Add original RefSeq IDs and accession numbers
all_genome_names <- read.delim("Metadata_NEXT/Dereplication_input_sequences_nodup_DB_origin.txt", header = F,
                               col.names = c("Virus_ID_original", "DB", "Virus_ID"))

vOTU_rep_RefSeq_genome <- vOTU_rep_RefSeq_genome %>%
  left_join(all_genome_names %>% 
              select(Virus_ID, Virus_ID_original),
            by = c("genome" = "Virus_ID"))

vOTU_rep_RefSeq_genome$Accession_number <- sub(" .*", "", vOTU_rep_RefSeq_genome$Virus_ID_original)

# Fetch taxonomy
accessions <- vOTU_rep_RefSeq_genome$Accession_number
summaries <- entrez_summary(db = "nuccore", id = accessions)
taxids <- map_chr(summaries, ~ as.character(.x$taxid))# Extract taxids

# Apply function to all taxids
RefSeq_tax <- map_dfr(taxids, get_taxonomy)

# Add accession numbers for reference
RefSeq_tax$Accession <- accessions

# Rename and reorder columns
ordered_cols <- c("realm", "superkingdom", "kingdom", "phylum", "class", 
                  "order", "family", "genus", "species", "Accession")
ordered_present <- ordered_cols[ordered_cols %in% colnames(RefSeq_tax)]
RefSeq_tax <- RefSeq_tax[, ordered_present]
colnames(RefSeq_tax) <- gsub("(^[a-z])", "\\U\\1", colnames(RefSeq_tax), perl = TRUE)
RefSeq_tax$Virus_ID <- vOTU_rep_RefSeq_genome$vOTU[match(RefSeq_tax$Accession, vOTU_rep_RefSeq_genome$Accession_number)]
RefSeq_tax$Accession <- NULL

# Process taxonomy to select only lowest common taxonomic rank for Virus_ID (vOTUs) with different predictions
# Define taxonomy ranks in desired order of specificity
tax_ranks <- c("Realm", "Kingdom", "Phylum", "Class", 
               "Order", "Family", "Genus", "Species")

# Collapse taxonomy to lowest common rank per Virus_ID
RefSeq_tax <- RefSeq_tax %>%
  group_by(Virus_ID) %>%
  reframe(across(
    all_of(tax_ranks),
    ~ {
      vals <- unique(na.omit(.))
      if (length(vals) == 1) vals[1] else NA_character_
    }
  )) %>%
  filter(!is.na(Virus_ID))

###################
#F.2. VITAP taxonomy
###################
# Load VITAP taxonomic assignment
VITAP_tax <- read.delim("4_TAXONOMY/best_determined_lineages.tsv", 
                        sep = "\t", header = T, col.names = c("Virus_ID", "Lineage", 
                                                              "Score", "Confidence"))  
VITAP_tax <- VITAP_tax[,c("Virus_ID", "Lineage")]

# Split viral taxonomy
tax_levels <- strsplit(as.character(VITAP_tax$Lineage), ";")

# Determine the maximum number of levels in the taxonomy column and fill empty values with NA
max_levels <- max(sapply(tax_levels, length))
tax_levels <- lapply(tax_levels, function(x) c(x, rep(NA, max_levels - length(x))))

# Bind the taxonomy vectors to the original data frame
VITAP_tax  <- cbind(VITAP_tax, do.call("rbind", tax_levels))
VITAP_tax$Lineage <- NULL

# Rename the new columns
colnames(VITAP_tax)[2:9] <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm")

# Remove viruses introduced by the tool (not present in the samples)
VITAP_tax <- VITAP_tax[!(VITAP_tax$Virus_ID %in% c("BK013902.1", "AJ428555.1")), ]

# Remove [Species]_, [Genus]_, [Family]_, [Order]_, [Class]_, [Phylum]_, [Kingdom]_ and [Realm]_ prefixes
columns_to_clean <- c("Species", "Genus", "Family", "Order", "Class", "Phylum", "Kingdom", "Realm")
VITAP_tax <- VITAP_tax %>%
  mutate(across(all_of(columns_to_clean), ~ str_remove(., "^\\[(Species|Genus|Family|Order|Class|Phylum|Kingdom|Realm)\\]_")))

# Substitute "-" values with NA
VITAP_tax[VITAP_tax == "-"] <- NA

###################
#F.3. geNomad taxonomy
###################

# Load geNomad taxonomic assignment
geNomad_tax <- read.delim("4_TAXONOMY/Present_vOTUs_taxonomy.tsv", 
                          sep = "\t", header = T, col.names = c("Virus_ID", "Number_genes_with_taxonomy", 
                                                                "Agreement", "Taxid", "Lineage"))  
geNomad_tax <- geNomad_tax[,c("Virus_ID", "Lineage")]

# Split viral taxonomy
tax_levels <- strsplit(as.character(geNomad_tax$Lineage), ";")

# Determine the maximum number of levels in the taxonomy column and fill empty values with NA
max_levels <- max(sapply(tax_levels, length))
tax_levels <- lapply(tax_levels, function(x) c(x, rep(NA, max_levels - length(x))))

# Bind the taxonomy vectors to the original data frame
geNomad_tax  <- cbind(geNomad_tax, do.call("rbind", tax_levels))
geNomad_tax$Lineage <- NULL

# Rename the new columns
colnames(geNomad_tax)[2:8] <- c("Superkingdom", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family")

# As there are some viruses without order assigned, some family-level taxonomy is assigned as class
# Move this taxonomic level to the correct column
assign_to_family <- grep("viridae", geNomad_tax$Order, value=T)
geNomad_tax$Family[grep("viridae", geNomad_tax$Order)] <- assign_to_family
geNomad_tax$Order[grep("viridae", geNomad_tax$Order)] <- NA

# Remove Superkingdom column
geNomad_tax$Superkingdom <- NULL

# Substitute empty values with NA
geNomad_tax[geNomad_tax == ""] <- NA

###################
#F.4. Combined taxonomy
###################

# Generated a combined taxonomy by using the following rules: 
# -- Add RefSeq taxonomy for vOTUs clustering with vOTU genomes (if multiple, select lowest common taxonomic rank)
# --If only one source has taxonomy -> Use it.
# --If both are available:
# Use always VITAP, unless both agree till Class level and geNomad offers more complete taxonomy (then use geNomad)

# Ensure all 3 taxonomies (geNomad, VITAP and RefSeq) have the same structure
tax_ranks <- c("Species","Genus","Family","Order","Class","Phylum","Kingdom","Realm")
RefSeq_ext <- RefSeq_tax %>% 
  select(Virus_ID, all_of(tax_ranks))

VITAP_ext <- VITAP_tax %>%
  select(Virus_ID, all_of(tax_ranks)) %>%
  rename_with(~ paste0(.x, "_vitap"), all_of(tax_ranks))

geNomad_ext <- geNomad_tax %>%
  mutate(Species = NA_character_, Genus = NA_character_) %>%  # fill missing levels
  select(Virus_ID, all_of(tax_ranks)) %>%
  rename_with(~ paste0(.x, "_genomad"), all_of(tax_ranks))

# Join VITAP and geNomad in a single data frame
combined_tax <- VITAP_ext %>%
  full_join(geNomad_ext, by = "Virus_ID") 

# Define function to check if VITAP and geNomad agree up to Class
agree_up_to_class <- function(row) {
  levels <- c("Kingdom", "Phylum", "Class")
  vitap <- row[paste0(levels, "_vitap")]
  genomad <- row[paste0(levels, "_genomad")]
  all(mapply(function(v, g) {
    if (is.na(v) || is.na(g)) return(TRUE)
    v == g
  }, vitap, genomad))
}

# Resolve final taxonomy
final_tax <- combined_tax %>%
  rowwise() %>%
  mutate(
    vitap_depth   = sum(!is.na(pick(paste0(tax_ranks, "_vitap")))),
    genomad_depth = sum(!is.na(pick(paste0(tax_ranks, "_genomad")))),
    agree_class   = agree_up_to_class(cur_data()),
    use_vitap = case_when(
      all(is.na(pick(ends_with("_genomad")))) ~ TRUE,
      all(is.na(pick(ends_with("_vitap")))) ~ FALSE,
      agree_class & genomad_depth > vitap_depth ~ FALSE,
      TRUE ~ TRUE
    ),
    Species = ifelse(use_vitap, Species_vitap, Species_genomad),
    Genus   = ifelse(use_vitap, Genus_vitap, Genus_genomad),
    Family  = ifelse(use_vitap, Family_vitap, Family_genomad),
    Order   = ifelse(use_vitap, Order_vitap, Order_genomad),
    Class   = ifelse(use_vitap, Class_vitap, Class_genomad),
    Phylum  = ifelse(use_vitap, Phylum_vitap, Phylum_genomad),
    Kingdom = ifelse(use_vitap, Kingdom_vitap, Kingdom_genomad),
    Realm   = ifelse(use_vitap, Realm_vitap, Realm_genomad)
  ) %>%
  ungroup() %>%
  select(Virus_ID, all_of(tax_ranks))


# Remove rows with only NAs
final_tax <- final_tax[rowSums(is.na(final_tax[ , -1])) != length(tax_ranks), ]

# Add RefSeq taxonomy when available
final_tax <- final_tax %>%
  filter(!Virus_ID %in% RefSeq_ext$Virus_ID) %>%
  bind_rows(RefSeq_ext) 

# Merge it with viral metadata
Virus_metadata <- left_join(Virus_metadata, final_tax, by = "Virus_ID")

# We identify 4 vOTUs assigned to Sinsheimervirus phiX174 species.
# We identify 6 vOTUs assigned to Riboviria realm. 
# We identify 2 vOTUs predicted to infect Phyllobacterium (by host assignment):LLNEXT_25294 and LLNEXT_79800
# We exclude them for our analysis (also from the Phage_IDs)
phiX174_vOTUs <- Virus_metadata$Virus_ID[which(Virus_metadata$Species == "Sinsheimervirus phiX174")]
Riboviria_vOTUs <- Virus_metadata$Virus_ID[which(Virus_metadata$Realm == "Riboviria")]
Virus_metadata <- Virus_metadata[-c(which(Virus_metadata$Species == "Sinsheimervirus phiX174")), ]
Virus_metadata <- Virus_metadata[-c(which(Virus_metadata$Realm == "Riboviria")), ]
Virus_metadata <- Virus_metadata[-c(which(Virus_metadata$Virus_ID %in% c("LLNEXT_25294", "LLNEXT_79800"))), ]

#**********************************
#F. Lifestyle assignment
#**********************************
## Note: Prediction was run for all viruses (including eukaryotic viruses)
# As this prediction is meant for bacteriophages, we keep only results for them
# 4 vOTUs were identified as eukaryotic (geNomad taxonomy):ELGV_15176, ELGV_15356, ELGV_16040 and RefSeq_98

# Prokaryotic viruses were identified based on RefSeq/geNomad/VITAP taxonomy
# When no taxonomy (or only high-level taxonomy) was available, family-level clustering was used to determine eukaryotic/prokaryotic vOTUs
# Within the family-cluster FC_92 there were vOTUs assigned as prokaryotic and eukaryotic: genus-level clustering was considered (GC_1190 and GC_1464 --> phage, GC_1189 --> not enough evidence:unassigned)
# 5 vOTUs had no taxonomy and did not cluster (family-level) with vOTUs with known taxonomy: MGV_80884, MGV_40026, MGV_39652, GPD_3864 and GPD_81693
Phage_IDs <- unlist(read.table("1_GENERAL_STATISTICS_AND_METADATA/All_viruses/list_LLNEXT_phages.txt"))

# Read BACPHLIP lifestyle prediction output
lifestyle_viruses <- read.delim("3_LIFESTYLE/Present_vOTUs.fa.bacphlip", 
                                col.names = c("Virus_ID", "Virulent", "Temperate"))

# Remove predictions for non-phage viruses
lifestyle_viruses <- lifestyle_viruses[lifestyle_viruses$Virus_ID %in% Phage_IDs,]

# Generate final lifestyle prediction:
## A) Viruses with Temperate > 0.5 --> Temperate
## B) Viruses with Virulent > 0.5:
## B.1) If Completeness > 90% (CheckV) --> Virulent
## B.2) If Completeness < 90% (CheckV) --> NA
quality <- Virus_metadata[,c("Virus_ID", "Quality")]
quality_lifestyle <- inner_join(lifestyle_viruses, quality, by = "Virus_ID")

quality_lifestyle <- quality_lifestyle %>%
  mutate(Lifestyle = case_when(
    Temperate > 0.5 ~ "Temperate",
    Virulent > 0.5 & (Quality == "High-quality" | Quality == "Complete") ~ "Virulent",
  )) 
quality_lifestyle <- quality_lifestyle[, c("Virus_ID", "Lifestyle")]

# Merge it with viral metadata
Virus_metadata <- left_join(Virus_metadata, quality_lifestyle, by = "Virus_ID")

#**********************************
#G. DGR identification
#**********************************
# Note: DGR identification was run for all viruses (including eukaryotic viruses)
# DGRs were only found in phage genomes

# Read DGR identification results
DGR_target_results <- read.delim("5_DGR_ANALYSIS/All_RT_1k_DGR_target_detection_filtered.tsv", check.names = F)

# Estimate A mutation bias and filter those DGRs with less than 75% of Adenine mutation rate
DGR_target_results$`Max A/T bias` <- sapply(strsplit(DGR_target_results$`Adjusted bias (A;T;C;G)`, ";"), function(x) {
  counts <- as.numeric(x)
  A <- counts[1]
  T <- counts[2]
  max(A, T) / sum(counts)
})

DGR_target_results <- DGR_target_results[DGR_target_results$`Max A/T bias` >= 0.75,]
DGR_target_results$Virus_ID <- sub("_[^_]*$", "", sub(".*:", "", DGR_target_results[[1]]))
Virus_metadata$DGR <- ifelse(Virus_metadata$Virus_ID %in% DGR_target_results$Virus_ID, "Yes", "No") #2548

#**********************************
#H. ARG annotation
#**********************************
# Note: ARG identification was run for all viruses (including eukaryotic viruses)

# Read DGR identification results
ARG_results <- read.delim("7_ARG_ANALYSIS/RGI_CARD.txt", check.names = F)
ARG_results$Virus_ID <- sub(" #.*", "", ARG_results$ORF_ID)
ARG_results$Virus_ID <- sub("_[^_]*$", "", ARG_results$Virus_ID)
ARG_results$Contig <-NULL
  
Virus_metadata$ARG <- ifelse(Virus_metadata$Virus_ID %in% ARG_results$Virus_ID, "Yes", "No")

#**********************************
#I. Genus and family-level clustering
#**********************************

# Load the genus_clusters.txt and family_clusters.txt files
genus_clusters <- read.delim("1_GENERAL_STATISTICS_AND_METADATA/All_viruses/genus_clusters.txt", header = F)
family_clusters <- read.delim("1_GENERAL_STATISTICS_AND_METADATA/All_viruses/family_clusters.txt", header = F)

# Add genus and family cluster representative as a variable in Virus_metadata
Virus_metadata$Genus_cluster <- NA
Virus_metadata$Family_cluster <- NA

# Estimate the origin of genomes from each vOTU at genus and family-level (vOTU_DB_composition extended)
results_genus_level <- estimate_origin("genus", genus_clusters, databases)
results_family_level <- estimate_origin("family", family_clusters, databases)

# Add a variable vOTU_DB_composition indicating vOTUs specific to LLNEXT
vOTUs_genus_with_external_genomes <- apply(genus_clusters, 1, function(row) {
  any(nchar(row) > 0 & !grepl("LLNEXT", row, ignore.case = TRUE))
})
vOTUs_family_with_external_genomes <- apply(family_clusters, 1, function(row) {
  any(nchar(row) > 0 & !grepl("LLNEXT", row, ignore.case = TRUE))
})

results_genus_level$vOTU_composition_genus <- ifelse(vOTUs_genus_with_external_genomes, "OTHER_DBs", "LLNEXT")
results_family_level$vOTU_composition_family <- ifelse(vOTUs_family_with_external_genomes, "OTHER_DBs", "LLNEXT")

# Add genus and family level clusters to Virus metadata
rownames(genus_clusters) <- paste("GC", seq_along(genus_clusters[, 1]), sep = "_")
rownames(family_clusters) <- paste("FC", seq_along(family_clusters[, 1]), sep = "_")

Virus_metadata$Genus_cluster <- sapply(Virus_metadata$Virus_ID, function(virus_id) {
  row_match <- which(genus_clusters == virus_id, arr.ind = TRUE)
  if (length(row_match) > 0) {
    return(rownames(genus_clusters)[row_match[1, 1]])
  } else {
    return(NA)  
  }
})

Virus_metadata$Family_cluster <- sapply(Virus_metadata$Virus_ID, function(virus_id) {
  row_match <- which(family_clusters == virus_id, arr.ind = TRUE)
  if (length(row_match) > 0) {
    return(rownames(family_clusters)[row_match[1, 1]])
  } else {
    return(NA)  
  }
})

#**********************************
#J. Prevalence and mean abundance 
#**********************************
# Generate abundance tables for mothers and infants
Abundance_table_mothers <- Abundance_table[, Sample_metadata_mothers$NG_ID]
Abundance_table_mothers <- Abundance_table_mothers[rowSums(Abundance_table_mothers) > 0, ]
Abundance_table_mothers <- Abundance_table_mothers[rownames(Abundance_table_mothers) %in% Virus_metadata$Virus_ID, ] #19,433
Abundance_table_mothers_rel <- sweep(Abundance_table_mothers, 2, colSums(Abundance_table_mothers, na.rm = TRUE), FUN = "/") * 100

Abundance_table_infants <- Abundance_table[, Sample_metadata_infants$NG_ID]
Abundance_table_infants <- Abundance_table_infants[rowSums(Abundance_table_infants) > 0, ]
Abundance_table_infants <- Abundance_table_infants[rownames(Abundance_table_infants) %in% Virus_metadata$Virus_ID, ] # 16,805
Abundance_table_infants_rel <- sweep(Abundance_table_infants, 2, colSums(Abundance_table_infants, na.rm = TRUE), FUN = "/") * 100

# Estimate the prevalence (per individual) and the mean abundance of each vOTU
Virus_metadata$N_mothers_present <- 0
Virus_metadata$Prevalence_mothers <- 0
Virus_metadata$Mean_abundance_mothers <- 0
Virus_metadata$Mean_rel_abundance_mothers <- 0

Virus_metadata$N_infants_present <- 0
Virus_metadata$Prevalence_infants <- 0
Virus_metadata$Mean_abundance_infants <- 0
Virus_metadata$Mean_rel_abundance_infants <- 0

# Mothers
mother_vOTUs <- Virus_metadata$Virus_ID[Virus_metadata$Virus_ID %in% rownames(Abundance_table_mothers)]
mother_present_individuals <- get_individual_presence(Abundance_table_mothers, Sample_metadata_mothers, mother_vOTUs)

Virus_metadata$N_mothers_present[match(mother_vOTUs, Virus_metadata$Virus_ID)] <- lengths(mother_present_individuals)
Virus_metadata$Prevalence_mothers[match(mother_vOTUs, Virus_metadata$Virus_ID)] <- 100 * lengths(mother_present_individuals) / length(unique(Sample_metadata_mothers$NEXT_ID))
Virus_metadata$Mean_abundance_mothers[match(mother_vOTUs, Virus_metadata$Virus_ID)] <- rowMeans(Abundance_table_mothers[mother_vOTUs, , drop = FALSE], na.rm = TRUE)
Virus_metadata$Mean_rel_abundance_mothers[match(mother_vOTUs, Virus_metadata$Virus_ID)] <- rowMeans(Abundance_table_mothers_rel[mother_vOTUs, , drop = FALSE], na.rm = TRUE)


# Infants
infant_vOTUs <- Virus_metadata$Virus_ID[Virus_metadata$Virus_ID %in% rownames(Abundance_table_infants)]
infant_present_individuals <- get_individual_presence(Abundance_table_infants, Sample_metadata_infants, infant_vOTUs)

Virus_metadata$N_infants_present[match(infant_vOTUs, Virus_metadata$Virus_ID)] <- lengths(infant_present_individuals)
Virus_metadata$Prevalence_infants[match(infant_vOTUs, Virus_metadata$Virus_ID)] <- 100 * lengths(infant_present_individuals) / length(unique(Sample_metadata_infants$NEXT_ID))
Virus_metadata$Mean_abundance_infants[match(infant_vOTUs, Virus_metadata$Virus_ID)] <- rowMeans(Abundance_table_infants[infant_vOTUs, , drop = FALSE], na.rm = TRUE)
Virus_metadata$Mean_rel_abundance_infants[match(infant_vOTUs, Virus_metadata$Virus_ID)] <- rowMeans(Abundance_table_infants_rel[infant_vOTUs, , drop = FALSE], na.rm = TRUE)


##************************************************************************
# 4. Generation of final phage and viral metadata and abundance tables
#*************************************************************************
# After generating the metadata table for all vOTUs, we now generate separated tables for phages and eukaryotic viruses
Virus_metadata_phages <- Virus_metadata[Virus_metadata$Virus_ID %in% Phage_IDs,]
Virus_metadata_eukaryotic_unnasigned <- Virus_metadata[!Virus_metadata$Virus_ID %in% Phage_IDs,]

# We generate the final abundance tables for mothers and infants
Abundance_table <- Abundance_table[rownames(Abundance_table) %in% Virus_metadata$Virus_ID,]
Abundance_table <- Abundance_table[match(Virus_metadata$Virus_ID, rownames(Abundance_table)), ]
Abundance_table_infants <- Abundance_table[, colnames(Abundance_table) %in% Sample_metadata_infants$NG_ID]
Absent_viruses_infants <- rownames(Abundance_table_infants[rowSums(Abundance_table_infants)==0, ]) 
Abundance_table_infants <- Abundance_table_infants[!rownames(Abundance_table_infants) %in% Absent_viruses_infants,] #16,805
Present_viruses_infants <- rownames(Abundance_table_infants) #16,805
Abundance_table_infants <- Abundance_table_infants[, Sample_metadata_infants$NG_ID]

Abundance_table_mothers <- Abundance_table[, colnames(Abundance_table) %in% Sample_metadata_mothers$NG_ID]
Absent_viruses_mothers <- rownames(Abundance_table_mothers[rowSums(Abundance_table_mothers)==0, ]) 
Abundance_table_mothers <- Abundance_table_mothers[!rownames(Abundance_table_mothers) %in% Absent_viruses_mothers,] #19,435
Present_viruses_mothers <- rownames(Abundance_table_mothers) #19,435
Abundance_table_mothers <- Abundance_table_mothers[, Sample_metadata_mothers$NG_ID]

Abundance_table_phages <- Abundance_table[rownames(Abundance_table) %in% Phage_IDs,] #30,845
Abundance_table_phages <- Abundance_table_phages[match(Virus_metadata_phages$Virus_ID, rownames(Abundance_table_phages)), ]
Abundance_table_phages_mothers <- Abundance_table_mothers[rownames(Abundance_table_mothers) %in% Phage_IDs,] #19,337
Abundance_table_phages_infants <- Abundance_table_infants[rownames(Abundance_table_infants) %in% Phage_IDs,] #16,668


##************************************************************************
# 5. Calculation of summary statistics of phage vOTUs
#*************************************************************************

#Estimate summary stats for all vOTUs
table(Virus_metadata$DB) #DB distribution of vOTU representatives
table(Virus_metadata$Quality)  #Quality distribution of vOTU representatives
summary_stats(Virus_metadata$Length) #Summary stats of the length of all vOTU representatives

# Estimate the summary stats for vOTUs with only LLNEXT genomes 
table(Virus_metadata$Quality[Virus_metadata$vOTU_DB_composition=="LLNEXT"])
summary_stats(Virus_metadata$Length[Virus_metadata$vOTU_DB_composition=="LLNEXT"])
summary_stats(Virus_metadata$GC[Virus_metadata$vOTU_DB_composition=="LLNEXT"])

# Estimate the number of genus and family-level vOTUs with only LLNEXT genomes
# (vOTUs with only LLNEXT representatives at genus or family level and only LLNEXT genomes at species level)
length(which(results_genus_level$Virus_ID[results_genus_level$vOTU_composition_genus == "LLNEXT"] %in%
        Virus_metadata$Virus_ID[Virus_metadata$vOTU_DB_composition == "LLNEXT"]))
length(which(results_family_level$Virus_ID[results_family_level$vOTU_composition_family == "LLNEXT"] %in%
               Virus_metadata$Virus_ID[Virus_metadata$vOTU_DB_composition == "LLNEXT"]))


##************************************************************************
# 5. Generation of summary plots
#*************************************************************************

#############################
# Origin of vOTUs
#############################

# Generate UpSet plot at species, genus and family levels
DB_distribution_species <- table(Virus_metadata$vOTU_DB_composition_extended)
DB_distribution_genus <- table(results_genus_level$vOTU_GENUS_composition_extended)
DB_distribution_family <- table(results_family_level$vOTU_FAMILY_composition_extended)

pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/Upset_DBs_species_vOTUs.pdf', width=23, height=11)
UpSetR::upset(fromExpression(DB_distribution_species), 
      nintersects = 30, 
      nsets = 12, 
      order.by = "freq", 
      decreasing = T, 
      mb.ratio = c(0.45, 0.55),
      text.scale = 3.9, 
      point.size = 5, 
      line.size = 1, 
      sets.bar.color = "darkblue",
      shade.color = "grey",
      shade.alpha = 0.15, show.numbers = F,
      sets.x.label = "Number vOTUs",
)
dev.off()


#############################
# Length  distribution
#############################
# Generate histogram for length distribution of vOTU representatives
pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/Length_distribution_log.pdf', width=6, height=1.5)
ggplot(Virus_metadata, aes(x = Length)) +
  geom_histogram(bins = 40, fill = alpha("steelblue", 0.7)) +
  labs(x = "Length (kbp)", y = "Frequency") +  # Modify the axis labels as needed
  theme_classic() +
  scale_x_log10(labels = scales::number_format(scale = 0.001)) +  
  theme(
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 13),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "top"
  )
dev.off()


#############################
# Distribution of viral taxonomy and lifestyle
#############################

# Change all values other than "Caudoviricetes" to "Other"
Viral_class <- Virus_metadata$Class[!is.na(Virus_metadata$Class)]
Viral_class[Viral_class != "Caudoviricetes"] <- "Other"
Viral_class_distribution <- table(Viral_class)

# Generate a dataframe with the counts and proportions
Viral_class_distribution <- data.frame(
  Class = names(Viral_class_distribution),
  Count = as.numeric(Viral_class_distribution) 
)
Viral_class_distribution$Description <- rep ("Viral class distribution", length(Viral_class_distribution$Class))
Viral_class_distribution <- Viral_class_distribution %>%
  mutate(Proportion = Count / sum(Count))

Viral_class_summary <- Viral_class_distribution %>%
  group_by(Class) %>%
  summarise(Count = sum(Count, na.rm = TRUE), .groups = "drop") %>%
  mutate(Proportion = Count / sum(Count))

palette <- c("#ADD8E6", "#E6E6FA")

# Generate plot
pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/Viral_classes_piechart.pdf', width=5, height=2.5)
ggplot(Viral_class_summary, aes(x = 1, y = Proportion, fill = Class)) +
  geom_col(width = 1) + 
  coord_polar(theta = "y") +
  scale_fill_manual(values = palette, labels = c("Caudoviricetes", "Other")) +
  labs(x = NULL, y = NULL, fill = "Class") +
  theme_void() +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.position = "right",
    legend.key.height = unit(1.4, "lines")
  )
dev.off()

# Generate a pie chart with the predicted viral lyfestyle
Viral_lifestyle_distribution <- table(Virus_metadata$Lifestyle, useNA = "ifany")

# Generate a dataframe with the counts and proportions
Viral_lifestyle_distribution <- data.frame(
  Lifestyle = names(Viral_lifestyle_distribution),
  Count = as.numeric(Viral_lifestyle_distribution) 
)
Viral_lifestyle_distribution$Description <- rep ("Viral lifestyle distribution", length(Viral_lifestyle_distribution$Lifestyle))
Viral_lifestyle_distribution <- Viral_lifestyle_distribution %>%
  mutate(Proportion = 100*(Count / sum(Count)))

Virus_metadata$Lifestyle[is.na(Virus_metadata$Lifestyle)] <- "Not-determined"
Virus_metadata$Lifestyle <- factor(Virus_metadata$Lifestyle, levels = c("Temperate", "Virulent",
                                                                    "Not-determined"))
# Generate a pie chart
palette <- c("#E0C9A8", "#C8D9BF", "#C0C0C0")

pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/Lifestyle_piechart.pdf', width=5, height=2.5)
ggplot(Virus_metadata, aes(x = "", fill = Lifestyle)) +
  geom_bar(size = 0.6) +
  coord_polar(theta = "y") +
  labs(x = "Number of genomes vOTU", y = NULL, fill = "Lifestyle") +
  theme_void() +
  scale_fill_manual(values = palette, labels = c("Temperate", "Virulent", "Not-determined")) +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.position = "right",
    legend.key.height = unit(1.4, "lines")
  )
dev.off()

#############################
# Completeness distribution
#############################

Virus_metadata$Quality <- factor(Virus_metadata$Quality, levels = c("Complete", "High-quality",
                                                                                "Medium-quality","Low-quality",
                                                                                "Not-determined"))
# Generate a pie chart
palette <- c("#D8BFD8", "#C0C0C0", "#B0E0E6", "#E6E6FA", "#F0E68C")

pdf('1_GENERAL_STATISTICS_AND_METADATA/Plots/Completeness_piechart.pdf', width=5, height=2.5)
ggplot(Virus_metadata, aes(x = "", fill = Quality)) +
  geom_bar(size=0.6) +
  coord_polar(theta = "y") +
  labs(x = "Number of genomes vOTU", y = NULL, fill = "Completeness distribution") +
  theme_void() +
  scale_fill_manual(values = palette, labels = c("Complete", "High quality",
                                                 "Medium quality","Low quality",
                                                 "Not determined")) +
  theme(legend.title = element_text(size = 18),
        legend.text = element_text(size = 15),
        legend.position = "right",
        legend.key.height = unit(1.4, "lines")  
  )
dev.off()

#****************
# Write results
#****************

# Abundance table and list of present and excluded viruses (with no presence in LLNEXT samples)
write.table(Abundance_table,"Abundance_table/LLNEXT_Viral_Abundance_Table_07052025.txt", sep = "\t", row.names = T, quote = FALSE)
write.table(Present_viruses,"Abundance_table/Present_viruses.txt", sep = "\n", row.names = F, col.names = F, quote = FALSE)
write.table(Abundance_table_infants,"Abundance_table/LLNEXT_Viral_Abundance_Table_INFANTS_07052025.txt", sep = "\t", row.names = T, quote = FALSE)
write.table(Present_viruses_infants,"Abundance_table/Present_viruses_INFANTS.txt", sep = "\n", row.names = F, col.names = F, quote = FALSE)
write.table(Abundance_table_mothers,"Abundance_table/LLNEXT_Viral_Abundance_Table_MOTHERS_07052025.txt", sep = "\t", row.names = T, quote = FALSE)
write.table(Present_viruses_mothers,"Abundance_table/Present_viruses_MOTHERS.txt", sep = "\n", row.names = F, col.names = F, quote = FALSE)

# Abundance tables of phages
write.table(Abundance_table_phages,"Abundance_table/LLNEXT_Phage_Abundance_Table_07052025.txt", sep = "\t", row.names = T, quote = FALSE)
write.table(Abundance_table_phages_mothers,"Abundance_table/LLNEXT_Phage_Abundance_Table_MOTHERS_07052025.txt", sep = "\t", row.names = T, quote = FALSE)
write.table(Abundance_table_phages_infants,"Abundance_table/LLNEXT_Phage_Abundance_Table_INFANTS_07052025.txt", sep = "\t", row.names = T, quote = FALSE)

# Sample metadata
write.table(Sample_metadata,"Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 

# Virus metadata
write.table(Virus_metadata,"Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Virus_metadata_phages,"Metadata_NEXT/LLNEXT_Phage_Metadata_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
write.table(Virus_metadata_eukaryotic_unnasigned,"Metadata_NEXT/LLNEXT_Viral_Metadata_Euk_NA_07052025.txt", sep = "\t", row.names = F, quote = FALSE) 
