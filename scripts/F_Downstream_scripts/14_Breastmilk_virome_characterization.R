################################################################################
##### LL-NEXT: Breastmilk results processing
### Author(s): Asier Fernández-Pato
### Last updated: 8th December, 2025
################################################################################

#****************
# Load libraries
#****************
library(dplyr)
library(data.table)
library(vegan)
library(rentrez)
library(xml2)
library(UpSetR)
library(ggplot2)
library(stringr)
library(tidyr)
library(purrr)

#****************
# Define functions
#****************

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

# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_Analysis/VIRUSES/BREASTMILK/")


##************************************************************************
# 1. Process abundance results for the BM LL-NEXT samples 
#*************************************************************************
# Load the metadata table 
Sample_metadata <- read.table("Metadata_BM_VG_Gut_20_05_2025.txt", header = T) 
Sample_metadata <- Sample_metadata [Sample_metadata$Type =="Milk",]

# Read abundance table
Abundance_table <- read.delim("BM_LLNEXT_Abundance_Table_RPKM.txt") 

# Select only the samples from the metadata
Abundance_table <- Abundance_table[, colnames(Abundance_table) %in% Sample_metadata$NG_ID]

# Remove from the table those viruses not present in any sample
Absent_viruses <- rownames(Abundance_table[rowSums(Abundance_table)==0, ]) 
Abundance_table <- Abundance_table[!rownames(Abundance_table) %in% Absent_viruses,] #250
Present_viruses <- rownames(Abundance_table) #250


##************************************************************************
# 2. Generation of metadata table for each vOTU representative
#*************************************************************************

#**********************************
#A. Add DB of origin and virus IDs
#**********************************
Virus_metadata <- read.delim("Dereplication_input_sequences_nodup_DB_origin.txt", header=F, col.names = c("Virus_ID_original", "DB", "Virus_ID"))
Virus_metadata <- Virus_metadata %>%
  filter(Virus_ID %in% Present_viruses)
rownames(Virus_metadata) <- Virus_metadata$Virus_ID

#**********************************
#B. Add sequence length and GC content
#**********************************
# Read file and merge it with viral metadata
Length_GC <- read.delim("Length_GC.txt", header=F, col.names = c("Virus_ID", "Length", "GC")) 
Virus_metadata <- left_join(Virus_metadata, Length_GC, by = "Virus_ID")

#**********************************
#C. Add CheckV quality
#**********************************
# Load the CheckV quality 
CheckV_quality_all <- read.delim("quality_summary.tsv", header=T)
CheckV_quality_all <- CheckV_quality_all[,c("contig_id", "viral_genes", "checkv_quality", "completeness_method")]
colnames(CheckV_quality_all) <- c("Virus_ID", "Viral_genes", "Quality", "Completeness_method")

Virus_metadata <- left_join(Virus_metadata, CheckV_quality_all, by = "Virus_ID")

#**********************************
#D. Add the DB of origin of the genomes within each vOTUs
#**********************************
# Read file with vOTU clusters
vOTU_clusters <- read.delim("viral_clusters.tsv", header = F)
vOTU_clusters_present <- vOTU_clusters[vOTU_clusters$V1 %in% Present_viruses,]
vOTU_clusters_present_processed <-  data.frame(vOTU_clusters_present[,-2], 
                                               tstrsplit(as.character(vOTU_clusters_present[,2]), ",", fixed = TRUE))
colnames(vOTU_clusters_present_processed) [1] <- c("rep_seq")

#############################
#D.1. All DBs from each vOTU
#############################
databases <- c("BM_LLNEXT", "LLNEXT","GPD", "MGV", "IMG_VR", "ELGV","Shah","Benler","RefSeq", "Gulyaeva", "Guerin","Yutin")

# Estimate the origin of genomes from each vOTU
results <- data.frame(matrix(ncol=1, nrow=(nrow(vOTU_clusters_present_processed))))
colnames(results) <- "vOTU_DB_composition_extended"

for (database in databases) {
  pattern <- paste0("^", database, "_") 
  match_rows <- apply(vOTU_clusters_present_processed, 1, function(row) {
    any(sapply(row, function(x) grepl(pattern, x, ignore.case = TRUE)))
  })
  
  results$vOTU_DB_composition_extended[match_rows] <- paste0(
    results$vOTU_DB_composition_extended[match_rows], "&", database
  )
}
results$vOTU_DB_composition_extended <- substr(results$vOTU_DB_composition_extended, 4, nchar(results$vOTU_DB_composition_extended))
results$Virus_ID <- vOTU_clusters_present_processed$rep_seq

# Merge it with viral metadata
Virus_metadata <- left_join(Virus_metadata, results, by = "Virus_ID")

# Estimate number of vOTUs exclusive to BM_LLNEXT samples
length(Virus_metadata$Virus_ID[Virus_metadata$vOTU_DB_composition_extended == "BM_LLNEXT"]) # 152
table(Virus_metadata$Quality)


#**********************************
#E. Taxonomic assignment
#**********************************

###################
#E.1. RefSeq taxonomy
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
all_genome_names <- read.delim("Dereplication_input_sequences_nodup_DB_origin.txt", header = F,
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
#E.2. VITAP taxonomy
###################
# Load VITAP taxonomic assignment
VITAP_tax <- read.delim("best_determined_lineages.tsv", 
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
#E.3. geNomad taxonomy
###################

# Load geNomad taxonomic assignment
geNomad_tax <- read.delim("Present_viruses_BM_taxonomy.tsv", 
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
#E.4. Combined taxonomy
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

#**********************************
#F. Lifestyle assignment
#**********************************
## Note: Prediction was run for all viruses (including eukaryotic viruses)
# As this prediction is meant for bacteriophages, we keep only results for them
# We identified 3 eukaryotic vOTUs (based on combined tax) (#ELGV_19307 no tax assigned - included here)

Euk_IDs <- c("ELGV_15176", "ELGV_16040", "RefSeq_98")

# Read BACPHLIP lifestyle prediction output
lifestyle_viruses <- read.delim("Present_viruses_BM.fa.bacphlip", 
                                col.names = c("Virus_ID", "Virulent", "Temperate"))

# Remove predictions for non-phage viruses
lifestyle_viruses <- lifestyle_viruses[!lifestyle_viruses$Virus_ID %in% Euk_IDs,]

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
#G. Host assignment
#**********************************

#################
# G.1. Default DB
#################

#A) Read iPHOP prediction output at genus level with default DB
viral_hosts_df <- read.delim("Host_prediction_to_genus_m90.csv", 
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
# - Predictions for 222/250 vOTUs

#################
# G.2. Enriched DB
#################

#B) Read iPHOP prediction output at genus level with enriched DB
viral_hosts_enr <- read.delim("Host_prediction_to_genus_m90_enriched.csv", 
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
# - Predictions for 229/250 vOTUs

#################
# G.3. Combined results
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
Virus_metadata <- left_join(Virus_metadata, viral_hosts_combined_to_merge, by = "Virus_ID")

# Set host columns as a factor (ordering them by occurence)
Virus_metadata$Bacterial_family_host <- factor(Virus_metadata$Bacterial_family_host)
Virus_metadata$Bacterial_genus_host <- factor(Virus_metadata$Bacterial_genus_host)

# Set to NA host prediction for the 3 eukaryotic vOTUs
Virus_metadata$Bacterial_family_host[Virus_metadata$Virus_ID %in% Euk_IDs] <- NA
Virus_metadata$Bacterial_genus_host[Virus_metadata$Virus_ID %in% Euk_IDs] <- NA

#**********************************
#H. Prevalence and mean abundance
#**********************************

# Estimate the prevalence (per individual) and the mean abundance of each vOTU
Virus_metadata$N_present <- 0
Virus_metadata$Prevalence <- 0
Virus_metadata$Mean_abundance <- 0

vOTUs <- Virus_metadata$Virus_ID[Virus_metadata$Virus_ID %in% rownames(Abundance_table)]
present_individuals <- get_individual_presence(Abundance_table, Sample_metadata, vOTUs)

Virus_metadata$N_present[match(vOTUs, Virus_metadata$Virus_ID)] <- lengths(present_individuals)
Virus_metadata$Prevalence[match(vOTUs, Virus_metadata$Virus_ID)] <- 100 * lengths(present_individuals) / length(unique(Sample_metadata$NEXT_ID))
Virus_metadata$Mean_abundance[match(vOTUs, Virus_metadata$Virus_ID)] <- rowMeans(Abundance_table[vOTUs, , drop = FALSE], na.rm = TRUE)

#############################
# PLOTS
#############################

#############################
# Origin of vOTU representatives
#############################

# Generate UpSet plot at species, genus and family levels
DB_distribution_species <- table(Virus_metadata$vOTU_DB_composition_extended)


pdf('BREASTMILK/Upset_DBs_species_vOTUs.pdf', width=23, height=11)
UpSetR::upset(fromExpression(DB_distribution_species), 
              nintersects = 30, 
              nsets = 12, 
              order.by = "freq", 
              decreasing = T, 
              mb.ratio = c(0.45, 0.55),
              text.scale = 2.5, 
              point.size = 3.5, 
              line.size = 1, 
              sets.bar.color = "darkblue",
              shade.color = "grey",
              shade.alpha = 0.15, show.numbers = F,
              sets.x.label = "Number vOTUs",
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

pdf('Completeness_piechart.pdf', width=5, height=2.5)
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

#############################
# Distribution of viral taxonomy 
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

pdf('Viral_classes_piechart.pdf', width=5, height=2.5)
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

#############################
# Bacterial host distribution
#############################

# Change all values not in the top 5 to Other
Virus_metadata <- Virus_metadata %>%
  mutate(Bacterial_genus_host = as.character(Bacterial_genus_host))

top5_hosts <- Virus_metadata %>%
  filter(!is.na(Bacterial_genus_host)) %>%  
  count(Bacterial_genus_host, sort = TRUE) %>%
  top_n(5, n) %>%
  pull(Bacterial_genus_host) 

Virus_metadata_plots <- Virus_metadata %>%
  mutate(Bacterial_genus_host = ifelse(Bacterial_genus_host %in% top5_hosts,
                                       Bacterial_genus_host, "Other"))

Virus_metadata$Bacterial_genus_host <- factor(Virus_metadata$Bacterial_genus_host,
                                              levels = c("Acinetobacter", "Streptococcus",
                                                         "Staphylococcus", "Pseudomonas_E",
                                                         "Stenotrophomonas", "Other"))

# Generate a pie chart
palette <- c("#A6CEE3","grey75","#B2DF8A", 
             "#FDBF6F","#FB9A99","#FFFFB3") 

pdf('Bacterial_hosts_piechart.pdf', width=5, height=2.5)
ggplot(Virus_metadata_plots, aes(x = "", fill = Bacterial_genus_host)) +
  geom_bar(width = 1, stat = "count") +
  coord_polar(theta = "y") +
  labs(x = NULL, y = NULL, fill = "Bacterial host") +
  theme_void() +
  scale_fill_manual(values = palette) +
  theme(
    legend.title = element_text(size = 18),
    legend.text = element_text(size = 15),
    legend.position = "right",
    legend.key.height = unit(1.4, "lines")
  )
dev.off()

#############################
# Lifestyle
#############################

Virus_metadata$Lifestyle[is.na(Virus_metadata$Lifestyle)] <- "Not-determined"
Virus_metadata$Lifestyle <- factor(Virus_metadata$Lifestyle, levels = c("Temperate", "Virulent",
                                                                        "Not-determined"))
# Generate a pie chart
palette <- c("#E0C9A8", "#C8D9BF", "#C0C0C0")

pdf('Lifestyle_piechart.pdf', width=5, height=2.5)
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


##************************************************************************
# 3. Comparison of BM vs FC: Diversity per sample
#*************************************************************************

# Note that some BM samples do not have any vOTU present
samples_no_viruses <- colnames(Abundance_table)[colSums(Abundance_table) == 0] 

# Load abundance tables and metadatas
Abundance_table_BM <- Abundance_table
Sample_metadata_BM <- Sample_metadata
Sample_metadata_BM <- Sample_metadata_BM [Sample_metadata_BM$Type == "Milk",]
Abundance_table_FC <- read.delim("../Abundance_table/LLNEXT_Viral_Abundance_Table_MOTHERS_07052025.txt")
Sample_metadata_FC <- read.delim("../Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_07052025.txt")

#Select families with milk samples
milk_families <- Sample_metadata_BM %>%
  filter(Type == "Milk") %>%
  pull(FAMILY)

# Subset those families from FC results
Sample_metadata_FC <- Sample_metadata_FC[(Sample_metadata_FC$FAMILY %in% milk_families) &
    (Sample_metadata_FC$Timepoint_categorical == "B"), ]
Abundance_table_FC <- Abundance_table_FC[, colnames(Abundance_table_FC) %in% Sample_metadata_FC$NG_ID]
Abundance_table_FC <- Abundance_table_FC[rowSums(Abundance_table_FC) > 0, ]

# Calculate richness and shannon diversity for breastmilk (M1) and fecal samples (M3)
diversity_BM <- data.frame(
  Sample = colnames(Abundance_table_BM),
  Type = "Breastmilk",
  Shannon = apply(Abundance_table_BM, 2, function(x) diversity(x, index = "shannon")),
  Richness = apply(Abundance_table_BM, 2, function(x) sum(x > 0))
)

diversity_FC <- data.frame(
  Sample = colnames(Abundance_table_FC),
  Type = "Feces",
  Shannon = apply(Abundance_table_FC, 2, function(x) diversity(x, index = "shannon")),
  Richness = apply(Abundance_table_FC, 2, function(x) sum(x > 0))
)

# Combine BM and FC
diversity_all <- bind_rows(diversity_BM, diversity_FC)
pairing_info_BM <- Sample_metadata_BM %>% select(NG_ID, FAMILY)
pairing_info_FC <- Sample_metadata_FC %>% select(NG_ID, FAMILY)
diversity_all <- diversity_all %>%
  left_join(bind_rows(pairing_info_BM %>% mutate(Type="Breastmilk"),
                      pairing_info_FC %>% mutate(Type="Feces")),
            by = c("Sample"="NG_ID", "Type"))

# Convert to long format 
diversity_long <- diversity_all %>%
  pivot_longer(cols = c(Shannon, Richness), names_to = "Metric", values_to = "Value")

# Generate the plot
pdf('Richness_Shannon_BM_FC_comparison.pdf', width = 2.6, height = 3.2)
ggplot(subset(diversity_long, Metric == "Shannon"),
       aes(x = Type, y = Value)) +
  geom_jitter(alpha = 0.4, color = "lightgrey", size = 1.8, width = 0.2) +
  geom_line(aes(group = FAMILY), color = "gray70", alpha = 0.5) +
  geom_boxplot(width = 0.6, aes(fill = Type), outlier.color = NA, size = 0.9) +
  scale_color_manual(values = c("Breastmilk" = "#D8A7B1", "Feces"= "#A1887F")) +
  scale_fill_manual(values = c("Breastmilk" = "#D8A7B1", "Feces"= "#A1887F")) +
  labs(x = NULL, y = "Shannon diversity") +
  scale_y_continuous(breaks = c(0,1,2,3,4,5), expand = expansion(add = c(0.3, 1.1))) +
  theme_bw(base_size = 16) +
  theme(
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA, size = 0.75),
    axis.text = element_text(size = 14),
    axis.title = element_text(size = 16),
    strip.text = element_text(size = 14),
    legend.position = "none"
  )
dev.off()

  
# Add statistical test
# Keep only 1 row per FAMILY and TYPE
diversity_filtered <- diversity_long[ave(diversity_long$Type,
                                         diversity_long$FAMILY,
                                         diversity_long$Type,
                                         FUN = seq_along) == 1,]

# Keep only families (mothers) with both breastmilk and fecal sample available
families_complete <- diversity_filtered %>%
  group_by(FAMILY) %>%
  summarize(n_types = n_distinct(Type)) %>%
  filter(n_types == 2) %>%
  pull(FAMILY)

diversity_filtered_complete <- diversity_filtered %>%
  filter(FAMILY %in% families_complete)

# Generate data in wide format
diversity_wide <- diversity_filtered_complete %>%
  filter(Metric == "Shannon") %>%  
  select(FAMILY, Type, Value) %>%
  pivot_wider(names_from = Type, values_from = Value, 
              names_prefix = "Value_")
colnames(diversity_wide) <- c("FAMILY", "Breastmilk_value", "Feces_value")

#Perform paired Wilcoxon test (p-value = 1.67e-10)
wilcox.test(diversity_wide$Breastmilk_value,
            diversity_wide$Feces_value,
            paired = TRUE)


#****************
# Write results
#****************
# For strain sharing, we exclude from the abundance and metadata table potential contaminants
contaminants <- c("BM_LLNEXT_1", "BM_LLNEXT_103")
Virus_metadata <- Virus_metadata[!Virus_metadata$Virus_ID %in% contaminants,]
Abundance_table <-Abundance_table[!rownames(Abundance_table) %in% contaminants,]
Present_viruses <- Present_viruses[!Present_viruses %in% contaminants]

# Abundance table, metadata and list of present viruses 
write.table(Present_viruses, file = "Present_viruses.txt", 
            quote = FALSE, row.names = FALSE, col.names = FALSE)
write.table(Abundance_table,"BM_LLNEXT_Viral_Abundance_Table_29072025.txt",
            sep = "\t", row.names = T, quote = FALSE)
write.table(Virus_metadata,"BM_LLNEXT_Viral_Metadata_29072025.txt",
            sep = "\t", row.names = F, quote = FALSE) 
