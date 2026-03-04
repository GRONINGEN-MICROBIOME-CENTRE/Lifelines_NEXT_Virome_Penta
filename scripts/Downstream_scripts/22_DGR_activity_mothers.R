################################################################################
##### LL-NEXT:  DGR analysis -  Activity in MOTHERS
### Author(s): Asier Fernández-Pato
### Last updated: 18th December, 2025
################################################################################

#****************
# Load libraries
#****************
library(dplyr)
library(tidyr)
library(data.table)
library(stringr)
library(Biostrings)
library(ggplot2)
library(scales)


#****************
# Define functions
#****************

# Set working directory 
setwd("~/Desktop/PhD/Projects/Virome LL-Next/Final_Analysis/VIRUSES/")

#Load metadata and abundance table for the LL-NEXT samples 

# Read metadata tables 
Sample_metadata <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_07052025.txt")
Sample_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_INFANTS_07052025.txt")
Sample_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Sample_Metadata_MOTHERS_07052025.txt")

# Read abundance table
Abundance_table <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_07052025.txt")
Abundance_table_infants <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_INFANTS_07052025.txt")
Abundance_table_mothers <- read.delim("Abundance_table/LLNEXT_Viral_Abundance_Table_MOTHERS_07052025.txt")

# Read viral metadata tables
Viral_metadata <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_07052025.txt")
Viral_metadata_infants <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_INFANTS_07052025.txt")
Viral_metadata_mothers <- read.delim("Metadata_NEXT/LLNEXT_Viral_Metadata_MOTHERS_07052025.txt")


#***********************************************
# 1. DGR activity analysis: data preparation
#***********************************************

###################################
# 1.1. Selection of DGR+ persistent phages in mothers and metadata generation
###################################

# We select mother vOTUs classified as persistent and with DGRs (n=1,076)
DGR_persistent_vOTUs <- Viral_metadata_mothers[Viral_metadata_mothers$Persistent == "Yes" &
                                                 Viral_metadata_mothers$DGR == "Yes",
                                               "Virus_ID"]

# Subset abundance table to only DGR+ persistent vOTUs
subset_abundance <- Abundance_table_mothers[rownames(Abundance_table_mothers) %in% DGR_persistent_vOTUs, ]

# Identify all samples where any DGR+ persistent vOTU is detected
samples_with_DGR_vOTUs <- colnames(subset_abundance)[colSums(subset_abundance > 0) > 0]

# Identify unique mothers in which any of the DGR+ vOTUs are detected (n=697)
mothers_with_DGR_vOTUs <- unique(Sample_metadata_mothers$NEXT_ID[
  Sample_metadata_mothers$NG_ID %in% samples_with_DGR_vOTUs
])

# Generate a table with the mothers in which each of the DGR+ vOTUs is detected as persistent
pregnancy_timepoints <- c("P12", "P28", "B")
post_pregnancy_timepoints  <- c("M3")

persistent_info <- data.frame(
  vOTU = rownames(subset_abundance),
  n_persistent_mothers = NA_integer_,
  persistent_mothers = NA_character_,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(persistent_info))) {
  votu <- persistent_info$vOTU[i]
  positive_samples <- colnames(subset_abundance)[subset_abundance[votu, ] > 0] # samples
  positive_meta <- Sample_metadata_mothers[Sample_metadata_mothers$NG_ID %in% positive_samples, ] # metadata
  
  # Split by mother
  mothers_list <- split(positive_meta, positive_meta$NEXT_ID)
  persistent_mothers <- c()
  
  for (mother_id in names(mothers_list)) {
    mother_data <- mothers_list[[mother_id]]
    timepoints <- mother_data$Timepoint
    has_early <- any(timepoints %in% pregnancy_timepoints)
    has_late  <- any(timepoints %in% post_pregnancy_timepoints)
    
    if (has_early && has_late) {
      persistent_mothers <- c(persistent_mothers, mother_id)
    }
  }
  
  # Store results
  persistent_info$n_persistent_mothers[i] <- length(persistent_mothers)
  persistent_info$persistent_mothers[i] <- paste(persistent_mothers, collapse = ",")
}

# Reshape to have persistence information per mother
persistent_info_mother <- persistent_info %>%
  filter(n_persistent_mothers > 0) %>%
  separate_rows(persistent_mothers, sep = ",") %>%
  group_by(persistent_mothers) %>%
  summarise(
    n_persistent_vOTUs = n(),
    persistent_vOTUs = paste(vOTU, collapse = ","),
    .groups = "drop"
  ) %>%
  dplyr::rename(mother = persistent_mothers)


###################################
# 1.2. Generation of input files for activity analysis
###################################

# Prepare simplified sample metadata file with 3 columns
Sample_metadata_mothers_filt <- Sample_metadata_mothers[,c("NG_ID", "NEXT_ID", "Timepoint_categorical")]
Sample_metadata_mothers_filt <- Sample_metadata_mothers_filt[Sample_metadata_mothers_filt$NG_ID %in% samples_with_DGR_vOTUs,]

# DGR_info_table and external_gene_calls are the same files prepared for infants
dgr_info_table <- read.delim("5_DGR_ANALYSIS/INFANTS/DGR_info_table.tsv")
external_gene_calls <- read.delim("5_DGR_ANALYSIS/external_gene_calls_0based.txt")

# Add to persistent_info the number of Variable regions assessed per vOTU
dgr_counts <- table(dgr_info_table$vOTU_ID)
persistent_info$n_VRs <- dgr_counts[match(persistent_info$vOTU, names(dgr_counts))]
persistent_info$n_VRs[is.na(persistent_info$n_DGRs)] <- 0

# Generate a DGR activity table (one row per VR assessed) to store results
# VR-start coordinates are 0-based (as in dgr_info_table)
dgr_activity_table <- dgr_info_table[dgr_info_table$vOTU_ID %in% persistent_info$vOTU,]
idx <- match(dgr_activity_table$vOTU_ID, persistent_info$vOTU)
dgr_activity_table$NEXT_ID <- persistent_info$persistent_mothers[idx]
dgr_activity_table$NEXT_ID <- gsub("\\s+", "", dgr_activity_table$NEXT_ID)
dgr_activity_table <- dgr_activity_table %>%
  mutate(NEXT_ID = strsplit(NEXT_ID, ",")) %>%
  unnest(NEXT_ID)

# Estimate total number of VRs and DGRs to be assessed
total_VRs_assessed <- nrow(dgr_activity_table) # n=3825
total_DGRs_assessed <- dgr_activity_table %>% # n=3340
  group_by(DGR_code) %>%
  filter(VR_ID %in% unique(VR_ID)[1]) %>%
  ungroup() %>%
  nrow()

# Read prodigal-gv gene predictions (all vOTUs)
prodigal_nt <- readDNAStringSet("5_DGR_ANALYSIS/vOTUs_genes.fna")
prodigal_info <- data.frame(
  gene_id = names(prodigal_nt),
  stringsAsFactors = FALSE
)

prodigal_info$contig <- sub("_[0-9]+$", "", sub(" #.*", "", prodigal_info$gene_id))
prodigal_info$start <- as.numeric(str_extract(prodigal_info$gene_id, "(?<=# )\\d+")) - 1
prodigal_info$end   <- as.numeric(str_extract_all(prodigal_info$gene_id, "(?<=# )\\d+", simplify = TRUE)[,2]) 
prodigal_info$NT_sequence <- as.character(prodigal_nt)

# Add NT sequence to external_gene_calls
external_gene_calls <- external_gene_calls %>%
  left_join(
    prodigal_info %>% 
      select(contig, start, end, NT_sequence),
    by = c("contig" = "contig", "start" = "start", "stop" = "end")
  )

# Write files
write.table(DGR_persistent_vOTUs,file = "5_DGR_ANALYSIS/MOTHERS/persistent_DGR_vOTUs.txt",
            quote = FALSE,row.names = FALSE, col.names = FALSE)
write.table(samples_with_DGR_vOTUs, file = "5_DGR_ANALYSIS/MOTHERS/samples_with_DGR_persistent_vOTUs.txt",
            quote = FALSE,row.names = FALSE, col.names = FALSE)
write.table(mothers_with_DGR_vOTUs, file = "5_DGR_ANALYSIS/MOTHERS/mothers_with_DGR_persistent_vOTUs.txt",
            quote = FALSE,row.names = FALSE, col.names = FALSE)
write.table(persistent_info, file = "5_DGR_ANALYSIS/MOTHERS/vOTU_persistence_mother_per_vOTU.txt",
            sep = "\t", quote = FALSE,row.names = FALSE)
write.table(persistent_info_mother, file = "5_DGR_ANALYSIS/MOTHERS/vOTU_persistence_per_mother.txt",
            sep = "\t", quote = FALSE,row.names = FALSE)
write.table(Sample_metadata_mothers_filt, file = "5_DGR_ANALYSIS/MOTHERS/Sample_metadata_mothers_filt.txt",
            sep = "\t", quote = FALSE,row.names = FALSE)


#***********************************************
# 2. DGR activity analysis: Anvio results
#***********************************************

# For the activity analysis we removed the following 3 cases (problems with activity estimation)
#---LLNEXT003593: LLNEXT_31063
#---LLNEXT011756: MGV_42317
#--LLNEXT202988: LLNEXT_111682

# --- A.1a. Read NT variability results and select relevant columns ---

# Define the folder containing all NT variability files
nt_folder <- "5_DGR_ANALYSIS/MOTHERS/2_VARIABILITY/"
nt_files <- list.files(path = nt_folder, pattern = "_NT_variability.txt$", 
                       full.names = TRUE, recursive = TRUE)

# Read all NT files into a named list
# Important: Note that files have been renamed to add LLNEXT ID as prefix (so each file as a unique name)
NT_results_list <- lapply(nt_files, function(file) {
  df <- read.delim(file, stringsAsFactors = FALSE) %>%
    mutate(across(everything(), as.character))  
  
  # Convert numeric columns back
  numeric_cols <- c("pos_in_contig", "codon_number", "base_pos_in_codon", "gene_length", "entropy")
  df[numeric_cols] <- lapply(df[numeric_cols], as.numeric)
  
  # Keep only relevant columns
  df[, c("pos_in_contig", "sample_id", "corresponding_gene_call",
         "codon_number", "base_pos_in_codon","gene_length", 
         "consensus", "entropy")]
})

# Name the list elements
names(NT_results_list) <- basename(nt_files)

# --- A.1b. Read AA variability results ---
aa_folder <- "5_DGR_ANALYSIS/MOTHERS/2_VARIABILITY/"   
aa_files <- list.files(path = aa_folder, pattern = "_AA_variability.txt$", 
                       full.names = TRUE, recursive = TRUE)

# Read all AA files into a named list
AA_results_list <- lapply(aa_files, function(file) {
  df <- read.delim(file, stringsAsFactors = FALSE) %>%
    mutate(across(everything(), as.character))  
  
  # Convert numeric columns back
  numeric_cols <- c("codon_number", "entropy")
  df[numeric_cols] <- lapply(df[numeric_cols], as.numeric)
  
  # Keep only relevant columns
  df[, c("codon_number", "sample_id", "corresponding_gene_call", "consensus", "entropy", "coverage")]
})

# Name the list elements
names(AA_results_list) <- basename(aa_files)

# --- Process each NT and AA file separately ---
for (i in seq_along(NT_results_list)) {
  
  nt_file_name <- names(NT_results_list)[i]
  NT_results_filt <- NT_results_list[[nt_file_name]]
  
  aa_file_name <- sub("_NT_variability.txt$", "_AA_variability.txt", nt_file_name)
  AA_results_filt <- AA_results_list[[aa_file_name]] # NULL is some cases
  
  cat("Processing NT file:", nt_file_name, 
      if(!is.null(AA_results_filt)) paste("and its AA file:", aa_file_name) else "(no AA file)", 
      "(", i, "/", length(NT_results_list), ")\n")
  
  # Anvio may report variability right outside the coding region (codon order in gene = -1)
  # Remove that information from the table and if the table is empty, skip analysis for this file
  NT_results_filt <- NT_results_filt %>%
    filter(codon_number > 0)
  
  if (nrow(NT_results_filt) == 0) {
    cat("No valid codon positions in file:", nt_file_name, "— skipping analysis.\n")
    next
  }
  
  # --- A.2. If different corresponding_gene_call values, we split into different tables ---
  gene_calls <- unique(NT_results_filt$corresponding_gene_call)
  
  # Create list of tables (each VR is a different table in the list)
  NT_results_by_gene <- lapply(gene_calls, function(gene_id) {
    NT_sub <- NT_results_filt[NT_results_filt$corresponding_gene_call == gene_id, ]
    NT_sub
  })
  names(NT_results_by_gene) <- gene_calls
  
  AA_results_by_gene <- lapply(gene_calls, function(gene_id) {
    AA_sub <- AA_results_filt[AA_results_filt$corresponding_gene_call == gene_id, ]
    AA_sub
  })
  
  names(AA_results_by_gene) <- gene_calls
  
  # Loop over each gene_call in the list
  for (gene_id in names(NT_results_by_gene)) {
    
    NT_results_filt_gene <- NT_results_by_gene[[gene_id]]
    AA_results_filt_gene <- AA_results_by_gene[[gene_id]]
    
    # Define mother ID, samples and gene_call_id
    gene_id_call <- gene_id
    all_samples <- unique(NT_results_filt_gene$sample_id)
    n_samples=length(all_samples)
    mother_ID <- unique(Sample_metadata_mothers[Sample_metadata_mothers$NG_ID %in% all_samples, "NEXT_ID"])
    
    cat("Processing gene_call:", gene_id_call, " in mother:", mother_ID, "\n")
    
    # --- A.3. Add to NT variability results also the non-variable positions ---
    gene_info <- external_gene_calls %>%
      filter(gene_callers_id == gene_id_call) %>%
      select(contig, start, stop, direction, NT_sequence, aa_sequence)
    
    NT_sequence <- strsplit(gene_info$NT_sequence, "")[[1]]
    codon_number <- (seq_along(NT_sequence) - 1) %/% 3 + 1
    
    NT_gene_pos_info <- data.frame(
      NT = rep(NT_sequence, times = length(all_samples)),
      codon_number = rep(codon_number, times = length(all_samples)),
      base_pos_in_codon = rep(((seq_along(NT_sequence) - 1) %% 3) + 1, times = length(all_samples)),
      sample_id = rep(all_samples, each = length(NT_sequence))
    )
    
    NT_gene_pos_info <- NT_gene_pos_info[order(NT_gene_pos_info$codon_number,
                                               NT_gene_pos_info$base_pos_in_codon,
                                               NT_gene_pos_info$sample_id), ]
    row.names(NT_gene_pos_info) <- NULL  
    
    NT_results_full <- merge(
      NT_gene_pos_info,
      NT_results_filt_gene,
      by = c("codon_number", "base_pos_in_codon", "sample_id"),
      all.x = TRUE  
    )
    
    NT_results_full$corresponding_gene_call <- gene_id_call
    NT_results_full$gene_length <- unique(NT_results_full$gene_length[!is.na(NT_results_full$gene_length)])
    NT_results_full$entropy[is.na(NT_results_full$entropy)] <- 0
    
    positions <- NT_results_full$pos_in_contig[!is.na(NT_results_full$pos_in_contig)]
    
    if (all(diff(positions) >= 0)) {
      strand <- "increasing"
      pos_seq <- gene_info$start:gene_info$stop
      pos_seq_trim <- pos_seq[-length(pos_seq)]
    } else if (all(diff(positions) <= 0)) {
      strand <- "decreasing"
      pos_seq <- gene_info$stop:gene_info$start
      pos_seq_trim <- pos_seq[-1] 
    } else {
      stop("Positions are mixed; cannot determine strand.")
    }
    
    NT_results_full$pos_in_contig_full <- rep(pos_seq_trim, each = n_samples)
    
    # --- A.4. Delimit the region corresponding to VR and add TR sequence ---
    vr_row <- dgr_info_table[dgr_info_table$gene_callers_id == gene_id_call, ]
    VR_start <- vr_row$VR_start + 1
    TR_seq   <- vr_row$TR_seq
    VR_length <- nchar(TR_seq)
    VR_end <- VR_start + VR_length - 1
    
    NT_results_full$VR_region <- ifelse(
      NT_results_full$pos_in_contig_full >= VR_start & 
        NT_results_full$pos_in_contig_full <= VR_end,
      "yes", "no"
    )
    
    TR_bases <- strsplit(TR_seq, "")[[1]]
    TR_bases_expanded <- rep(TR_bases, each = n_samples)
    
    NT_results_full$TR_base <- NA
    NT_results_full$TR_base[NT_results_full$VR_region == "yes"] <- TR_bases_expanded
    
    # --- A.5 Longitudinal VR diversification ---
  
    NT_results_VR <- NT_results_full[NT_results_full$VR_region == "yes", ] # Subset VR positions
    
    # Check for variation at each pos_in_contig_full (ignore NAs in consensus)
    VR_diversification <- NT_results_VR %>%
      group_by(pos_in_contig_full) %>%
      summarise(consensus_var = n_distinct(consensus[!is.na(consensus)]) > 1,
                .groups = "drop")
    
    # Number of VR positions showing variation
    num_var_positions <- sum(VR_diversification$consensus_var)
    
    # If variation in any position, flag DGR as showing longitudinal VR diversification 
    long_NT_diversification_flag <- ifelse(num_var_positions > 0, "yes", "no")
    
    # Add columns to dgr_activity_table
    if (!"long_NT_diversification" %in% colnames(dgr_activity_table)) {
      dgr_activity_table$long_NT_diversification <- NA_character_
    }
    if (!"num_var_VR_positions" %in% colnames(dgr_activity_table)) {
      dgr_activity_table$num_var_VR_positions <- NA_integer_
    }
    
    dgr_activity_table$long_NT_diversification[dgr_activity_table$gene_callers_id == gene_id_call &
                                        dgr_activity_table$NEXT_ID == mother_ID] <- long_NT_diversification_flag
    
    dgr_activity_table$num_var_VR_positions[dgr_activity_table$gene_callers_id == gene_id_call &
                                              dgr_activity_table$NEXT_ID == mother_ID] <- num_var_positions

    # --- A.5b. Determine AA-level first activation, early/late changes, and long diversification ---
    
    timepoint_levels <- c("P12", "P28", "B", "M3")
    
    # Ensure columns exist first
    for(col in c("long_AA_diversification", "first_month_AA_change", "pregnancy_AA_change", "post_pregnancy_AA_change")) {
      if(!col %in% colnames(dgr_activity_table)) dgr_activity_table[[col]] <- NA
    }
    
    # Only proceed if AA results exist for this gene
    if (!is.null(AA_results_filt_gene) && nrow(AA_results_filt_gene) > 0) {
      
      # Mark AA positions in VR
      AA_results_filtered <- AA_results_filt_gene %>%
        mutate(VR_region = ifelse(codon_number %in%
                                    NT_results_VR$codon_number[NT_results_VR$VR_region == "yes"],
                                  "yes", "no")) %>%
        filter(VR_region == "yes")  # keep all positions in VR, ignore coverage (checked consensus AA is provided for ALL)
      
      if (nrow(AA_results_filtered) == 0) {
        # No VR positions, then mark all changes as FALSE/NA
        long_div <- "no"
        first_month_AA_change <- NA_character_
        pregnancy_AA_change <- FALSE
        post_pregnancy_AA_change <- FALSE
      } else {
        # Add timepoints
        VR_results_AA <- AA_results_filtered %>%
          left_join(Sample_metadata_mothers %>% select(NG_ID, Timepoint_categorical), 
                    by = c("sample_id" = "NG_ID")) %>%
          mutate(Timepoint_categorical = factor(Timepoint_categorical,
                                                levels = timepoint_levels, ordered = TRUE))
        
        # Build full VR consensus per sample
        sample_vr_aa <- VR_results_AA %>%
          group_by(sample_id, Timepoint_categorical) %>%
          summarise(vr_aa_seq = paste0(consensus, collapse = ""), .groups = "drop") %>%
          arrange(Timepoint_categorical)
        
        # --- Long AA diversification across all samples ---
        long_div <- if(length(unique(sample_vr_aa$vr_aa_seq)) > 1) "yes" else "no"
        
        # --- First activation month ---
        if(nrow(sample_vr_aa) < 2) {
          first_month_AA_change <- NA_character_
        } else {
          ref_seq <- sample_vr_aa$vr_aa_seq[1]  # first sample as reference
          diff_idx <- which(sample_vr_aa$vr_aa_seq != ref_seq)
          first_month_AA_change <- if(length(diff_idx) == 0) NA_character_ else as.character(sample_vr_aa$Timepoint_categorical[diff_idx[1]])
        }
        
        # --- Pregnancy and post-pregnancy changes ---
        pregnancy_tps <- c("P12","P28","B")
        
        pregnancy_seqs <- sample_vr_aa %>% filter(Timepoint_categorical %in% pregnancy_tps) %>% pull(vr_aa_seq)
        pregnancy_AA_change <- if(length(pregnancy_seqs) < 2) NA else any(outer(pregnancy_seqs, pregnancy_seqs, `!=`))
        
        post_preg_seqs <- sample_vr_aa %>% filter(Timepoint_categorical %in% c("B","M3")) %>% pull(vr_aa_seq)
        post_pregnancy_AA_change <- if(length(post_preg_seqs) < 2) NA else any(outer(post_preg_seqs, post_preg_seqs, `!=`))
      }
      
      # Store results in dgr_activity_table
      dgr_activity_table$long_AA_diversification[dgr_activity_table$gene_callers_id == gene_id_call &
                                                   dgr_activity_table$NEXT_ID == mother_ID] <- long_div
      dgr_activity_table$first_month_AA_change[dgr_activity_table$gene_callers_id == gene_id_call &
                                                 dgr_activity_table$NEXT_ID == mother_ID] <- first_month_AA_change
      dgr_activity_table$pregnancy_AA_change[dgr_activity_table$gene_callers_id == gene_id_call &
                                           dgr_activity_table$NEXT_ID == mother_ID] <- pregnancy_AA_change
      dgr_activity_table$post_pregnancy_AA_change[dgr_activity_table$gene_callers_id == gene_id_call &
                                          dgr_activity_table$NEXT_ID == mother_ID] <- post_pregnancy_AA_change
    }
    
    
    # --- A.6 Test VR vs non-VR per-site variability ---
    VR_variability_test_data <- NT_results_full %>%
      group_by(pos_in_contig_full, VR_region) %>%
      summarise(mean_dep = mean(as.numeric(entropy), na.rm = TRUE)) #mean across longitudinal samples
    
    VR_variability_test_data$VR_region <- factor(VR_variability_test_data$VR_region, levels = c("yes", "no"))
    var_test <- wilcox.test(mean_dep ~ VR_region, data = VR_variability_test_data,
                            alternative = "greater", exact = FALSE)
    
    if (!"P_value_VR_var" %in% colnames(dgr_activity_table)) {
      dgr_activity_table$P_value_VR_var <- NA_real_
    }
    dgr_activity_table$P_value_VR_var[dgr_activity_table$gene_callers_id == gene_id_call &
                                        dgr_activity_table$NEXT_ID == mother_ID] <- var_test$p.value
    
    
    # --- A.7 Identify codon (AA) and base-level changes ---
    if (!is.null(AA_results_filt) && nrow(AA_results_filt) > 0) {
      
      if (long_NT_diversification_flag == "yes") {
        # Merge NT variability results with AA results
        NT_results_VR <- NT_results_VR %>%
          mutate(corresponding_gene_call = as.character(corresponding_gene_call))
        
        AA_results <- AA_results_filt_gene %>%
          mutate(corresponding_gene_call = as.character(corresponding_gene_call)) %>%
          select(codon_number, sample_id, corresponding_gene_call, consensus, entropy) 
        colnames(AA_results)[colnames(AA_results) == "consensus"] <- "consensus_AA"
        colnames(AA_results)[colnames(AA_results) == "entropy"] <- "entropy_AA"
        
        NY_AA_results <- NT_results_VR %>%
          left_join(
            AA_results,
            by = c("codon_number", "sample_id", "corresponding_gene_call")
          )
        
        # Select only codons with changes in NT consensus (in any base of the codon) across samples (in the VR region)
        NY_AA_results_changes <- NY_AA_results %>%
          filter(!is.na(consensus)) %>% 
          group_by(codon_number, base_pos_in_codon) %>%
          filter(n_distinct(consensus[!is.na(consensus)]) > 1) %>%
          ungroup()
        
        # (1) Codons with NT consensus changes
        codons_with_NT_changes <- NY_AA_results_changes %>% distinct(codon_number)
        n_codons_NT <- nrow(codons_with_NT_changes)
        
        # (2) Codons with AA changes (restricted to NT-changing codons)
        codons_with_AA_changes <- NY_AA_results %>%
          filter(codon_number %in% codons_with_NT_changes$codon_number) %>%
          group_by(codon_number) %>%
          filter(n_distinct(consensus_AA[!is.na(consensus_AA)]) > 1) %>%
          ungroup() %>%
          distinct(codon_number)
        
        n_codons_AA <- nrow(codons_with_AA_changes)
        
        # Position bias of mutations in codons (fraction of changes at 3rd position)
        codon_position_counts <- NY_AA_results_changes %>%
          select(codon_number, base_pos_in_codon) %>%
          distinct() %>%  
          group_by(base_pos_in_codon) %>%
          summarise(change_count = n(), .groups = "drop")
        
        # Extract counts for each position
        pos1_changes <- codon_position_counts %>% filter(base_pos_in_codon == 1) %>% pull(change_count)
        pos2_changes <- codon_position_counts %>% filter(base_pos_in_codon == 2) %>% pull(change_count)
        pos3_changes <- codon_position_counts %>% filter(base_pos_in_codon == 3) %>% pull(change_count)
        
        # If no changes in a position, set to 0
        if (length(pos1_changes) == 0) pos1_changes <- 0
        if (length(pos2_changes) == 0) pos2_changes <- 0
        if (length(pos3_changes) == 0) pos3_changes <- 0
        
        # Fraction of codons with non-synonymous changes
        fraction_non_syn <- n_codons_AA / n_codons_NT
        if (length(fraction_non_syn) == 0) fraction_non_syn <- 0
        
      } else {
        # If no VR diversification, then set codon metrics to NA
        n_codons_NT <- 0 #if long_NT_diversification_flag == "no", then no codons with NT changes
        n_codons_AA <- 0 #if long_NT_diversification_flag == "no", then no codons with AA changes
        fraction_non_syn <- NA_real_ # NA (not 0 - as it cannot be estimated)
        pos1_changes <- NA_real_ # NA (not 0 - as it cannot be estimated)
        pos2_changes <- NA_real_ # NA (not 0 - as it cannot be estimated)
        pos3_changes <- NA_real_ # NA (not 0 - as it cannot be estimated)
      }
      
      # Store results in dgr_activity_table
      if (!"codons_with_NT_changes" %in% colnames(dgr_activity_table)) {
        dgr_activity_table$codons_with_NT_changes <- NA_integer_
      }
      if (!"codons_with_AA_changes" %in% colnames(dgr_activity_table)) {
        dgr_activity_table$codons_with_AA_changes <- NA_integer_
      }
      if (!"fraction_non_synonymous" %in% colnames(dgr_activity_table)) {
        dgr_activity_table$fraction_non_synonymous <- NA_real_
      }
      if (!"codon_pos1_changes" %in% colnames(dgr_activity_table)) {
        dgr_activity_table$codon_pos1_changes <- NA_integer_
      }
      if (!"codon_pos2_changes" %in% colnames(dgr_activity_table)) {
        dgr_activity_table$codon_pos2_changes <- NA_integer_
      }
      if (!"codon_pos3_changes" %in% colnames(dgr_activity_table)) {
        dgr_activity_table$codon_pos3_changes <- NA_integer_
      }
      
      dgr_activity_table$codons_with_NT_changes[dgr_activity_table$gene_callers_id == gene_id_call &
                                                  dgr_activity_table$NEXT_ID == mother_ID] <- n_codons_NT

      dgr_activity_table$codons_with_AA_changes[dgr_activity_table$gene_callers_id == gene_id_call &
                                                  dgr_activity_table$NEXT_ID == mother_ID] <- n_codons_AA
      
      dgr_activity_table$fraction_non_synonymous[dgr_activity_table$gene_callers_id == gene_id_call &
                                                   dgr_activity_table$NEXT_ID == mother_ID] <- fraction_non_syn
      
      dgr_activity_table$codon_pos1_changes[dgr_activity_table$gene_callers_id == gene_id_call &
                                              dgr_activity_table$NEXT_ID == mother_ID] <- pos1_changes
      dgr_activity_table$codon_pos2_changes[dgr_activity_table$gene_callers_id == gene_id_call &
                                              dgr_activity_table$NEXT_ID == mother_ID] <- pos2_changes
      dgr_activity_table$codon_pos3_changes[dgr_activity_table$gene_callers_id == gene_id_call &
                                              dgr_activity_table$NEXT_ID == mother_ID] <- pos3_changes
    } else {
      cat("Skipping codon/AA-level analysis (no AA file) for gene_call:", gene_id_call, "\n")
    }
  }
}

# Add BH FDR correction
dgr_activity_table$P_value_VR_var_adj <- p.adjust(dgr_activity_table$P_value_VR_var, method = "BH")

# --- A.9. Add missing results to VRs with no variability ---

# Add 2 columns indicating in any within-sample variability was detected at NT or AA level (anvio output present)
dgr_activity_table$Within_sample_NT_var_target <- ifelse(is.na(dgr_activity_table$long_NT_diversification),
                                                         "no","yes")

dgr_activity_table$Within_sample_AA_var_target <- ifelse(is.na(dgr_activity_table$long_AA_diversification),
                                                         "no","yes")

# VRs with no within-sample variability detected (in the target gene) do not generate an output file in anvio
# This is due to the lack of variability (and no other reason) and should be counted as non-variable VRs
# Set to 0 (instead of NA) longitudinal variability results for those VRs
no_var_rows_NT <- which(is.na(dgr_activity_table$long_NT_diversification))
no_var_rows_AA <- which(is.na(dgr_activity_table$long_AA_diversification))
dgr_activity_table$long_NT_diversification[no_var_rows_NT] <- "no"
dgr_activity_table$num_var_VR_positions[no_var_rows_NT] <- 0
dgr_activity_table$codons_with_NT_changes[no_var_rows_NT] <- 0
dgr_activity_table$codons_with_AA_changes[no_var_rows_NT] <- 0
dgr_activity_table$long_AA_diversification[no_var_rows_AA] <- "no"
dgr_activity_table$pregnancy_AA_change[no_var_rows_AA] <- FALSE #NAs left in this column: No enough pregnancy timepoints to check 
dgr_activity_table$post_pregnancy_AA_change[no_var_rows_AA] <- FALSE #NAs left in this column: No enough post-pregnancy timepoints to check (B and M3)

# Reorder colums
dgr_activity_table <- dgr_activity_table[, c("vOTU_ID", "DGR_code", "DGR_target", "VR_start", "VR_end","TR_seq",
                                             "Adjusted_VR_sequence", "Adjusted_VR_amino_acid_seq","Adjusted_bias_ATCG",
                                             "VR_position_within_Target", "VR_ID","gene_callers_id", "NEXT_ID",
                                             "Within_sample_NT_var_target", "Within_sample_AA_var_target", "P_value_VR_var",
                                             "P_value_VR_var_adj","long_NT_diversification","num_var_VR_positions",
                                             "long_AA_diversification", "first_month_AA_change", "pregnancy_AA_change",
                                             "post_pregnancy_AA_change", "codons_with_NT_changes","codons_with_AA_changes",
                                             "fraction_non_synonymous", "codon_pos1_changes","codon_pos2_changes","codon_pos3_changes")]

# Change FALSE/TRUE columns to yes no (early_AA_change, late_AA_change, early_to_late_AA_change)
dgr_activity_table$pregnancy_AA_change       <- ifelse(dgr_activity_table$pregnancy_AA_change, "yes", "no")
dgr_activity_table$post_pregnancy_AA_change  <- ifelse(dgr_activity_table$post_pregnancy_AA_change, "yes", "no")

# --- A.10. Add annotation of target proteins ---
dgr_activity_table <- dgr_activity_table %>%
  left_join(
    Annotation_DGR_targets %>%
      select(Protein_ID, 
             DGR_target_ann_description = DESCRIPTION, 
             DGR_target_ann_category = CATEGORY),
    by = c("DGR_target" = "Protein_ID")
  )

# Summarize annotations for each of the unique targets of active DGRs
active_DGR_target_summary <- dgr_activity_table %>%
  filter(long_NT_diversification == "yes") %>%    
  select(DGR_target, DGR_target_ann_category) %>%  
  distinct() %>%                                    
  arrange(DGR_target)                      

#****************
# Write results
#****************
write.table(dgr_activity_table,"5_DGR_ANALYSIS/MOTHERS/DGR_activity_table_MOTHERS.txt", sep = "\t", row.names = F, quote = FALSE) 

