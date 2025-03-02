#' Export Simulated Zero-Inflated Virome Data to Common File Formats
#'
#' Converts simulated virome data with variable zero-inflation to common file formats
#' used in microbiome and virome analysis, facilitating integration with external tools
#' and workflows.
#'
#' @param sim_data A list object returned by \code{simulate_zero_inflated_virome}
#' @param output_dir Directory to save output files (default: current working directory)
#' @param base_filename Base name for output files (default: "virome_simulation")
#' @param formats Vector of output formats to generate (default: c("csv", "biom", "phyloseq", "qiime2"))
#' @param include_metadata Whether to include the metadata in the exports (default: TRUE)
#' @param include_zi_parameters Whether to include zero-inflation parameters in the exports (default: TRUE)
#' @param normalize Whether to normalize counts before export (default: FALSE)
#' @param normalization_method Method for normalization: "relative", "TMM", "CSS", "CPM" (default: "relative")
#' @param taxonomy_levels Number of taxonomic levels to simulate (default: 7)
#'
#' @return A list containing:
#'   \item{exported_files}{Character vector of paths to exported files}
#'   \item{formats}{Summary of formats and their status}
#'   \item{summary}{Statistics about the exported data}
#'
#' @details
#' This function facilitates the use of simulated data with variable zero-inflation in
#' external analysis tools by converting it to standard formats used in microbiome research.
#' It supports the following output formats:
#'
#' \itemize{
#'   \item CSV: Simple count matrix and metadata as CSV files
#'   \item TSV: Tab-separated format compatible with various tools
#'   \item BIOM: Biological Observation Matrix format (requires \code{biomformat} package)
#'   \item phyloseq: R object for the phyloseq package (requires \code{phyloseq} package)
#'   \item QIIME2: Files compatible with QIIME2 import (requires \code{qiime2R} package)
#'   \item Mothur: Files compatible with Mothur (shared and metadata files)
#' }
#'
#' The function generates simulated taxonomy information to make the exports compatible
#' with standard tools. It also provides options to normalize the count data and to
#' include zero-inflation parameters as part of the metadata.
#'
#' @examples
#' \dontrun{
#' # Generate simulated data with variable zero-inflation
#' sim_data <- simulate_zero_inflated_virome(
#'   n_samples = 30,
#'   n_viruses = 100,
#'   variable_zi_rates = TRUE,
#'   zi_alpha = 2,
#'   zi_beta = 5
#' )
#' 
#' # Export to CSV and BIOM formats
#' export_files <- export_zinb_simulation(
#'   sim_data = sim_data,
#'   formats = c("csv", "biom"),
#'   output_dir = "output",
#'   base_filename = "variable_zi_simulation"
#' )
#' 
#' # Export normalized data for visualization
#' export_norm <- export_zinb_simulation(
#'   sim_data = sim_data,
#'   formats = c("csv"),
#'   normalize = TRUE,
#'   normalization_method = "CPM",
#'   base_filename = "normalized_sim"
#' )
#' }
#'
#' @export
export_zinb_simulation <- function(sim_data,
                                  output_dir = getwd(),
                                  base_filename = "virome_simulation",
                                  formats = c("csv", "biom", "phyloseq", "qiime2"),
                                  include_metadata = TRUE,
                                  include_zi_parameters = TRUE,
                                  normalize = FALSE,
                                  normalization_method = "relative",
                                  taxonomy_levels = 7) {
  
  # Check if the input is a valid simulation result
  if (!is.list(sim_data) || !all(c("counts", "metadata") %in% names(sim_data))) {
    stop("Invalid sim_data. Must be a list containing 'counts' and 'metadata' as returned by simulate_zero_inflated_virome().")
  }
  
  # Create output directory if it doesn't exist
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }
  
  # Initialize results
  exported_files <- character()
  format_status <- list()
  
  # Extract data from simulation
  counts <- sim_data$counts
  metadata <- sim_data$metadata
  
  # Check if we're dealing with variable ZI data
  has_variable_zi <- !is.null(sim_data$parameters$variable_zi_rates) && 
    sim_data$parameters$variable_zi_rates
  
  # Add ZI parameters to metadata if requested
  if (include_zi_parameters && has_variable_zi) {
    metadata$mean_zi_rate <- sim_data$parameters$structural_zeros
    if (!is.null(sim_data$parameters$zi_alpha)) {
      metadata$zi_alpha <- sim_data$parameters$zi_alpha
      metadata$zi_beta <- sim_data$parameters$zi_beta
    }
  }
  
  # Generate simulated taxonomy information
  n_viruses <- nrow(counts)
  tax_levels <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_levels <- tax_levels[1:min(taxonomy_levels, length(tax_levels))]
  
  # Create taxonomy data frame
  taxonomy <- data.frame(
    OTU_ID = rownames(counts),
    stringsAsFactors = FALSE
  )
  
  # Fill in taxonomy levels with simulated names
  set.seed(42)  # For reproducible taxonomy
  
  # Number of unique taxa at each level (decreases at higher resolution)
  n_taxa_per_level <- round(n_viruses / c(1, 2, 4, 8, 16, 32, 64)[1:length(tax_levels)])
  n_taxa_per_level <- pmin(n_taxa_per_level, n_viruses)
  
  # Generate taxonomy
  for (i in seq_along(tax_levels)) {
    level <- tax_levels[i]
    n_taxa <- n_taxa_per_level[i]
    
    # For viruses, use viral taxonomy
    prefixes <- switch(level,
                      Kingdom = c("Virus"),
                      Phylum = c("DNA_virus", "RNA_virus", "Retrovirus"),
                      Class = c("ss", "ds", "circular", "linear"),
                      Order = c("Caudoviral", "Herpesvi", "Circo", "Parvo", "Corona"),
                      Family = c("viridae", "virus", "phage"),
                      Genus = letters[1:10],
                      Species = letters[1:10])
    
    # Generate names for this level
    taxa_names <- paste0(
      sample(prefixes, n_taxa, replace = TRUE),
      "_",
      sample(1:100, n_taxa, replace = FALSE)
    )
    
    # Assign to taxonomy dataframe
    taxonomy[[level]] <- sample(taxa_names, n_viruses, replace = TRUE)
  }
  
  # Normalize counts if requested
  if (normalize) {
    if (normalization_method == "relative") {
      # Simple relative abundance (proportions)
      norm_counts <- sweep(counts, 2, colSums(counts), FUN = "/")
      norm_counts <- norm_counts * 100  # Convert to percentages
    } else if (normalization_method == "CPM") {
      # Counts per million
      norm_counts <- sweep(counts, 2, colSums(counts), FUN = "/") * 1e6
    } else if (normalization_method == "TMM") {
      # Trimmed Mean of M-values
      if (!requireNamespace("edgeR", quietly = TRUE)) {
        warning("Package 'edgeR' is required for TMM normalization. Falling back to relative abundance.")
        norm_counts <- sweep(counts, 2, colSums(counts), FUN = "/") * 100
      } else {
        dge <- edgeR::DGEList(counts = counts)
        dge <- edgeR::calcNormFactors(dge, method = "TMM")
        norm_counts <- edgeR::cpm(dge)
      }
    } else if (normalization_method == "CSS") {
      # Cumulative Sum Scaling
      if (!requireNamespace("metagenomeSeq", quietly = TRUE)) {
        warning("Package 'metagenomeSeq' is required for CSS normalization. Falling back to relative abundance.")
        norm_counts <- sweep(counts, 2, colSums(counts), FUN = "/") * 100
      } else {
        mgseq <- metagenomeSeq::newMRexperiment(counts)
        mgseq <- metagenomeSeq::cumNorm(mgseq)
        norm_counts <- metagenomeSeq::MRcounts(mgseq, normalized = TRUE)
      }
    } else {
      stop("Invalid normalization_method. Choose from: 'relative', 'TMM', 'CSS', 'CPM'")
    }
    
    # Replace counts with normalized values for export
    counts_for_export <- norm_counts
  } else {
    counts_for_export <- counts
  }
  
  # Export to requested formats
  
  # 1. CSV format
  if ("csv" %in% formats) {
    # Write count matrix
    count_file <- file.path(output_dir, paste0(base_filename, "_counts.csv"))
    utils::write.csv(counts_for_export, file = count_file, row.names = TRUE)
    exported_files <- c(exported_files, count_file)
    
    # Write metadata
    if (include_metadata) {
      meta_file <- file.path(output_dir, paste0(base_filename, "_metadata.csv"))
      utils::write.csv(metadata, file = meta_file, row.names = FALSE)
      exported_files <- c(exported_files, meta_file)
    }
    
    # Write taxonomy
    tax_file <- file.path(output_dir, paste0(base_filename, "_taxonomy.csv"))
    utils::write.csv(taxonomy, file = tax_file, row.names = FALSE)
    exported_files <- c(exported_files, tax_file)
    
    format_status$csv <- "Success"
  }
  
  # 2. TSV format
  if ("tsv" %in% formats) {
    # Write count matrix
    count_file <- file.path(output_dir, paste0(base_filename, "_counts.tsv"))
    utils::write.table(counts_for_export, file = count_file, 
                      sep = "\t", row.names = TRUE, col.names = TRUE, quote = FALSE)
    exported_files <- c(exported_files, count_file)
    
    # Write metadata
    if (include_metadata) {
      meta_file <- file.path(output_dir, paste0(base_filename, "_metadata.tsv"))
      utils::write.table(metadata, file = meta_file, 
                        sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      exported_files <- c(exported_files, meta_file)
    }
    
    # Write taxonomy
    tax_file <- file.path(output_dir, paste0(base_filename, "_taxonomy.tsv"))
    utils::write.table(taxonomy, file = tax_file, 
                      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    exported_files <- c(exported_files, tax_file)
    
    format_status$tsv <- "Success"
  }
  
  # 3. BIOM format
  if ("biom" %in% formats) {
    if (!requireNamespace("biomformat", quietly = TRUE)) {
      warning("Package 'biomformat' is required for BIOM format export. Skipping BIOM export.")
      format_status$biom <- "Failed - biomformat package not available"
    } else {
      tryCatch({
        # Prepare taxonomy in the format expected by biomformat
        # Convert taxonomy data frame to list of character vectors
        tax_list <- lapply(1:nrow(taxonomy), function(i) {
          as.character(taxonomy[i, -1])  # Exclude OTU_ID
        })
        names(tax_list) <- taxonomy$OTU_ID
        
        # Create sample metadata
        sample_md <- metadata
        rownames(sample_md) <- metadata$sample_id
        
        # Create BIOM object
        biom_obj <- biomformat::make_biom(
          data = as.matrix(counts_for_export),
          observation_metadata = tax_list,
          sample_metadata = as.data.frame(sample_md)
        )
        
        # Write BIOM file
        biom_file <- file.path(output_dir, paste0(base_filename, ".biom"))
        biomformat::write_biom(biom_obj, biom_file)
        exported_files <- c(exported_files, biom_file)
        
        format_status$biom <- "Success"
      }, error = function(e) {
        warning("Error creating BIOM file: ", e$message)
        format_status$biom <- paste("Failed -", e$message)
      })
    }
  }
  
  # 4. phyloseq format
  if ("phyloseq" %in% formats) {
    if (!requireNamespace("phyloseq", quietly = TRUE)) {
      warning("Package 'phyloseq' is required for phyloseq format export. Skipping phyloseq export.")
      format_status$phyloseq <- "Failed - phyloseq package not available"
    } else {
      tryCatch({
        # Create OTU table
        otu_table <- phyloseq::otu_table(as.matrix(counts_for_export), taxa_are_rows = TRUE)
        
        # Create sample data
        sample_data <- phyloseq::sample_data(metadata)
        rownames(sample_data) <- metadata$sample_id
        
        # Create taxonomy table
        tax_matrix <- as.matrix(taxonomy[, -1, drop = FALSE])
        rownames(tax_matrix) <- taxonomy$OTU_ID
        colnames(tax_matrix) <- tax_levels
        tax_table <- phyloseq::tax_table(tax_matrix)
        
        # Create phyloseq object
        physeq <- phyloseq::phyloseq(otu_table, tax_table, sample_data)
        
        # Save as RDS file
        physeq_file <- file.path(output_dir, paste0(base_filename, "_phyloseq.rds"))
        saveRDS(physeq, file = physeq_file)
        exported_files <- c(exported_files, physeq_file)
        
        format_status$phyloseq <- "Success"
      }, error = function(e) {
        warning("Error creating phyloseq object: ", e$message)
        format_status$phyloseq <- paste("Failed -", e$message)
      })
    }
  }
  
  # 5. QIIME2 format
  if ("qiime2" %in% formats) {
    # QIIME2 format requires specific file structure
    # We'll create the necessary files for QIIME2 import
    
    # Feature table
    qiime_feature_file <- file.path(output_dir, paste0(base_filename, "_feature-table.tsv"))
    # QIIME2 expects the first column to be the feature IDs, so we transpose
    qiime_counts <- t(counts_for_export)
    utils::write.table(
      cbind(feature_id = colnames(qiime_counts), qiime_counts),
      file = qiime_feature_file,
      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE
    )
    exported_files <- c(exported_files, qiime_feature_file)
    
    # Taxonomy
    qiime_tax_file <- file.path(output_dir, paste0(base_filename, "_taxonomy.tsv"))
    # Format taxonomy string like "k__Virus; p__DNA_virus; c__ss; ..."
    tax_strings <- apply(taxonomy[, -1, drop = FALSE], 1, function(row) {
      paste0(
        sapply(seq_along(tax_levels), function(i) {
          paste0(substr(tax_levels[i], 1, 1), "__", row[i])
        }),
        collapse = "; "
      )
    })
    
    qiime_tax <- data.frame(
      Feature.ID = taxonomy$OTU_ID,
      Taxon = tax_strings,
      Confidence = round(runif(nrow(taxonomy), 0.8, 1.0), 2)
    )
    
    utils::write.table(qiime_tax, file = qiime_tax_file,
                      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    exported_files <- c(exported_files, qiime_tax_file)
    
    # Metadata
    if (include_metadata) {
      qiime_meta_file <- file.path(output_dir, paste0(base_filename, "_metadata.tsv"))
      # QIIME2 metadata requires sample IDs in the first column, named #SampleID
      qiime_meta <- metadata
      colnames(qiime_meta)[colnames(qiime_meta) == "sample_id"] <- "#SampleID"
      
      utils::write.table(qiime_meta, file = qiime_meta_file,
                        sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
      exported_files <- c(exported_files, qiime_meta_file)
    }
    
    # Create bash script with import commands
    qiime_script_file <- file.path(output_dir, paste0(base_filename, "_qiime2_import.sh"))
    qiime_script <- paste0(
      "#!/bin/bash\n\n",
      "# QIIME2 import commands for simulated virome data\n\n",
      
      "# Import feature table\n",
      "qiime tools import \\\n",
      "  --input-path ", basename(qiime_feature_file), " \\\n",
      "  --output-path ", basename(base_filename), "_table.qza \\\n",
      "  --type 'FeatureTable[Frequency]' \\\n",
      "  --input-format TSVTaxonomyFormat\n\n",
      
      "# Import taxonomy\n",
      "qiime tools import \\\n",
      "  --input-path ", basename(qiime_tax_file), " \\\n",
      "  --output-path ", basename(base_filename), "_taxonomy.qza \\\n",
      "  --type 'FeatureData[Taxonomy]' \\\n",
      "  --input-format TSVTaxonomyFormat\n\n",
      
      "# Import metadata\n",
      "# No import needed, use with --m-metadata-file parameter\n"
    )
    
    cat(qiime_script, file = qiime_script_file)
    exported_files <- c(exported_files, qiime_script_file)
    
    format_status$qiime2 <- "Success"
  }
  
  # 6. Mothur format
  if ("mothur" %in% formats) {
    # Mothur shared file
    mothur_shared_file <- file.path(output_dir, paste0(base_filename, ".mothur.shared"))
    
    # Format for mothur shared file: label group numOtus Otu001 Otu002 ...
    mothur_header <- c("label", "Group", "numOtus", paste0("Otu", sprintf("%03d", 1:ncol(counts_for_export))))
    
    mothur_data <- data.frame(
      label = rep("1", nrow(counts_for_export)),
      Group = rownames(counts_for_export),
      numOtus = rep(ncol(counts_for_export), nrow(counts_for_export)),
      stringsAsFactors = FALSE
    )
    
    # Add count data
    mothur_data <- cbind(mothur_data, t(counts_for_export))
    colnames(mothur_data) <- mothur_header
    
    # Write mothur shared file
    utils::write.table(mothur_data, file = mothur_shared_file,
                      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    exported_files <- c(exported_files, mothur_shared_file)
    
    # Mothur taxonomy
    mothur_tax_file <- file.path(output_dir, paste0(base_filename, ".mothur.taxonomy"))
    
    # Format taxonomy for mothur
    tax_strings <- apply(taxonomy[, -1, drop = FALSE], 1, function(row) {
      paste0(
        sapply(seq_along(tax_levels), function(i) {
          paste0(row[i], "(", round(runif(1, 0, 100)), ")")
        }),
        collapse = ";"
      )
    })
    
    mothur_tax <- data.frame(
      OTU = paste0("Otu", sprintf("%03d", 1:nrow(taxonomy))),
      Taxonomy = paste0(tax_strings, ";"),
      stringsAsFactors = FALSE
    )
    
    utils::write.table(mothur_tax, file = mothur_tax_file,
                      sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE)
    exported_files <- c(exported_files, mothur_tax_file)
    
    format_status$mothur <- "Success"
  }
  
  # Create a summary of the export
  export_summary <- list(
    n_samples = ncol(counts),
    n_taxa = nrow(counts),
    formats_exported = names(format_status)[format_status == "Success"],
    total_files = length(exported_files),
    normalized = normalize,
    normalization_method = if (normalize) normalization_method else "None",
    zero_inflation_type = if (has_variable_zi) "Variable (Beta distribution)" else "Fixed"
  )
  
  # Return a list with paths to exported files and status information
  return(list(
    exported_files = exported_files,
    formats = format_status,
    summary = export_summary
  ))
}