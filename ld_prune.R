#!/usr/bin/env Rscript
# ld_prune.R
#
# This script reads a compressed PLINK LD matrix (.txt.gz) without a header
# along with a corresponding .bim file to obtain variant IDs. A lead SNP is provided
# via optparse and is forced to be selected as an independent variant regardless of its LD.
#
# The algorithm works as follows:
# 1. Read the .bim file to get the list of variant IDs.
# 2. Read the compressed LD matrix and assign the variant IDs as row and column names.
# 3. If a lead SNP is provided:
#    - Force the lead SNP into the independent set.
#    - Exclude any variant with high LD (|r²| > threshold) with the lead SNP.
#    - From the remaining variants, perform clumping (in original .bim order) to select additional independent variants.
# 4. If no lead SNP is provided, perform standard clumping on all variants.
#
# Usage:
#   Rscript ld_prune.R --ld ld_matrix.txt.gz --bim variants.bim --lead <lead_SNP_ID> --threshold 0.1 --output independent_variants.txt

suppressWarnings(suppressMessages(library(optparse)))

# Define command-line options.
option_list <- list(
  make_option(c("-l", "--ld"), type = "character", default = NULL,
              help = "Input compressed PLINK LD matrix file (.txt.gz) without header", metavar = "character"),
  make_option(c("-b", "--bim"), type = "character", default = NULL,
              help = "Input .bim file containing variant information", metavar = "character"),
  make_option(c("--lead"), type = "character", default = NULL,
              help = "Lead SNP ID that must always be selected", metavar = "character"),
  make_option(c("-t", "--threshold"), type = "double", default = 0.1,
              help = "LD threshold for pruning (default = 0.1)", metavar = "double"),
  make_option(c("-o", "--output"), type = "character", default = NULL,
              help = "Output file to save independent variant IDs", metavar = "character")
)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# Check if required inputs are provided.
if (is.null(opt$ld) || is.null(opt$bim) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("Please specify --ld, --bim, and --output.", call. = FALSE)
}

# Read the .bim file.
# Expected format: chromosome, variant_ID, genetic distance, base-pair position, allele1, allele2 (no header).
bim_data <- read.table(opt$bim, header = FALSE, stringsAsFactors = FALSE)
if (ncol(bim_data) < 2) {
  stop("The .bim file must have at least two columns.", call. = FALSE)
}
variant_ids <- bim_data[, 2]

# Read the compressed LD matrix file (no header).
ld_matrix <- as.matrix(read.table(gzfile(opt$ld), header = FALSE, stringsAsFactors = FALSE))
num_variants <- length(variant_ids)
if (nrow(ld_matrix) != num_variants || ncol(ld_matrix) != num_variants) {
  stop("Mismatch between number of variants in .bim file and dimensions of the LD matrix.", call. = FALSE)
}
rownames(ld_matrix) <- variant_ids
colnames(ld_matrix) <- variant_ids

# Initialize independent variants list.
independent <- c()

if (!is.null(opt$lead)) {
  if (!(opt$lead %in% variant_ids)) {
    stop(sprintf("The specified lead SNP '%s' is not found in the .bim file.", opt$lead), call. = FALSE)
  }
  
  # Force the lead SNP into the independent set.
  independent <- c(opt$lead)
  
  # Create a logical vector indicating candidate variants (all variants except the lead SNP).
  candidate_flags <- rep(TRUE, num_variants)
  names(candidate_flags) <- variant_ids
  candidate_flags[opt$lead] <- FALSE
  
  # Exclude any variant in high LD with the lead SNP.
  for (var in variant_ids[candidate_flags]) {
    if (abs(ld_matrix[opt$lead, var]) > opt$threshold) {
      candidate_flags[var] <- FALSE
    }
  }
  
  # Get the remaining variants from the original .bim order.
  remaining_variants <- variant_ids[candidate_flags]
  
  # Perform standard clumping on the remaining variants.
  if (length(remaining_variants) > 0) {
    additional_independents <- c()
    remain_flags <- rep(TRUE, length(remaining_variants))
    names(remain_flags) <- remaining_variants
    for (i in seq_along(remaining_variants)) {
      if (!remain_flags[i]) next
      current_var <- remaining_variants[i]
      additional_independents <- c(additional_independents, current_var)
      if (i < length(remaining_variants)) {
        for (j in (i + 1):length(remaining_variants)) {
          if (remain_flags[j]) {
            ld_val <- abs(ld_matrix[current_var, remaining_variants[j]])
            if (ld_val > opt$threshold) {
              remain_flags[j] <- FALSE
            }
          }
        }
      }
    }
    independent <- c(independent, additional_independents)
  }
  
} else {
  # If no lead SNP is specified, perform standard clumping on all variants.
  candidate_flags <- rep(TRUE, num_variants)
  names(candidate_flags) <- variant_ids
  for (i in seq_along(variant_ids)) {
    if (!candidate_flags[i]) next
    current_var <- variant_ids[i]
    independent <- c(independent, current_var)
    if (i < num_variants) {
      for (j in (i + 1):num_variants) {
        if (candidate_flags[j]) {
          ld_val <- abs(ld_matrix[current_var, variant_ids[j]])
          if (ld_val > opt$threshold) {
            candidate_flags[j] <- FALSE
          }
        }
      }
    }
  }
}

cat(sprintf("Selected %d independent variants out of %d total variants.\n", length(independent), num_variants))
writeLines(independent, con = opt$output)
cat(sprintf("Independent variant IDs saved to '%s'.\n", opt$output))
