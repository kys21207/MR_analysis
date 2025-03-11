#!/usr/bin/env Rscript
# ld_prune.R
#
# This script reads in a PLINK LD matrix from a compressed file (.txt.gz) that does not contain header information.
# It also reads a corresponding .bim file to obtain the variant IDs.
# The script then performs LD pruning using a specified r² threshold and outputs a list of independent variant IDs.
#
# Usage:
#   Rscript ld_prune.R --ld ld_matrix.txt.gz --bim variants.bim --threshold 0.1 --output independent_variants.txt

suppressWarnings(suppressMessages(library(optparse)))

# Define command-line options using optparse.
option_list <- list(
  make_option(c("-l", "--ld"), type="character", default=NULL,
              help="Input compressed PLINK LD matrix file (.txt.gz) with no header", metavar="character"),
  make_option(c("-b", "--bim"), type="character", default=NULL,
              help="Input .bim file containing variant information; will be used to supply row/column names", metavar="character"),
  make_option(c("-t", "--threshold"), type="double", default=0.1,
              help="LD threshold for pruning (default = 0.1)", metavar="double"),
  make_option(c("-o", "--output"), type="character", default=NULL,
              help="Output file to save independent variant IDs", metavar="character")
)

opt_parser <- OptionParser(option_list=option_list)
opt <- parse_args(opt_parser)

# Check if required files are provided.
if (is.null(opt$ld) || is.null(opt$bim) || is.null(opt$output)) {
  print_help(opt_parser)
  stop("LD matrix file, bim file, and output file must be specified.", call.=FALSE)
}

# Read in the .bim file.
# The .bim file is assumed to be whitespace-delimited with no header.
# Typical columns: chromosome, variant ID, genetic distance, base-pair position, allele1, allele2.
bim_data <- read.table(opt$bim, header=FALSE, stringsAsFactors=FALSE)
if (ncol(bim_data) < 2) {
  stop("The .bim file does not have at least two columns. Please check the file format.")
}
variant_ids <- bim_data[, 2]

# Read in the compressed LD matrix file.
# Since there is no header, set header = FALSE.
ld_matrix <- as.matrix(read.table(gzfile(opt$ld), header=FALSE, stringsAsFactors=FALSE))
num_variants <- length(variant_ids)
if (nrow(ld_matrix) != num_variants || ncol(ld_matrix) != num_variants) {
  stop("Number of variants in the bim file does not match the dimensions of the LD matrix.")
}

# Assign row and column names to the LD matrix using the variant IDs from the bim file.
rownames(ld_matrix) <- variant_ids
colnames(ld_matrix) <- variant_ids

# Function to perform LD pruning based on a specified r² threshold.
ld_pruning <- function(ld_matrix, threshold=0.1) {
  independent_variants <- c()
  variants <- rownames(ld_matrix)
  
  for (variant in variants) {
    keep <- TRUE
    # Check pairwise LD with already selected independent variants.
    for (sel in independent_variants) {
      ld_value <- abs(ld_matrix[variant, sel])
      if (ld_value > threshold) {
        keep <- FALSE
        break
      }
    }
    if (keep) {
      independent_variants <- c(independent_variants, variant)
    }
  }
  
  return(independent_variants)
}

# Perform LD pruning.
indep_variants <- ld_pruning(ld_matrix, threshold=opt$threshold)
cat(sprintf("Found %d independent variants out of %d total variants.\n", length(indep_variants), num_variants))

# Write the independent variant IDs to the output file.
writeLines(indep_variants, con = opt$output)
cat(sprintf("Independent variant IDs saved to '%s'.\n", opt$output))
