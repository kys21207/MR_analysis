# Install TwoSampleMR
install.packages("TwoSampleMR", repos = c("https://mrcieu.r-universe.dev", "https://cloud.r-project.org"))

# Install MRPRESSO from GitHub
devtools::install_github("rondolab/MR-PRESSO")

library(TwoSampleMR)
library(MRPRESSO)
library(data.table)

# Set file paths
eqtl_file <- "/mnt/project/publically_available_supporting_files/rosmap_brain/celltype-eqtl-sumstats.Mic.tsv.gz"  # Update with actual file path
gwas_file <- "/mnt/project/approved_results/pd_w_age_bsize_400_firth_parkinsons_all_comers_v3.regenie.tsv.gz"    # Update with actual file path
ref_file <- "/mnt/project/publically_available_supporting_files/pops_data/demo_file_rsID4regenie.gwaslab.gz"
geneID <- "ENSG00000157601"          # Update with the gene you want to test
gene_symbol <- "MX1"                  # Add gene symbol if available

# Step 1: Read and Filter eQTL Data (Selecting Specific Gene)
eqtl_data <- fread(cmd = paste("zcat", eqtl_file), sep = "\t")

# Select only eQTLs for the gene of interest in microglia
eqtl_filtered <- eqtl_data[gene_id == geneID, 
                           .(SNP = snps, 
                             beta.exposure = beta, 
                             se.exposure = se, 
                             pval.exposure = pvalue, 
                             effect_allele.exposure = ALT, 
                             other_allele.exposure = REF, 
                             eaf.exposure = ALT_AF
                             )]

# Add required columns for TwoSampleMR
eqtl_filtered[, `:=`(id.exposure = gene_symbol,  # Use gene symbol or ID
                     exposure = gene_symbol)]   # Name the exposure column

# Step 2: Read and Format GWAS Data
gwas_data <- fread(cmd = paste("zcat", gwas_file), sep = " ")


ref <- fread(cmd = paste("gunzip -c", shQuote(ref_file)), sep = "\t", header = TRUE)

gwas_data <- merge(gwas_data, ref, by.x ="ID",by.y="SNPID")

# Convert LOG10P to p-value
gwas_data[, pval.outcome := 10^(-LOG10P)]

# Perform the merge operation
merged_data <- merge(
  gwas_data, 
  eqtl_filtered[, .(SNP)], 
  by.x = "rsID", 
  by.y = "SNP"
)

# Format the merged data
gwas_formatted <- merged_data[ , .(
  SNP = rsID, 
  beta.outcome = BETA, 
  se.outcome = SE, 
  pval.outcome = pval.outcome, 
  effect_allele.outcome = ALLELE1, 
  other_allele.outcome = ALLELE0
)]

# Add required identifier columns for TwoSampleMR
gwas_formatted[, `:=`(id.outcome = "GWAS_trait", outcome = "GWAS_trait")]

# Step 3: Harmonize the Data
harmonized_data <- harmonise_data(eqtl_filtered, gwas_formatted)

# Step 4: Perform MR using IVW (Primary Method)
mr_results <- mr(harmonized_data, method_list = c("mr_ivw"))

# Step 4: Sensitivity Analyses
mr_egger_results <- mr(harmonized_data, method_list = c("mr_egger_regression"))
mr_weighted_median <- mr(harmonized_data, method_list = c("mr_weighted_median"))

# Step 6: MR-PRESSO (Detect and Correct for Pleiotropy)
mr_presso_result <- tryCatch(
  {
    result <- mr_presso(
      BetaOutcome = beta_outcome, 
      BetaExposure = beta_exposure, 
      SdOutcome = se_outcome, 
      SdExposure = se_exposure, 
      data = harmonized_data, 
      NbDistribution = 1000, 
      SignifThreshold = 0.05
    )
    # Check if outlier detection output is present.
    if (is.null(result$`Distortion Test`)) {
      cat("MR-PRESSO completed without detecting outliers.\n")
    }
    result
  },
  error = function(e) {
    cat("MR-PRESSO error caught:\n")
    cat(e$message, "\n")
    NULL
  }
)

if (is.null(mr_presso_result)) {
  cat("MR-PRESSO did not complete successfully. This error may occur if no outliers are detected or due to internal indexing issues.\n")
  cat("Proceeding with other MR methods (IVW, MR-Egger, weighted median) for interpretation.\n")
} else {
  cat("MR-PRESSO results:\n")
  print(mr_presso_result)
}


# Step 7: Display Results
print(mr_results)          # Main MR (IVW) results
print(mr_egger_results)    # Egger regression results
print(mr_weighted_median)  # Weighted median MR results
print(mr_presso_result)    # MR-PRESSO pleiotropy results

# Step 8: Generate MR Plots
mr_forest_plot(harmonized_data)      
mr_scatter_plot(mr_results, harmonized_data) 

