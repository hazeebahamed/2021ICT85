
# Rice Genetic Diversity Analysis
# Step 5: SNP Quality Control (Revised)

# 1. Load processed datasets
genotype_matrix <- readRDS(
  "data/genotype_matrix.rds"
)

rice_metadata <- readRDS(
  "data/rice_metadata.rds"
)

# 2. Define QC thresholds
sample_threshold <- 0.15
snp_threshold <- 0.10

# 3. Identify missing genotype calls
missing_matrix <- is.na(genotype_matrix) |
  genotype_matrix == "--"

# 4. Calculate missingness per sample
sample_missing_rate <- colMeans(missing_matrix)

# 5. Filter samples with <= 15% missingness
keep_samples <- sample_missing_rate <= sample_threshold

sample_filtered <- genotype_matrix[
  ,
  keep_samples,
  drop = FALSE
]

# 6. Recalculate SNP missingness
filtered_missing <- is.na(sample_filtered) |
  sample_filtered == "--"

snp_missing_rate <- rowMeans(filtered_missing)

# 7. Filter SNPs with <= 10% missingness
keep_snps <- snp_missing_rate <= snp_threshold

clean_genotypes <- sample_filtered[
  keep_snps,
  ,
  drop = FALSE
]

# 8. Align accession metadata
clean_metadata <- rice_metadata[
  match(
    colnames(clean_genotypes),
    rice_metadata$Accession
  ),
  ,
  drop = FALSE
]

# 9. Verify metadata alignment
stopifnot(
  identical(
    colnames(clean_genotypes),
    clean_metadata$Accession
  )
)

# 10. Display QC results
cat("Original SNPs:", nrow(genotype_matrix), "\n")
cat("Original samples:", ncol(genotype_matrix), "\n")

cat(
  "Samples retained:",
  ncol(clean_genotypes),
  "\n"
)

cat(
  "SNPs retained:",
  nrow(clean_genotypes),
  "\n"
)

cat(
  "Remaining missing percentage:",
  round(
    mean(
      is.na(clean_genotypes) |
        clean_genotypes == "--"
    ) * 100,
    2
  ),
  "%\n"
)

# 11. Save cleaned data
saveRDS(
  clean_genotypes,
  "data/clean_genotypes.rds"
)

saveRDS(
  clean_metadata,
  "data/clean_metadata.rds"
)

cat("Revised QC completed successfully!\n")
