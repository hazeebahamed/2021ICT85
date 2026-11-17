
# Rice Genetic Diversity Analysis
# Step 10: Minor Allele Frequency Analysis

library(ggplot2)

# 1. Load encoded genotype data
numeric_genotypes <- readRDS(
  "data/numeric_genotypes.rds"
)

# 2. Calculate allele frequency
# Rows = SNPs, columns = rice accessions
allele_frequency <- rowMeans(
  numeric_genotypes,
  na.rm = TRUE
) / 2

# 3. Calculate Minor Allele Frequency
maf <- pmin(
  allele_frequency,
  1 - allele_frequency
)

# 4. Remove SNPs with no observed genotype calls
maf[is.nan(maf)] <- NA_real_

# 5. Identify monomorphic and polymorphic SNPs
monomorphic <- sum(maf == 0, na.rm = TRUE)
polymorphic <- sum(maf > 0, na.rm = TRUE)

cat("Total SNPs:", length(maf), "\n")
cat("Monomorphic SNPs:", monomorphic, "\n")
cat("Polymorphic SNPs:", polymorphic, "\n")

# 6. Display MAF summary
print(summary(maf))

# 7. Create MAF results table
maf_results <- data.frame(
  SNP_ID = rownames(numeric_genotypes),
  MAF = maf
)

# 8. Save results
write.csv(
  maf_results,
  "results/snp_maf_results.csv",
  row.names = FALSE
)

# 9. Create MAF distribution histogram
maf_plot <- ggplot(
  maf_results,
  aes(x = MAF)
) +
  geom_histogram(
    bins = 40,
    fill = "steelblue",
    color = "white"
  ) +
  labs(
    title = "Minor Allele Frequency Distribution in Rice",
    x = "Minor Allele Frequency",
    y = "Number of SNP Markers"
  ) +
  theme_minimal(base_size = 12)

print(maf_plot)

# 10. Save figure
ggsave(
  "figures/maf_distribution.png",
  plot = maf_plot,
  width = 8,
  height = 5,
  dpi = 300
)

# 11. Retain polymorphic SNPs for PCA
polymorphic_genotypes <- numeric_genotypes[
  !is.na(maf) & maf > 0,
  ,
  drop = FALSE
]

saveRDS(
  polymorphic_genotypes,
  "data/polymorphic_genotypes.rds"
)

cat(
  "SNPs retained for PCA:",
  nrow(polymorphic_genotypes),
  "\n"
)

cat("MAF analysis completed successfully!\n")
