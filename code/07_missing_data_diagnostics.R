
# Rice Genetic Diversity Analysis
# Step 7: Missing Data Diagnostics

genotype_matrix <- readRDS(
  "data/genotype_matrix.rds"
)

rice_metadata <- readRDS(
  "data/rice_metadata.rds"
)

# Identify missing genotype calls
missing_matrix <- is.na(genotype_matrix) |
  genotype_matrix == "--"

# Calculate sample and SNP missingness
sample_missing_rate <- colMeans(missing_matrix)
snp_missing_rate <- rowMeans(missing_matrix)

# Summarise sample missingness
cat("\nSample missingness summary (%):\n")
print(summary(sample_missing_rate * 100))

# Summarise SNP missingness
cat("\nSNP missingness summary (%):\n")
print(summary(snp_missing_rate * 100))

# Compare different sample QC thresholds
thresholds <- c(0.10, 0.15, 0.20, 0.25, 0.30)

qc_comparison <- data.frame(
  Threshold = paste0(thresholds * 100, "%"),
  Samples_Retained = sapply(
    thresholds,
    function(t) sum(sample_missing_rate <= t)
  ),
  SNPs_Retained_Before_Sample_Filter = sapply(
    thresholds,
    function(t) sum(snp_missing_rate <= t)
  )
)

cat("\nQC threshold comparison:\n")
print(qc_comparison)

# Missingness by rice group
group_missing <- data.frame(
  Accession = names(sample_missing_rate),
  Missing_Percentage = as.numeric(
    sample_missing_rate * 100
  )
)

group_missing$Group <- rice_metadata$`Group name`[
  match(
    group_missing$Accession,
    rice_metadata$Accession
  )
]

cat("\nMissingness by rice group:\n")

print(
  aggregate(
    Missing_Percentage ~ Group,
    data = group_missing,
    FUN = function(x) round(median(x), 2)
  )
)
