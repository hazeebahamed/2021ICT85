
# Rice Genetic Diversity Analysis
# Step 6: QC Visualization

library(ggplot2)

# Load the original genotype matrix
genotype_matrix <- readRDS(
  "data/genotype_matrix.rds"
)

# Calculate missingness
missing_matrix <- is.na(genotype_matrix) |
  genotype_matrix == "--"

snp_missing_rate <- rowMeans(missing_matrix)

# Prepare plotting data
qc_data <- data.frame(
  Missing_Percentage = snp_missing_rate * 100
)

# Create histogram
p <- ggplot(
  qc_data,
  aes(x = Missing_Percentage)
) +
  geom_histogram(
    binwidth = 2,
    fill = "#2E7D32",
    color = "white"
  ) +
  geom_vline(
    xintercept = 10,
    color = "red",
    linetype = "dashed",
    linewidth = 1
  ) +
  labs(
    title = "Missing Genotype Calls Across Rice SNP Markers",
    x = "Missing Genotypes (%)",
    y = "Number of SNP Markers",
    caption = "Dashed line indicates the 10% missingness threshold"
  ) +
  theme_minimal(base_size = 12)

# Display figure
print(p)

# Save figure
ggsave(
  "figures/snp_missingness_histogram.png",
  plot = p,
  width = 9,
  height = 6,
  dpi = 300
)
