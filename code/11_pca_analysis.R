
# Rice Genetic Diversity Analysis
# Step 11: Principal Component Analysis

library(ggplot2)

# 1. Load data
polymorphic_genotypes <- readRDS(
  "data/polymorphic_genotypes.rds"
)

clean_metadata <- readRDS(
  "data/clean_metadata.rds"
)

# 2. Check sample alignment
stopifnot(
  identical(
    colnames(polymorphic_genotypes),
    as.character(clean_metadata$Accession)
  )
)

# 3. Transpose matrix
# Rows = rice accessions
# Columns = SNP markers
pca_matrix <- t(polymorphic_genotypes)

# 4. Replace missing genotypes with
# the mean dosage of each SNP
for (j in seq_len(ncol(pca_matrix))) {
  
  missing <- is.na(pca_matrix[, j])
  
  if (any(missing)) {
    pca_matrix[missing, j] <- mean(
      pca_matrix[, j],
      na.rm = TRUE
    )
  }
}

# 5. Perform PCA
pca_result <- prcomp(
  pca_matrix,
  center = TRUE,
  scale. = FALSE
)

# 6. Calculate explained variance
variance_explained <- 100 *
  pca_result$sdev^2 /
  sum(pca_result$sdev^2)

# 7. Create PCA results table
pca_data <- data.frame(
  Accession = rownames(pca_matrix),
  PC1 = pca_result$x[, 1],
  PC2 = pca_result$x[, 2],
  Group = clean_metadata$`Group name`
)

# 8. Create PCA plot
pca_plot <- ggplot(
  pca_data,
  aes(
    x = PC1,
    y = PC2,
    color = Group
  )
) +
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  labs(
    title = "PCA of Rice Genetic Diversity",
    x = paste0(
      "PC1 (",
      round(variance_explained[1], 2),
      "%)"
    ),
    y = paste0(
      "PC2 (",
      round(variance_explained[2], 2),
      "%)"
    ),
    color = "Rice Group"
  ) +
  theme_minimal(base_size = 12)

# 9. Display plot
print(pca_plot)

# 10. Save figure
ggsave(
  "figures/rice_pca_plot.png",
  plot = pca_plot,
  width = 10,
  height = 6,
  dpi = 300
)

# 11. Save PCA results
write.csv(
  pca_data,
  "results/pca_coordinates.csv",
  row.names = FALSE
)

saveRDS(
  pca_result,
  "results/pca_result.rds"
)

# 12. Display results
cat("PCA completed successfully!\n")
cat(
  "PC1 variance:",
  round(variance_explained[1], 2),
  "%\n"
)
cat(
  "PC2 variance:",
  round(variance_explained[2], 2),
  "%\n"
)
