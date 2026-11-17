
# Rice Genetic Diversity Analysis
# Step 12: Hierarchical Clustering

library(ggplot2)

# 1. Load data
genotypes <- readRDS(
  "data/polymorphic_genotypes.rds"
)

metadata <- readRDS(
  "data/clean_metadata.rds"
)

# 2. Verify sample order
stopifnot(
  identical(
    colnames(genotypes),
    as.character(metadata$Accession)
  )
)

# 3. Transpose genotype matrix
# Rows = rice accessions
# Columns = SNP markers
cluster_matrix <- t(genotypes)

# 4. Impute missing values using SNP mean dosage
for (j in seq_len(ncol(cluster_matrix))) {
  
  missing <- is.na(cluster_matrix[, j])
  
  if (any(missing)) {
    cluster_matrix[missing, j] <- mean(
      cluster_matrix[, j],
      na.rm = TRUE
    )
  }
}

stopifnot(!anyNA(cluster_matrix))

# 5. Calculate Euclidean genetic distances
genetic_distance <- dist(
  cluster_matrix,
  method = "euclidean"
)

# 6. Perform hierarchical clustering
hc <- hclust(
  genetic_distance,
  method = "ward.D2"
)

# 7. Save dendrogram
png(
  "figures/rice_hierarchical_clustering.png",
  width = 3000,
  height = 1800,
  res = 300
)

plot(
  hc,
  labels = FALSE,
  main = "Hierarchical Clustering of Rice Accessions",
  xlab = "Rice Accessions",
  ylab = "Clustering Height",
  sub = ""
)

dev.off()

# 8. Create six exploratory clusters
cluster_assignments <- cutree(
  hc,
  k = 6
)

# 9. Create clustering results table
cluster_results <- data.frame(
  Accession = rownames(cluster_matrix),
  Group = metadata$`Group name`,
  Cluster = cluster_assignments
)

# 10. Save results
write.csv(
  cluster_results,
  "results/hierarchical_cluster_assignments.csv",
  row.names = FALSE
)

saveRDS(
  hc,
  "results/hierarchical_clustering.rds"
)

# 11. Display cluster counts
cat("Cluster sizes:\n")
print(table(cluster_results$Cluster))

cat("\nGroups within each cluster:\n")
print(table(
  cluster_results$Group,
  cluster_results$Cluster
))

cat("\nHierarchical clustering completed!\n")
