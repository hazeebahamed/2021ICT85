
# Rice Genetic Diversity Analysis
# Step 13: Clustering Evaluation

# 1. Load clustering results
cluster_results <- read.csv(
  "results/hierarchical_cluster_assignments.csv"
)

# 2. Create contingency table
cluster_table <- table(
  cluster_results$Group,
  cluster_results$Cluster
)

print(cluster_table)

# 3. Calculate cluster purity
cluster_purity <- sum(
  apply(cluster_table, 2, max)
) / sum(cluster_table)

cat(
  "\nCluster purity:",
  round(cluster_purity * 100, 2),
  "%\n"
)

# 4. Calculate group-wise dominant cluster percentage
group_consistency <- apply(
  cluster_table,
  1,
  function(x) {
    max(x) / sum(x) * 100
  }
)

cat("\nGroup consistency (%):\n")
print(round(group_consistency, 2))

# 5. Save evaluation results
write.csv(
  data.frame(
    Group = names(group_consistency),
    Dominant_Cluster_Percentage =
      as.numeric(group_consistency)
  ),
  "results/group_cluster_consistency.csv",
  row.names = FALSE
)

cat("\nCluster evaluation completed!\n")
