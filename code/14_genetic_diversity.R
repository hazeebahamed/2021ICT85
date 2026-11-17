
# Rice Genetic Diversity Analysis
# Step 14: Genetic Diversity Comparison

library(ggplot2)

# 1. Load data
genotypes <- readRDS(
  "data/polymorphic_genotypes.rds"
)

metadata <- readRDS(
  "data/clean_metadata.rds"
)

# 2. Verify sample alignment
stopifnot(
  identical(
    colnames(genotypes),
    as.character(metadata$Accession)
  )
)

# 3. Identify rice groups
groups <- unique(metadata$`Group name`)

# 4. Calculate expected heterozygosity
diversity_results <- data.frame(
  Group = character(),
  Sample_Size = integer(),
  Mean_Expected_Heterozygosity = numeric()
)

for (group in groups) {
  
  sample_ids <- metadata$Accession[
    metadata$`Group name` == group
  ]
  
  group_genotypes <- genotypes[
    ,
    sample_ids,
    drop = FALSE
  ]
  
  # Allele frequency for each SNP
  p <- rowMeans(
    group_genotypes,
    na.rm = TRUE
  ) / 2
  
  # Expected heterozygosity
  he <- 2 * p * (1 - p)
  
  he[is.nan(he)] <- NA_real_
  
  diversity_results <- rbind(
    diversity_results,
    data.frame(
      Group = group,
      Sample_Size = length(sample_ids),
      Mean_Expected_Heterozygosity =
        mean(he, na.rm = TRUE)
    )
  )
}

# 5. Display results
print(diversity_results)

# 6. Save results
write.csv(
  diversity_results,
  "results/genetic_diversity_comparison.csv",
  row.names = FALSE
)

# 7. Create comparison plot
diversity_plot <- ggplot(
  diversity_results,
  aes(
    x = reorder(
      Group,
      Mean_Expected_Heterozygosity
    ),
    y = Mean_Expected_Heterozygosity,
    fill = Group
  )
) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Genetic Diversity Among Rice Groups",
    x = "Rice Group",
    y = "Mean Expected Heterozygosity"
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "none")

print(diversity_plot)

# 8. Save figure
ggsave(
  "figures/genetic_diversity_comparison.png",
  plot = diversity_plot,
  width = 9,
  height = 6,
  dpi = 300
)

cat("\nGenetic diversity analysis completed!\n")
