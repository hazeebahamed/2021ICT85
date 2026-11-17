
# Rice Genetic Diversity Analysis
# Step 4: SNP Preprocessing

library(readxl)

# 1. Import datasets
accession_data <- read_excel(
  "data/S1 Table.  List of 544 Oryza accessions .xlsx",
  skip = 1
)

snp_raw <- read_excel(
  "data/S4 Table.  SNP Genotype Data Matrix .xlsx",
  skip = 1,
  col_types = "text"
)

# 2. Filter Oryza sativa accessions
rice_accessions <- accession_data[
  accession_data$Species == "O.sativa",
]

# 3. Remove the non-SNP Height row
snp_genotypes <- snp_raw[
  !is.na(snp_raw$rs) &
    snp_raw$rs != "Height",
]

# 4. Identify genotype sample columns
snp_sample_ids <- names(snp_genotypes)[12:ncol(snp_genotypes)]

# 5. Find matching rice accessions
matched_ids <- intersect(
  rice_accessions$Accession,
  snp_sample_ids
)

# 6. Select matched accession information
rice_metadata <- rice_accessions[
  match(matched_ids, rice_accessions$Accession),
]

# 7. Extract genotype matrix
genotype_data <- snp_genotypes[
  ,
  matched_ids,
  drop = FALSE
]

# 8. Convert to matrix
genotype_matrix <- as.matrix(genotype_data)

# 9. Check dimensions
print(dim(genotype_matrix))
print(nrow(rice_metadata))

# 10. Check SNP genotype values
print(table(
  as.vector(genotype_matrix),
  useNA = "ifany"
))

# 11. Save processed data
saveRDS(
  genotype_matrix,
  "data/genotype_matrix.rds"
)

saveRDS(
  rice_metadata,
  "data/rice_metadata.rds"
)

print("SNP preprocessing completed successfully!")
