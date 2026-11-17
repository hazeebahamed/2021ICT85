
# Rice Genetic Diversity Analysis
# Step 3: Data Cleaning

library(readxl)

# Import accession information
accession_data <- read_excel(
  "data/S1 Table.  List of 544 Oryza accessions .xlsx",
  skip = 1
)

# Import SNP genotype data
snp_raw <- read_excel(
  "data/S4 Table.  SNP Genotype Data Matrix .xlsx",
  skip = 1,
  col_types = "text"
)

# Check dimensions before cleaning
print(dim(accession_data))
print(dim(snp_raw))

# Inspect species distribution
print(table(accession_data$Species, useNA = "ifany"))

# Inspect rice group distribution
print(table(accession_data$`Group name`, useNA = "ifany"))

# Check accession identifiers
print(head(accession_data$Accession))

# Check SNP identifiers
print(head(snp_raw$rs))

# Check missing accession names
print(sum(is.na(accession_data$Accession)))

# Check duplicate accession names
print(sum(duplicated(accession_data$Accession)))


# Step 4: Filter Oryza sativa accessions

rice_accessions <- accession_data[
  accession_data$Species == "O.sativa",
]

# Check number of selected accessions
print(nrow(rice_accessions))

# Check group distribution
print(table(rice_accessions$`Group name`))

# Check first six accessions
print(head(rice_accessions$Accession))


# Step 5: Inspect SNP sample identifiers

# First 11 columns contain SNP metadata
snp_metadata_columns <- names(snp_raw)[1:11]

# Remaining columns contain sample genotypes
snp_sample_ids <- names(snp_raw)[12:ncol(snp_raw)]

# Number of sample columns
print(length(snp_sample_ids))

# First 20 SNP sample IDs
print(head(snp_sample_ids, 20))

# First 20 accession names
print(head(accession_data$Accession, 20))

# Count exact matches between the identifiers
print(sum(snp_sample_ids %in% accession_data$Accession))


# Step 6: Identify unmatched sample identifiers

# SNP sample IDs not found in accession information
unmatched_snp <- setdiff(
  snp_sample_ids,
  accession_data$Accession
)

# Accession IDs not found in SNP matrix
unmatched_accessions <- setdiff(
  accession_data$Accession,
  snp_sample_ids
)

# Display unmatched identifiers
print(unmatched_snp)
print(unmatched_accessions)

# Check duplicate SNP sample identifiers
print(sum(duplicated(snp_sample_ids)))

# Count matched Oryza sativa accessions
matched_rice <- rice_accessions[
  rice_accessions$Accession %in% snp_sample_ids,
]

print(nrow(matched_rice))


