
# Rice Genetic Diversity Analysis
# Step 8: SNP Encoding

library(readxl)

# 1. Load quality-controlled genotypes
clean_genotypes <- readRDS(
  "data/clean_genotypes.rds"
)

clean_metadata <- readRDS(
  "data/clean_metadata.rds"
)

# 2. Load original SNP dataset
snp_raw <- read_excel(
  "data/S4 Table.  SNP Genotype Data Matrix .xlsx",
  skip = 1,
  col_types = "text"
)

# 3. Remove Height row
snp_raw <- snp_raw[
  !is.na(snp_raw$rs) &
    snp_raw$rs != "Height",
]

# 4. Identify matching SNP rows
# QC preserved the original SNP row order
original_ids <- rownames(clean_genotypes)

# If rownames are default row numbers,
# we must NOT use them as biological SNP IDs.
print(head(original_ids))

# 5. Inspect SNP allele definitions
print(head(snp_raw[, c("rs", "alleles")]))

# 6. Check whether SNP identifiers are unique
print(sum(duplicated(snp_raw$rs)))

# 7. Inspect allele definition formats
print(table(snp_raw$alleles, useNA = "ifany"))
