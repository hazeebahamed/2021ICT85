
# Rice Genetic Diversity Analysis
# Step 9: Numerical SNP Encoding

library(readxl)

# 1. Load original and QC datasets
genotype_matrix <- readRDS("data/genotype_matrix.rds")
clean_genotypes <- readRDS("data/clean_genotypes.rds")
clean_metadata <- readRDS("data/clean_metadata.rds")

# 2. Load SNP information
snp_raw <- read_excel(
  "data/S4 Table.  SNP Genotype Data Matrix .xlsx",
  skip = 1,
  col_types = "text"
)

snp_raw <- snp_raw[
  !is.na(snp_raw$rs) & snp_raw$rs != "Height",
]

# 3. Reconstruct the SNP QC selection
sample_filtered <- genotype_matrix[
  ,
  colnames(clean_genotypes),
  drop = FALSE
]

missing_matrix <- is.na(sample_filtered) |
  sample_filtered == "--"

keep_snps <- rowMeans(missing_matrix) <= 0.10

# Verify the same SNPs were retained
stopifnot(
  sum(keep_snps) == nrow(clean_genotypes)
)

# Verify that filtered genotypes match
stopifnot(
  identical(
    unname(sample_filtered[keep_snps, , drop = FALSE]),
    unname(clean_genotypes)
  )
)

# 4. Extract matching SNP metadata
snp_info <- snp_raw[keep_snps, c("rs", "alleles")]

# 5. Prepare numerical genotype matrix
numeric_genotypes <- matrix(
  NA_real_,
  nrow = nrow(clean_genotypes),
  ncol = ncol(clean_genotypes),
  dimnames = list(
    snp_info$rs,
    colnames(clean_genotypes)
  )
)


# 6. Encode genotypes with strand-aware allele matching

complement_base <- function(base) {
  complement <- c(
    A = "T",
    T = "A",
    C = "G",
    G = "C"
  )
  
  unname(complement[base])
}

for (i in seq_len(nrow(clean_genotypes))) {
  
  alleles <- strsplit(
    snp_info$alleles[i],
    "/",
    fixed = TRUE
  )[[1]]
  
  if (length(alleles) != 2) {
    stop(paste("Invalid alleles at SNP", snp_info$rs[i]))
  }
  
  a1 <- alleles[1]
  a2 <- alleles[2]
  
  # Original allele coding
  forward_codes <- c(
    paste0(a1, a1),
    paste0(a1, a2),
    paste0(a2, a1),
    paste0(a2, a2)
  )
  
  # Complementary allele coding
  c1 <- complement_base(a1)
  c2 <- complement_base(a2)
  
  reverse_codes <- c(
    paste0(c1, c1),
    paste0(c1, c2),
    paste0(c2, c1),
    paste0(c2, c2)
  )
  
  observed <- as.character(
    unlist(clean_genotypes[i, ])
  )
  
  valid_observed <- unique(
    observed[
      !is.na(observed) &
        observed != "--"
    ]
  )
  
  # Choose encoding matching the observed genotypes
  if (all(valid_observed %in% forward_codes)) {
    
    selected_codes <- forward_codes
    
  } else if (all(valid_observed %in% reverse_codes)) {
    
    selected_codes <- reverse_codes
    
  } else {
    
    stop(paste(
      "Unrecognised genotype coding at SNP",
      snp_info$rs[i]
    ))
  }
  
  dosage_values <- c(0, 1, 1, 2)
  
  numeric_genotypes[i, ] <- dosage_values[
    match(observed, selected_codes)
  ]
}

# 7. Display encoding summary
cat("Encoded SNPs:", nrow(numeric_genotypes), "\n")
cat("Encoded samples:", ncol(numeric_genotypes), "\n")

cat("\nNumerical genotype counts:\n")
print(table(
  as.vector(numeric_genotypes),
  useNA = "ifany"
))

# 8. Save encoded data and SNP information
saveRDS(
  numeric_genotypes,
  "data/numeric_genotypes.rds"
)

saveRDS(
  snp_info,
  "data/filtered_snp_info.rds"
)

cat("\nSNP encoding completed successfully!\n")
