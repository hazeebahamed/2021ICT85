
# Rice Genetic Diversity Analysis
# Step 2: Data Import and Inspection

# Load required package
library(readxl)

# Define dataset file paths
accession_file <- file.path(
  "data",
  "S1 Table.  List of 544 Oryza accessions .xlsx"
)

snp_file <- file.path(
  "data",
  "S4 Table.  SNP Genotype Data Matrix .xlsx"
)

# Check that both files exist
print(file.exists(accession_file))
print(file.exists(snp_file))

# Display Excel sheet names
print(excel_sheets(accession_file))
print(excel_sheets(snp_file))

# Import both datasets
accession_data <- read_excel(accession_file)
snp_data <- read_excel(snp_file)

# Display dataset dimensions
print(dim(accession_data))
print(dim(snp_data))

# Display first six rows
print(head(accession_data))
print(head(snp_data))


# Step 3: Inspect original Excel table structure

# Read first 6 rows of accession file
accession_preview <- read_excel(
  accession_file,
  col_names = FALSE,
  n_max = 6
)

# Read first 6 rows and 15 columns of SNP file
snp_preview <- read_excel(
  snp_file,
  col_names = FALSE,
  range = cell_limits(c(1, 1), c(6, 15))
)

# Display previews
print(accession_preview, width = Inf)
print(snp_preview, width = Inf)



# Step 4: Import SNP data with correct headers

# Skip the first descriptive row
snp_raw <- read_excel(
  snp_file,
  skip = 1,
  col_types = "text"
)

# Check corrected column names
print(head(names(snp_raw), 15))

# Check dimensions
print(dim(snp_raw))

# Inspect first five rows
print(snp_raw[1:5, 1:15])



# Step 5: Import accession data with correct headers

accession_raw <- read_excel(
  accession_file,
  skip = 1
)

# Display corrected column names
print(names(accession_raw))

# Display dataset dimensions
print(dim(accession_raw))

# Preview first six rows
print(head(accession_raw))

