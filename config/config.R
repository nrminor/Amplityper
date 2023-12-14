# name of the FASTA containing gene sequences
gene_fasta <- "sars-cov-2_genes.fasta"

# desired output name
output_name <- "genes.bed"

# set the columns for the data frame that will be made by splitting
# up the header
columns <- c(
  "gene number",
  "gene",
  "tag",
  "db_xref",
  "location",
  "gbkey"
)