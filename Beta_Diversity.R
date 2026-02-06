# Data Preprocessing
tab <- read.table(
  "metaphlan_merged.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  comment.char = "#",
  check.names = FALSE
)

# keep species only
tab_species <- tab[grep("\\|s__", rownames(tab)), ]

dim(tab)
dim(tab_species)
head(rownames(tab_species), 5)

rownames(tab_species) <- sub(".*\\|s__", "s__", rownames(tab_species))
head(rownames(tab_species), 5)


colnames(tab_species)

library(readr)
metadata <- read_delim("metadata.csv", delim = ";", 
                       escape_double = FALSE, trim_ws = TRUE)
View(metadata)

rownames(metadata)

metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$Sampleid

common_samples <- intersect(
  colnames(tab_species),
  rownames(metadata)
)


length(common_samples)


tab_species <- tab_species[, common_samples, drop = FALSE]
metadata    <- metadata[common_samples, , drop = FALSE]

all(colnames(tab_species) == rownames(metadata))

dim(tab_species)

library(phyloseq)

OTU <- otu_table(
  as.matrix(tab_species),
  taxa_are_rows = TRUE
)

ps <- phyloseq(
  OTU,
  sample_data(metadata)
)


library(phyloseq)
library(vegan)

ps <- phyloseq(
  otu_table(as.matrix(tab_species), taxa_are_rows = TRUE),
  sample_data(metadata)
)

# Bray–Curtis distance
dist_bc <- phyloseq::distance(ps, method = "bray")
bc_mat <- as.matrix(dist_bc)

# Healthy controls
meta_df <- as(sample_data(ps), "data.frame")
healthy_ids <- rownames(meta_df[meta_df$Group == "Healthy", ])