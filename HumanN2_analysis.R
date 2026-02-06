# HumanN2 Analysis
library(tidyverse)
library(pheatmap)
library(ggplot2)

# Load HUMAnN2 pathway abundance table
pathway <- read.table(
  "gene_pathabundance_relab.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE,
  row.names = 1,
  comment.char = ""   # <-- CRITICAL FIX
)

head(colnames(pathway))

colnames(pathway) <- sub("^HumanN2_", "", colnames(pathway))
head(colnames(pathway))


metadata <- read.table(
  "metadata.tsv",
  header = TRUE,
  sep = "\t",
  check.names = FALSE
)

metadata <- as.data.frame(metadata)
rownames(metadata) <- metadata$Sampleid
metadata$Sampleid <- NULL

common_samples <- intersect(colnames(pathway), rownames(metadata))

length(common_samples)     # MUST be > 0
common_samples

pathway  <- pathway[, common_samples, drop = FALSE]
metadata <- metadata[common_samples, , drop = FALSE]

all(colnames(pathway) == rownames(metadata))

library(pheatmap)

pathway <- read.table(
  "gene_pathabundance_relab.tsv",
  header = TRUE,
  sep = "\t",
  row.names = 1,
  check.names = FALSE,
  comment.char = ""   # CRITICAL: keep header line starting with '#'
)

colnames(pathway) <- sub("^HumanN2_", "", colnames(pathway))

pathway <- pathway[!grepl("^UNMAPPED|^UNINTEGRATED", rownames(pathway)), ]

top_features <- names(
  sort(apply(pathway, 1, var), decreasing = TRUE)
)[1:30]

#khusus yang gene family
# Remove unknown features (case-insensitive)
pathway_clean <- pathway[!grepl("unknown", rownames(pathway), ignore.case = TRUE), ]

# Prevalence filter
prev <- rowSums(pathway_clean > 0) / ncol(pathway_clean)
pathway_filt <- pathway_clean[prev >= 0.2, ]

# Log transform
pathway_log <- log10(pathway_filt + 1e-6)

# Select top variable pathways
top_features <- names(
  sort(apply(pathway_log, 1, var), decreasing = TRUE)
)[1:30]

pathway_clean <- pathway_clean[
  !grepl("^UNMAPPED|^UNINTEGRATED", rownames(pathway_clean)),
]

mat_scaled <- mat_scaled[!grepl("unknown", rownames(mat_scaled), ignore.case = TRUE), ]


pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_cols = FALSE,
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  show_colnames = FALSE,
  fontsize_row = 7,
  border_color = NA,
  main = "HUMAnN2 Pathway Heatmap (Annotated Features Only)"
)


#######

mat <- as.matrix(pathway[top_features, ])
storage.mode(mat) <- "numeric"

# Row-wise Z-score
mat_scaled <- t(scale(t(mat)))

# Remove rows that become NA after scaling
mat_scaled <- mat_scaled[complete.cases(mat_scaled), ]

stopifnot(nrow(mat_scaled) > 1)

annotation_col <- data.frame(
  Group = factor(metadata$Group)
)

rownames(annotation_col) <- rownames(metadata)
annotation_col <- annotation_col[colnames(mat_scaled), , drop = FALSE]

stopifnot(all(rownames(annotation_col) == colnames(mat_scaled)))

ann_colors <- list(
  Group = c(
    "Healthy"  = "#4DAF4A",
    "Diarrhea" = "#E41A1C"
  )
)

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_colnames = FALSE,
  border_color = NA,
  fontsize_row = 8,
  main = "HUMAnN2 Pathway Heatmap (Top Variable Pathways)"
)

prev <- rowSums(pathway > 0) / ncol(pathway)
pathway_filt <- pathway[prev >= 0.2, ]

pathway_log <- log10(pathway_filt + 1e-6)

top_features <- names(
  sort(apply(pathway_log, 1, var), decreasing = TRUE)
)[1:30]

mat <- as.matrix(pathway_log[top_features, ])
mat_scaled <- t(scale(t(mat)))
mat_scaled <- mat_scaled[complete.cases(mat_scaled), ]

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_colnames = FALSE,
  show_rownames = TRUE,
  fontsize_row = 7,
  border_color = NA,
  legend = TRUE,
  main = "HUMAnN2 Pathway Heatmap (Filtered & Log-transformed)"
)

####

metadata$Group <- factor(
  metadata$Group,
  levels = c("Healthy", "Diarrhea")
)

sample_order <- rownames(metadata)[order(metadata$Group)]

sample_order

mat_scaled <- mat_scaled[, sample_order, drop = FALSE]
annotation_col <- annotation_col[sample_order, , drop = FALSE]

all(colnames(mat_scaled) == rownames(annotation_col))
# must be TRUE

pheatmap(
  mat_scaled,
  annotation_col = annotation_col,
  annotation_colors = ann_colors,
  cluster_cols = FALSE,      # IMPORTANT: preserve group order
  clustering_distance_rows = "euclidean",
  clustering_method = "complete",
  show_colnames = FALSE,
  fontsize_row = 7,
  border_color = NA,
  main = "HUMAnN2 Pathway Heatmap (Samples Ordered by Group)"
)


# Look at the first few lines of the HUMAnN file
readLines("humann_pathabundance_relab.tsv", n = 5)

# Look at the first few lines of metadata
readLines("metadata.tsv", n = 5)


pathway <- pathway[!grepl("UNMAPPED|UNINTEGRATED", rownames(pathway)), ]

# Load Metadata
metadata <- read.table(
  "metadata.tsv",
  header = TRUE,
  row.names = 1,
  sep = "\t"
)

common_samples <- intersect(colnames(pathway), rownames(metadata))
pathway <- pathway[, common_samples]
metadata <- metadata[common_samples, ]

colnames(pathway)
rownames(metadata)

# Stacked barplot (TOP pathways)
topN <- 20

path_long <- pathway %>%
  as.data.frame() %>%
  rownames_to_column("Pathway") %>%
  pivot_longer(-Pathway, names_to = "Sample", values_to = "Abundance")

top_pathways <- path_long %>%
  group_by(Pathway) %>%
  summarise(mean_abundance = mean(Abundance)) %>%
  arrange(desc(mean_abundance)) %>%
  slice_head(n = topN) %>%
  pull(Pathway)

path_long %>%
  filter(Pathway %in% top_pathways) %>%
  left_join(
    metadata %>% rownames_to_column("Sample"),
    by = "Sample"
  ) %>%
  ggplot(aes(x = Sample, y = Abundance, fill = Pathway)) +
  geom_bar(stat = "identity") +
  facet_grid(~ Group, scales = "free_x", space = "free_x") +
  theme_classic() +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank()
  ) +
  labs(
    title = "Top Functional Pathways (HUMAnN2)",
    y = "Relative abundance",
    x = "Samples"
  )

#Heatmap (top variable pathways)
top_var <- apply(pathway, 1, var)
top_features <- names(sort(top_var, decreasing = TRUE))[1:30]

annotation_col <- data.frame(
  Group = metadata$Group
)
rownames(annotation_col) <- rownames(metadata)

pheatmap(
  pathway[top_features, ],
  scale = "row",
  annotation_col = annotation_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "bray",
  clustering_method = "complete",
  show_colnames = FALSE
)

pheatmap(
  pathway[top_features, ],
  scale = "row",
  annotation_col = annotation_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_colnames = FALSE
)


mat <- pathway[top_features, ]
# Remove rows with zero variance
mat <- mat[apply(mat, 1, sd, na.rm = TRUE) > 0, ]

dim(mat)

mat <- as.matrix(mat)
storage.mode(mat) <- "numeric"

annotation_col <- data.frame(
  Group = metadata$Group
)
rownames(annotation_col) <- rownames(metadata)
all(rownames(annotation_col) == colnames(mat))

pheatmap(
  mat,
  scale = "row",
  annotation_col = annotation_col,
  clustering_distance_rows = "euclidean",
  clustering_distance_cols = "euclidean",
  clustering_method = "complete",
  show_colnames = FALSE,
  na_col = "grey90"
)


mat0 <- pathway[top_features, ]

# Force numeric
mat0 <- as.matrix(mat0)
storage.mode(mat0) <- "numeric"

mat_scaled <- t(scale(t(mat0)))

mat_scaled <- mat_scaled[complete.cases(mat_scaled), ]

dim(mat_scaled)

annotation_col <- data.frame(
  Group = factor(metadata$Group)
)

rownames(annotation_col) <- rownames(metadata)

# Align order
annotation_col <- annotation_col[colnames(mat_scaled), , drop = FALSE]

# Sanity check
all(rownames(annotation_col) == colnames(mat_scaled))

# Make sure order matches matrix columns
annotation_col <- annotation_col[colnames(mat), , drop = FALSE]

#Boxplot for specific pathways (Disease vs Healthy)
target_pathways <- c(
  "PWY-6572",   # example
  "PWY-6386"
)

path_long %>%
  filter(Pathway %in% target_pathways) %>%
  left_join(
    metadata %>% rownames_to_column("Sample"),
    by = "Sample"
  ) %>%
  ggplot(aes(x = Group, y = Abundance, fill = Group)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.6) +
  facet_wrap(~ Pathway, scales = "free_y") +
  theme_classic() +
  labs(
    y = "Relative abundance",
    x = ""
  )

#PCA / PCoA on functional profiles (beta diversity)
library(vegan)

dist_func <- vegdist(t(pathway), method = "bray")
ord_func <- cmdscale(dist_func, k = 2)

ord_df <- data.frame(
  ord_func,
  Group = metadata$Group
)

ggplot(ord_df, aes(X1, X2, color = Group)) +
  geom_point(size = 3) +
  theme_classic() +
  labs(
    x = "PCoA1",
    y = "PCoA2",
    title = "Functional Beta Diversity (HUMAnN2)"
  )
