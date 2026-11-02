# ============================================================
# RSCU CODON USAGE ANALYSIS
# ============================================================
#
# PURPOSE:
# Compare synonymous codon usage across five genomic sources:
#   1. Mega
#   2. Plasmid
#   3. Chromosome
#   4. Fungi
#   5. Phage (subdivided into four groups)
#
# ANALYSES:
#   - PCA: overall similarity and separation among samples
#   - Correspondence analysis: relationships between codons
#     and samples
#   - Group-level RSCU summaries (mean, SD, CV)
#   - Codon contributions and quality of representation
#   - Individual codon comparisons
#   - RSCU distribution histograms
#
# INPUT:
#   Five CSV files containing gene/genome-level RSCU values.
#
# OUTPUT:
#   PCA_RSCU_plot.svg
#   RSCU_biplot.svg
#   RSCU_codon_usage_histograms.svg
#
# NOTE:
#   All analyses use codons shared across all five datasets.
#
# ============================================================


# ============================================================
# 1. LOAD REQUIRED PACKAGES
# ============================================================

library(ggplot2)       # Plotting
library(svglite)       # Export editable SVG figures
library(FactoMineR)    # Correspondence analysis
library(factoextra)    # Visualization of PCA and CA
library(dplyr)         # Data manipulation
library(gplots)        # Balloon plots
library(tidyr)         # Reshape data for plotting
library(MetBrewer)     # Color palettes


# ============================================================
# 2. INPUT FILES
# ============================================================

# Named list linking each biological group to its RSCU CSV.
#
# The names on the left are used throughout the analysis.
# The filenames on the right are the actual input files.

files <- list(
  Mega       = "bacterial_mega_RSCU_edited.csv",
  Plasmid    = "bacterial_plasmid_RSCU_edited.csv",
  Chromosome = "bacterial_chromsome_RSCU_edited.csv",
  Fungi      = "fungal_RSCU_table_edited.csv",
  Phage      = "phage_edited_RSCU_names.csv"
)

# Check that all input files exist before proceeding.

missing_files <- files[!file.exists(unlist(files))]

if (length(missing_files) > 0) {
  stop(
    "The following input files are missing: ",
    paste(missing_files, collapse = ", ")
  )
}


# ============================================================
# 3. LOAD AND STANDARDIZE RSCU DATA
# ============================================================

# Function to read an RSCU CSV and standardize its structure.
#
# Standard datasets:
#   - First column contains the sample identifier.
#   - Remaining columns contain codon RSCU values.
#
# Phage dataset:
#   - First column contains the phage group (1-4).
#   - Second column contains the phage/sample identifier.
#   - Remaining columns contain codon RSCU values.
#
# All datasets are standardized to include:
#   Sample = sample identifier
#   Group  = biological group

load_rscu <- function(file, group_name, phage = FALSE) {
  
  df <- read.csv(
    file,
    stringsAsFactors = FALSE,
    check.names = FALSE
  )
  
  if (phage) {
    
    # Phage-specific formatting:
    # First column identifies the phage subgroup.
    # Second column identifies the individual phage.
    
    df$Sample <- df[[2]]
    df$Group <- paste0("Phage_", df[[1]])
    
    # Remove the original group and sample columns.
    # The standardized Sample and Group columns are retained.
    
    df <- df[, -(1:2)]
    
  } else {
    
    # Rename the first column to Sample.
    
    colnames(df)[1] <- "Sample"
    
    # Assign the dataset's biological group.
    
    df$Group <- group_name
  }
  
  return(df)
}


# Load all five datasets into a named list.

data_list <- list(
  Mega       = load_rscu(files$Mega, "Chromid"),
  Plasmid    = load_rscu(files$Plasmid, "Plasmid"),
  Chromosome = load_rscu(files$Chromosome, "Chromosome"),
  Fungi      = load_rscu(files$Fungi, "Fungi"),
  Phage      = load_rscu(files$Phage, "Phage", phage = TRUE)
)


# ============================================================
# 4. IDENTIFY CODONS SHARED ACROSS DATASETS
# ============================================================

# Identify codon columns in each dataset.
#
# Sample and Group are metadata columns and are excluded.

codon_cols <- lapply(
  data_list,
  function(df) {
    setdiff(colnames(df), c("Sample", "Group"))
  }
)

# Find the intersection of codon columns across all datasets.
#
# Only codons present in EVERY dataset are retained.
# This ensures that all samples are compared using
# the same set of codons.

common_codons <- Reduce(intersect, codon_cols)

# Report the number of shared codons.

cat("Number of shared codons:", length(common_codons), "\n")

if (length(common_codons) == 0) {
  stop("No common codon columns were found.")
}


# ============================================================
# 5. COMBINE DATASETS
# ============================================================

# Define the columns to retain:
#   Sample, all shared codons, and Group.

final_columns <- c("Sample", common_codons, "Group")

# Subset each dataset to the shared columns and combine them.
#
# bind_rows() stacks the datasets vertically.
# Each row represents one sample.

combined <- bind_rows(
  lapply(data_list, function(df) {
    df[, final_columns, drop = FALSE]
  })
)

# Convert codon columns to numeric.
# This ensures they can be used in statistical analyses.

combined[common_codons] <- lapply(
  combined[common_codons],
  as.numeric
)

# Set Group as a categorical variable.

group_levels <- c(
  "Chromid", "Plasmid", "Chromosome", "Fungi",
  "Phage_1", "Phage_2", "Phage_3", "Phage_4"
)

combined$Group <- factor(
  combined$Group,
  levels = group_levels
)

# Assign sample identifiers as row names.
# This helps keep sample identities consistent in PCA and CA.

rownames(combined) <- make.unique(as.character(combined$Sample))

# Inspect the combined dataset.

cat("\nCombined dataset dimensions:\n")
print(dim(combined))

cat("\nNumber of samples per group:\n")
print(table(combined$Group))


# ============================================================
# 6. PREPARE THE RSCU MATRIX
# ============================================================

# Extract only the codon columns.
#
# Rows    = individual samples
# Columns = codons
#
# Metadata (Sample and Group) is excluded.

pca_matrix <- combined[, common_codons, drop = FALSE]

# Remove codons with zero variance.
#
# A codon with zero variance has the same RSCU value in every
# sample and therefore cannot contribute to PCA separation.
#
# The same filtered codon matrix will be used for CA.

variable_codons <- apply(
  pca_matrix,
  2,
  function(x) var(x, na.rm = TRUE) > 0
)

# Remove zero-variance codons.

pca_matrix <- pca_matrix[, variable_codons, drop = FALSE]

# Ensure that all retained values are finite.
# Missing or non-finite values must be addressed before PCA/CA.

if (any(!is.finite(as.matrix(pca_matrix)))) {
  stop("The RSCU matrix contains missing or non-finite values.")
}

if (ncol(pca_matrix) < 2) {
  stop("Fewer than two variable codons remain.")
}


# ============================================================
# 7. COLOR AND SHAPE ASSIGNMENTS
# ============================================================

# Generate eight colors from the Thomas palette.
#
# Each biological group receives a consistent color
# across all figures.

thomas_colors <- met.brewer(
  "Thomas",
  n = 8
)

color_map <- c(
  "Chromid"    = thomas_colors[1],
  "Plasmid"    = thomas_colors[2],
  "Chromosome" = thomas_colors[3],
  "Fungi"      = thomas_colors[4],
  "Phage_1"    = thomas_colors[5],
  "Phage_2"    = thomas_colors[6],
  "Phage_3"    = thomas_colors[7],
  "Phage_4"    = thomas_colors[8]
)

# Display the palette.

scales::show_col(thomas_colors)

# Assign point shapes to distinguish biological groups.
#
# Bacterial groups: filled circles
# Fungi:            diamonds
# Phage groups:     open circles

shape_map <- c(
  "Chromid"    = 16,
  "Plasmid"    = 16,
  "Chromosome" = 16,
  "Fungi"      = 18,
  "Phage_1"    = 1,
  "Phage_2"    = 1,
  "Phage_3"    = 1,
  "Phage_4"    = 1
)


# ============================================================
# 8. PRINCIPAL COMPONENT ANALYSIS (PCA)
# ============================================================

# PURPOSE:
# PCA reduces the dimensionality of the RSCU dataset.
#
# Each sample is represented by its codon usage profile.
# Samples with similar codon usage should be positioned
# closer together in the PCA plot.
#
# Rows    = samples
# Columns = codon RSCU values
#
# center = TRUE:
#   Centers each codon around its mean.
#
# scale. = TRUE:
#   Scales codons to unit variance so that codons with
#   different variability contribute comparably.

pca <- prcomp(
  pca_matrix,
  center = TRUE,
  scale. = TRUE
)

# Calculate the percentage of variance explained by each PC.

pca_var <- round(
  100 * pca$sdev^2 / sum(pca$sdev^2),
  1
)

# Prepare a plotting dataframe.
#
# Each point represents one sample.
# PC1 and PC2 are the sample coordinates along the first
# two principal components.

pca_df <- data.frame(
  Sample = rownames(pca_matrix),
  PC1    = pca$x[, 1],
  PC2    = pca$x[, 2],
  Group  = combined$Group
)

# Generate PCA plot.

plot_pca <- ggplot(
  pca_df,
  aes(
    x = PC1,
    y = PC2,
    color = Group,
    shape = Group
  )
) +
  geom_point(
    size = 3,
    alpha = 0.8
  ) +
  scale_color_manual(values = color_map) +
  scale_shape_manual(values = shape_map) +
  theme_classic() +
  labs(
    title = "PCA of Codon Usage (RSCU)",
    x = paste0("PC1 (", pca_var[1], "%)"),
    y = paste0("PC2 (", pca_var[2], "%)")
  )

# Display PCA plot.

print(plot_pca)

# Save PCA figure as an editable SVG.

ggsave(
  "PCA_RSCU_plot.svg",
  plot = plot_pca,
  device = svglite,
  width = 7,
  height = 6
)


# ============================================================
# 9. GROUP-LEVEL RSCU SUMMARIES
# ============================================================

# PURPOSE:
# Calculate average codon usage for each biological group.
#
# These averages summarize the overall codon usage profile
# of each group.
#
# Rows    = biological groups
# Columns = codons
# Values  = mean RSCU

cont_table <- combined %>%
  group_by(Group) %>%
  summarise(
    across(
      all_of(common_codons),
      ~ mean(.x, na.rm = TRUE)
    ),
    .groups = "drop"
  ) %>%
  as.data.frame()

# Use group names as row names.

rownames(cont_table) <- cont_table$Group
cont_table$Group <- NULL

# Inspect group-level mean RSCU values.

head(cont_table)


# ------------------------------------------------------------
# 9a. Calculate mean, standard deviation, and CV
# ------------------------------------------------------------

# For every codon within every group, calculate:
#
# Mean = average RSCU across samples.
#
# SD = standard deviation of RSCU across samples.
#
# CV = coefficient of variation (SD / mean).
#
# CV describes relative variability in codon usage.
# Higher CV indicates greater variation relative to the mean.

group_stats <- combined %>%
  group_by(Group) %>%
  summarise(
    across(
      all_of(common_codons),
      list(
        mean = ~ mean(.x, na.rm = TRUE),
        sd   = ~ sd(.x, na.rm = TRUE),
        cv   = ~ sd(.x, na.rm = TRUE) / mean(.x, na.rm = TRUE)
      ),
      .names = "{col}_{fn}"
    ),
    .groups = "drop"
  ) %>%
  as.data.frame()

# Display group-level statistics.

print(group_stats)


# ============================================================
# 10. BALLOON PLOT OF GROUP-LEVEL RSCU
# ============================================================

# Transpose the group mean table so that:
#
# Rows    = codons
# Columns = biological groups
#
# This format is used to visualize the average RSCU
# of each codon across the different groups.

cont_table_t <- as.data.frame(t(cont_table))

# Convert to a matrix for the balloon plot.

dt <- as.table(as.matrix(cont_table_t))

# Create balloon plot.
#
# This provides a visual comparison of codon usage
# across biological groups.

balloonplot(
  t(dt),
  main = "RSCU Codon Usage",
  xlab = "Groups",
  ylab = "Codons",
  label = FALSE,
  show.margins = FALSE
)


# ============================================================
# 11. CORRESPONDENCE ANALYSIS (CA)
# ============================================================

# PURPOSE:
# Examine associations between codons and samples.
#
# Unlike PCA, correspondence analysis represents the
# relationships between rows and columns of a nonnegative
# matrix.
#
# For this analysis:
#   Rows    = codons
#   Columns = individual samples
#
# The CA biplot displays both codons and samples.
#
# Samples positioned near one another have similar
# profiles in the CA representation.
#
# Codons positioned in similar directions can indicate
# similar patterns of relative usage across samples.
#
# NOTE:
# RSCU values are relative frequencies, not raw codon counts.
# CA here describes structure in the RSCU matrix rather
# than a contingency table of observed codon counts.

# Transpose the sample x codon matrix.
#
# Original: samples x codons
# Transposed: codons x samples

codon_matrix <- t(as.matrix(pca_matrix))

# Run correspondence analysis.

res.ca <- CA(
  codon_matrix,
  graph = FALSE
)


# ------------------------------------------------------------
# 11a. Eigenvalues and explained inertia
# ------------------------------------------------------------

# Extract eigenvalues and percentage of inertia explained
# by each CA dimension.

eig.val <- get_eigenvalue(res.ca)

print(eig.val)


# ------------------------------------------------------------
# 11b. Codon coordinates and contributions
# ------------------------------------------------------------

# Extract information for codon rows.

row <- get_ca_row(res.ca)

# Coordinates:
#   Position of each codon along each CA dimension.

head(row$coord)

# Cos2:
#   Quality of representation of each codon on each axis.
#   Higher values indicate better representation.

head(row$cos2)

# Contributions:
#   Contribution of each codon to the construction
#   of the CA dimensions.

head(row$contrib)


# ------------------------------------------------------------
# 11c. Visualize codon positions
# ------------------------------------------------------------

# Plot codon coordinates.

fviz_ca_row(
  res.ca,
  repel = TRUE
)

# Color codons by their contribution to the CA dimensions.

fviz_ca_row(
  res.ca,
  col.row = "contrib",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE
)

# Adjust label transparency according to contribution.

fviz_ca_row(
  res.ca,
  alpha.row = "contrib",
  repel = TRUE
)


# ------------------------------------------------------------
# 11d. Sample coordinates and contributions
# ------------------------------------------------------------

# Extract information for sample columns.

col <- get_ca_col(res.ca)

# Coordinates of individual samples.

head(col$coord)

# Quality of sample representation (cos2).

head(col$cos2)

# Sample contributions to the CA dimensions.

head(col$contrib)


# Plot individual samples.

fviz_ca_col(res.ca)

# Color samples according to their quality of representation.

fviz_ca_col(
  res.ca,
  col.col = "cos2",
  gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
  repel = TRUE
)

# Examine the quality of representation of samples
# on the first two CA dimensions.

fviz_cos2(
  res.ca,
  choice = "col",
  axes = 1:2
)


# ------------------------------------------------------------
# 11e. CA biplot
# ------------------------------------------------------------

# Create a biplot showing both:
#   - Codons (rows)
#   - Samples (columns)
#
# label = "row" displays codon labels but removes sample labels.
# This makes the figure easier to interpret when there are
# many individual samples.

RSCU_biplot <- fviz_ca_biplot(
  res.ca,
  repel = TRUE,
  label = "row"
)

# Display the biplot.

print(RSCU_biplot)

# Save as editable SVG.

ggsave(
  "RSCU_biplot.svg",
  plot = RSCU_biplot,
  device = svglite,
  width = 12,
  height = 10
)


# ------------------------------------------------------------
# 11f. Additional CA visualizations
# ------------------------------------------------------------

# Row principal map: emphasizes codon positions.

fviz_ca_biplot(
  res.ca,
  map = "rowprincipal",
  arrow = c(TRUE, TRUE),
  repel = TRUE
)

# Column Greenacre map: emphasizes sample positions.

fviz_ca_biplot(
  res.ca,
  map = "colgreen",
  arrow = c(TRUE, FALSE),
  repel = TRUE
)

# Identify variables associated with each CA dimension.

res.desc <- dimdesc(
  res.ca,
  axes = c(1, 2)
)

# Display the strongest associations with Dimension 1.

head(res.desc[[1]]$row, 4)


# ============================================================
# 12. INDIVIDUAL CODON COMPARISONS
# ============================================================

# PURPOSE:
# Examine differences in individual codon RSCU values
# among biological groups.
#
# AGA and CTT are used here as examples.
#
# Boxplots display:
#   - Median
#   - Interquartile range
#   - Distribution of individual sample values

# ------------------------------------------------------------
# 12a. AGA
# ------------------------------------------------------------

# Boxplot of AGA RSCU values.

plot_AGA <- ggplot(
  combined,
  aes(x = Group, y = AGA, fill = Group)
) +
  geom_boxplot() +
  theme_classic() +
  labs(
    title = "AGA Codon Usage",
    x = NULL,
    y = "RSCU"
  )

print(plot_AGA)


# Calculate descriptive statistics for AGA.

stats_AGA <- combined %>%
  group_by(Group) %>%
  summarise(
    n      = sum(!is.na(AGA)),
    mean   = mean(AGA, na.rm = TRUE),
    median = median(AGA, na.rm = TRUE),
    min    = min(AGA, na.rm = TRUE),
    max    = max(AGA, na.rm = TRUE),
    .groups = "drop"
  )

print(stats_AGA)


# ------------------------------------------------------------
# 12b. CTT
# ------------------------------------------------------------

# Boxplot of CTT RSCU values.

plot_CTT <- ggplot(
  combined,
  aes(x = Group, y = CTT, fill = Group)
) +
  geom_boxplot() +
  theme_classic() +
  labs(
    title = "CTT Codon Usage",
    x = NULL,
    y = "RSCU"
  )

print(plot_CTT)


# Calculate descriptive statistics for CTT.

stats_CTT <- combined %>%
  group_by(Group) %>%
  summarise(
    n      = sum(!is.na(CTT)),
    mean   = mean(CTT, na.rm = TRUE),
    median = median(CTT, na.rm = TRUE),
    min    = min(CTT, na.rm = TRUE),
    max    = max(CTT, na.rm = TRUE),
    .groups = "drop"
  )

print(stats_CTT)


# ============================================================
# 13. RSCU DISTRIBUTION HISTOGRAMS
# ============================================================

# PURPOSE:
# Examine the distribution of RSCU values across individual
# codons and biological groups.
#
# The data are reshaped from wide to long format:
#
# Original:
#   One row per sample, one column per codon.
#
# Long:
#   One row per sample-codon combination.
#
# Only selected groups are included in this figure.
#
# RSCU values between 0 and 3 are displayed to focus
# on the central distribution.
#
# The dashed vertical line at RSCU = 1 represents
# equal synonymous codon usage.

# Select groups to include in the histogram.

plot_groups <- c(
  "Phage_1",
  "Phage_2",
  "Chromid",
  "Plasmid",
  "Chromosome"
)


# ------------------------------------------------------------
# 13a. Reshape data
# ------------------------------------------------------------

# Convert codon columns into two columns:
#   Codon = codon identity
#   RSCU  = corresponding RSCU value

rscu_hist_long <- combined %>%
  filter(Group %in% plot_groups) %>%
  pivot_longer(
    cols = all_of(common_codons),
    names_to = "Codon",
    values_to = "RSCU"
  ) %>%
  filter(
    is.finite(RSCU),
    RSCU >= 0,
    RSCU <= 3
  ) %>%
  mutate(
    Group = factor(Group, levels = plot_groups)
  )


# ------------------------------------------------------------
# 13b. Generate histogram
# ------------------------------------------------------------

rscu_hist_plot <- ggplot(
  rscu_hist_long,
  aes(x = RSCU, fill = Group)
) +
  
  # Histogram bins represent the frequency of RSCU values.
  
  geom_histogram(
    binwidth = 0.10,
    boundary = 0,
    color = "white",
    alpha = 0.85
  ) +
  
  # Reference line representing equal synonymous codon usage.
  
  geom_vline(
    xintercept = 1,
    linetype = "dashed",
    color = "black",
    linewidth = 0.5
  ) +
  
  # Separate panels for each biological group.
  
  facet_wrap(
    ~ Group,
    ncol = 3,
    scales = "free_y"
  ) +
  
  # Display RSCU values from 0 to 3.
  
  scale_x_continuous(
    limits = c(0, 3)
  ) +
  
  # Apply consistent group colors.
  
  scale_fill_manual(
    values = color_map
  ) +
  
  theme_classic() +
  
  theme(
    legend.position = "none",
    strip.background = element_rect(
      fill = "grey95",
      color = "black"
    ),
    strip.text = element_text(
      face = "bold"
    )
  ) +
  
  labs(
    title = "RSCU Distributions by Genomic Source",
    subtitle = "Dashed line: RSCU = 1 (equal synonymous-codon use)",
    x = "RSCU",
    y = "Number of RSCU observations"
  )


# Display histogram.

print(rscu_hist_plot)

# Save histogram as SVG.

ggsave(
  "RSCU_codon_usage_histograms.svg",
  plot = rscu_hist_plot,
  device = svglite,
  width = 12,
  height = 8
)


# ============================================================
# 14. EXPORT SUMMARY TABLES
# ============================================================

# Save group-level mean, SD, and CV for all codons.

write.csv(
  group_stats,
  "RSCU_group_statistics.csv",
  row.names = FALSE
)

# Save group-level mean RSCU values.

write.csv(
  cont_table,
  "RSCU_group_means.csv"
)

# Save descriptive statistics for AGA and CTT.

write.csv(
  stats_AGA,
  "RSCU_AGA_statistics.csv",
  row.names = FALSE
)

write.csv(
  stats_CTT,
  "RSCU_CTT_statistics.csv",
  row.names = FALSE
)


# ============================================================
# END OF ANALYSIS
# ============================================================

cat("\nRSCU analysis completed successfully.\n")
cat("Samples analyzed:", nrow(combined), "\n")
cat("Codons analyzed:", ncol(pca_matrix), "\n")
cat("Groups analyzed:", length(unique(combined$Group)), "\n")