# =============================================================================
# CAI (Codon Adaptation Index) analysis of phage genomes
# -----------------------------------------------------------------------------
# Compares CAI between two phage clusters (Cluster1 vs Cluster2) using two
# reference sets:
#   (1) C/M/P-referenced CAI     : per-gene CAI for replicons C, M, and P
#   (2) Ribosomal-referenced CAI : per-gene CAI from the InterProScan-derived
#                                  reference set
#
# Outputs
#   Figures   : Mean_CAI_per_Genome.svg
#               CAI_functional_categories.svg
#   Tables    : cmp_highest_CAI_genes.csv, cmp_lowest_CAI_genes.csv,
#               ribo_highest_CAI_genes.csv, ribo_lowest_CAI_genes.csv,
#               genes_by_functional_category.csv,
#               genes_per_category_summary.csv
#   Sequences : CAI_outlier_genes_nt_annotated.fasta (input to the external
#               annotation step)
#
# Workflow
#   1. Load per-gene CAI values
#   2. Per-genome mean CAI, figure, and statistical tests
#   3. Attach Prokka annotations
#   4. Select top/bottom CAI genes per cluster and write CSVs
#   5. Write nucleotide FASTA of outlier genes for annotation
#   6. Assign functional categories, plot composition, write category tables
#   7. Save session info for reproducibility
# =============================================================================

library(dplyr)
library(readr)
library(tidyr)
library(stringr)
library(purrr)
library(ggplot2)
library(scales)
library(svglite)
library(Biostrings)

# -----------------------------------------------------------------------------
# Parameters
# -----------------------------------------------------------------------------

pct             <- 0.10   # Fraction of genes at each tail (0.10 = top/bottom 10%)
replicons       <- c("C", "M", "P")
prokka_dir      <- "prokka_phage_outputs"
annotation_file <- "final_annotation_summary_10.tsv"   # Pharokka/PHROG/BLAST output

# Consistent colour palette for functional categories
cat_palette_fine <- c(
  "Structural"             = "#6A3D9A",
  "Packaging"              = "#B15FD3",
  "Replication"            = "#1B9E77",
  "Recombination"          = "#66C2A5",
  "Integration"            = "#E7298A",
  "Regulation"             = "#F781BF",
  "Mobile Genetic Element" = "#66A61E",
  "Toxin-Antitoxin"        = "#E6AB02",
  "Lysis"                  = "#A6761D",
  "DNA Modification"       = "#8C564B",
  "Host/Metabolic"         = "#7570B3",
  "Hypothetical/Unknown"   = "#D9D9D9",
  "Other"                  = "#444444"
)

# -----------------------------------------------------------------------------
# Helper functions
# -----------------------------------------------------------------------------

# Read one EMBOSS-style file with lines of the form
#   "Sequence: <gene_id> CAI: <value>"
# and return a tibble with Gene_ID, CAI, and the source file name.
read_cai_file <- function(path) {
  lines <- readLines(path, warn = FALSE)
  lines <- lines[grepl("CAI:", lines)]
  
  tibble(
    Gene_ID = trimws(sub("Sequence:\\s*", "", sub("CAI:.*", "", lines))),
    CAI     = as.numeric(trimws(sub(".*CAI:\\s*", "", lines))),
    File    = basename(path)
  )
}

# Read all files in a directory matching `pattern` and stack them.
# Always returns the same columns, even when the directory is empty,
# so downstream code never fails on a missing column.
read_cai_dir <- function(dir, pattern = "\\.txt$") {
  files <- list.files(dir, pattern = pattern, full.names = TRUE)
  
  if (length(files) == 0) {
    warning("No files matching '", pattern, "' found in: ", dir)
    return(tibble(Gene_ID = character(), CAI = numeric(), File = character()))
  }
  
  map(files, read_cai_file) %>% bind_rows()
}

# Compute high/low CAI cutoffs within each cluster and flag each gene.
# Cutoffs are the (1 - pct) and pct quantiles of CAI within the cluster.
add_cluster_cutoffs <- function(df, pct) {
  df %>%
    group_by(Cluster) %>%
    mutate(
      high_cutoff = quantile(CAI, 1 - pct, na.rm = TRUE),
      low_cutoff  = quantile(CAI, pct,     na.rm = TRUE)
    ) %>%
    ungroup()
}

# -----------------------------------------------------------------------------
# 1. Load per-gene CAI values
# -----------------------------------------------------------------------------

stopifnot(
  "cluster1_CAI / cluster2_CAI directories not found in the working directory" =
    all(dir.exists(c("cluster1_CAI", "cluster2_CAI")))
)

# C/M/P-referenced CAI: one folder per cluster, with C, M, and P subfolders.
read_cmp_cluster <- function(cluster_dir, cluster_label) {
  map(replicons, function(rep) {
    read_cai_dir(file.path(cluster_dir, rep), pattern = "\\.txt$") %>%
      mutate(Replicon = rep)
  }) %>%
    bind_rows() %>%
    mutate(Cluster = cluster_label)
}

cmp_genes <- bind_rows(
  read_cmp_cluster("cluster1_CAI", "Cluster1"),
  read_cmp_cluster("cluster2_CAI", "Cluster2")
) %>%
  filter(!is.na(CAI))

# Ribosomal-referenced CAI: InterProScan-derived folders.
interpro_dirs <- list.files(
  "updated_cai_C_from_interproscan",
  pattern    = "^[Cc]luster[12]_interpro_C_cai$",
  full.names = TRUE
)

stopifnot(
  "No cluster1/cluster2 InterProScan directories found" = length(interpro_dirs) > 0
)

ribo_genes <- map(interpro_dirs, function(dir) {
  read_cai_dir(dir, pattern = "\\.cai$") %>%
    mutate(Cluster = paste0("Cluster", str_extract(basename(dir), "[12]")))
}) %>%
  bind_rows() %>%
  filter(!is.na(CAI))

cat("C/M/P genes loaded:", nrow(cmp_genes), "\n")
cat("Ribosomal genes loaded:", nrow(ribo_genes), "\n")

# -----------------------------------------------------------------------------
# 2. Per-genome mean CAI (C/M/P-referenced): figure and statistics
# -----------------------------------------------------------------------------

# One mean CAI value per genome file, per cluster and replicon.
cmp_genome_stats <- cmp_genes %>%
  group_by(File, Cluster, Replicon) %>%
  summarise(Mean_CAI = mean(CAI), .groups = "drop")

p_mean_cai <- ggplot(cmp_genome_stats,
                     aes(x = Cluster, y = Mean_CAI, fill = Cluster)) +
  geom_boxplot() +
  facet_wrap(~ Replicon) +
  theme_bw() +
  labs(title = "Mean CAI per Genome", x = NULL, y = "Mean CAI")

print(p_mean_cai)

ggsave("Mean_CAI_per_Genome.svg",
       plot   = p_mean_cai,
       device = svglite,
       width  = 8,
       height = 6,
       units  = "in")

# Statistics for C/M/P-referenced CAI, per replicon:
#   Shapiro-Wilk normality within each cluster, then Wilcoxon rank-sum
#   comparison between clusters.
for (rep in replicons) {
  sub_data <- cmp_genes %>% filter(Replicon == rep)
  cat("\n---", rep, "---\n")
  print(shapiro.test(sub_data$CAI[sub_data$Cluster == "Cluster1"]))
  print(shapiro.test(sub_data$CAI[sub_data$Cluster == "Cluster2"]))
  print(wilcox.test(CAI ~ Cluster, data = sub_data))
}

# Statistics for ribosomal-referenced CAI.
cat("\n--- Ribosomal-referenced ---\n")
print(shapiro.test(ribo_genes$CAI[ribo_genes$Cluster == "Cluster1"]))
print(shapiro.test(ribo_genes$CAI[ribo_genes$Cluster == "Cluster2"]))
print(wilcox.test(CAI ~ Cluster, data = ribo_genes))

# -----------------------------------------------------------------------------
# 3. Attach Prokka annotations (locus_tag -> gene, product)
# -----------------------------------------------------------------------------

annotation_lookup <- list.files(prokka_dir,
                                pattern   = "\\.tsv$",
                                recursive = TRUE,
                                full.names = TRUE) %>%
  map(~ read_tsv(.x, show_col_types = FALSE) %>%
        select(locus_tag, gene, product)) %>%
  bind_rows() %>%
  distinct(locus_tag, .keep_all = TRUE)

cmp_genes  <- cmp_genes  %>% left_join(annotation_lookup, by = c("Gene_ID" = "locus_tag"))
ribo_genes <- ribo_genes %>% left_join(annotation_lookup, by = c("Gene_ID" = "locus_tag"))

# -----------------------------------------------------------------------------
# 4. Select top and bottom CAI genes per cluster and write CSVs
# -----------------------------------------------------------------------------

cmp_ranked  <- add_cluster_cutoffs(cmp_genes,  pct)
ribo_ranked <- add_cluster_cutoffs(ribo_genes, pct)

cmp_highest  <- cmp_ranked  %>% filter(CAI >= high_cutoff)
cmp_lowest   <- cmp_ranked  %>% filter(CAI <= low_cutoff)
ribo_highest <- ribo_ranked %>% filter(CAI >= high_cutoff)
ribo_lowest  <- ribo_ranked %>% filter(CAI <= low_cutoff)

write_csv(cmp_highest,  "cmp_highest_CAI_genes.csv")
write_csv(cmp_lowest,   "cmp_lowest_CAI_genes.csv")
write_csv(ribo_highest, "ribo_highest_CAI_genes.csv")
write_csv(ribo_lowest,  "ribo_lowest_CAI_genes.csv")

# -----------------------------------------------------------------------------
# 5. Nucleotide FASTA of outlier genes for functional annotation
# -----------------------------------------------------------------------------

# Combine all outlier sets, labelling each row with its origin.
genes_of_interest <- bind_rows(
  cmp_highest  %>% mutate(Source = paste0("cmp_",  Cluster, "_highest")),
  cmp_lowest   %>% mutate(Source = paste0("cmp_",  Cluster, "_lowest")),
  ribo_highest %>% mutate(Source = paste0("ribo_", Cluster, "_highest")),
  ribo_lowest  %>% mutate(Source = paste0("ribo_", Cluster, "_lowest"))
)

# Load all Prokka nucleotide gene files and strip headers to the locus tag.
ffn_files <- list.files(prokka_dir,
                        pattern   = "\\.ffn$",
                        recursive = TRUE,
                        full.names = TRUE)

all_ffn <- do.call(c, lapply(ffn_files, readDNAStringSet))
names(all_ffn) <- sub(" .*", "", names(all_ffn))

# Match each requested gene to its sequence. A gene can appear in more than
# one outlier set, so each row becomes its own FASTA entry with a unique header.
idx   <- match(genes_of_interest$Gene_ID, names(all_ffn))
found <- !is.na(idx)

outlier_seqs <- all_ffn[idx[found]]
names(outlier_seqs) <- paste0(
  genes_of_interest$Gene_ID[found], "|",
  genes_of_interest$Source[found],  "|",
  "CAI=", round(genes_of_interest$CAI[found], 3)
)

writeXStringSet(outlier_seqs, filepath = "CAI_outlier_genes_nt_annotated.fasta")

cat("Genes requested:", nrow(genes_of_interest), "\n")
cat("Sequences found:", sum(found), "\n")
cat("Missing:", paste(genes_of_interest$Gene_ID[!found], collapse = ", "), "\n")

# -----------------------------------------------------------------------------
# 6. Functional categories
# -----------------------------------------------------------------------------
# Requires the annotation summary produced by Pharokka/PHROG/BLAST.
# Parsing and category rules are carried over unchanged from the original
# working scripts.

annotation_cat <- read_tsv(annotation_file, show_col_types = FALSE) %>%
  mutate(target_clean = str_remove(target, "_\\d+$")) %>%
  separate(target_clean, into = c("Gene_ID", "Source", "CAI"), sep = "\\|") %>%
  mutate(
    CAI       = as.numeric(str_remove(CAI, "CAI=")),
    CAI_Group = ifelse(str_detect(Source, "highest"), "High CAI", "Low CAI"),
    Reference = ifelse(str_detect(Source, "^cmp"),
                       "C/M/P-referenced", "Ribosomal-referenced"),
    Cluster   = str_extract(Source, "(?<=_)[Cc]luster[12]"),
    Category  = case_when(
      str_detect(product, regex("tail|capsid|head|spike|portal|baseplate|fiber|virion structural|tape measure",
                                ignore_case = TRUE)) ~ "Structural",
      str_detect(product, regex("terminase",
                                ignore_case = TRUE)) ~ "Packaging",
      str_detect(product, regex("primase|polymerase|helicase|ligase|nuclease|replicat|single strand DNA binding",
                                ignore_case = TRUE)) ~ "Replication",
      str_detect(product, regex("ssDNA annealing|recombinase|recombination",
                                ignore_case = TRUE)) ~ "Recombination",
      str_detect(product, regex("integrase",
                                ignore_case = TRUE)) ~ "Integration",
      str_detect(product, regex("lysin|holin|lysozyme|spanin",
                                ignore_case = TRUE)) ~ "Lysis",
      str_detect(product, regex("methyltransferase",
                                ignore_case = TRUE)) ~ "DNA Modification",
      str_detect(product, regex("transcriptional regulator|transcriptional repressor|helix-turn-helix|HTH",
                                ignore_case = TRUE)) ~ "Regulation",
      str_detect(product, regex("toxin.antitoxin|antitoxin|RelE|RelBE|YafQ|DinJ|HicB",
                                ignore_case = TRUE)) ~ "Toxin-Antitoxin",
      str_detect(product, regex("transposase|IS110|insertion sequence",
                                ignore_case = TRUE)) ~ "Mobile Genetic Element",
      str_detect(product, regex("galactosidase|metabolic|auxiliary|nitrogen regulation|hydroxylase",
                                ignore_case = TRUE)) ~ "Host/Metabolic",
      str_detect(product, regex("hypothetical|uncharacterized|protein of unknown function|domain of unknown function|DUF",
                                ignore_case = TRUE)) | is.na(product) ~ "Hypothetical/Unknown",
      TRUE ~ "Other"
    )
  )

# Check which products fell through to "Other"
annotation_cat %>%
  filter(Category == "Other") %>%
  count(product, sort = TRUE)

# Composition plot. Hypothetical/Unknown genes are excluded so the plot
# reflects known function.
annotation_cat_known <- annotation_cat %>%
  filter(Category != "Hypothetical/Unknown")

p_categories <- ggplot(annotation_cat_known,
                       aes(x = CAI_Group, fill = Category)) +
  geom_bar(position = "fill", color = "grey20", linewidth = 0.15) +
  facet_grid(Cluster ~ Reference) +
  scale_y_continuous(labels = percent) +
  scale_fill_manual(values = cat_palette_fine) +
  theme_bw() +
  labs(
    title    = "Functional category composition of high vs low CAI genes",
    subtitle = paste0("Top/bottom ", pct * 100, "% CAI genes per cluster (annotated genes only)"),
    x        = NULL,
    y        = "Proportion of genes",
    fill     = "Category"
  )

print(p_categories)

ggsave("CAI_functional_categories.svg",
       plot   = p_categories,
       device = svglite,
       width  = 9,
       height = 6,
       units  = "in")

# Gene-level table, sorted by category and then by CAI (highest first).
gene_category_table <- annotation_cat %>%
  select(Gene_ID, product, Category, Cluster, Reference, CAI_Group, CAI) %>%
  arrange(Category, desc(CAI))

write_csv(gene_category_table, "genes_by_functional_category.csv")

# One row per category, with all member genes listed.
genes_per_category <- annotation_cat %>%
  group_by(Category) %>%
  summarise(
    n_genes  = n(),
    Gene_IDs = paste(Gene_ID, collapse = ", "),
    .groups  = "drop"
  ) %>%
  arrange(desc(n_genes))

write_csv(genes_per_category, "genes_per_category_summary.csv")

# -----------------------------------------------------------------------------
# 7. Session info for reproducibility
# -----------------------------------------------------------------------------

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")