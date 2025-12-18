W1118_r1_a_mRNA = read.delim("annotated_modkit_outfiles/W1118_induro_mRNA_2025_102.sorted_a.tsv") %>% dplyr::filter(score >= 5, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
W1118_r1_m_mRNA = read.delim("annotated_modkit_outfiles/W1118_induro_mRNA_2025_102.sorted_m.tsv") %>% dplyr::filter(score >= 5, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
W1118_r1_pseu_mRNA = read.delim("annotated_modkit_outfiles/W1118_induro_mRNA_2025_102.sorted_17802.tsv") %>% dplyr::filter(score >= 5, percent_modified > 15, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)

W1118_r2_a_mRNA = read.delim("annotated_modkit_outfiles/W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.sorted_a.tsv") %>% dplyr::filter(score >= 5, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
W1118_r2_m_mRNA = read.delim("annotated_modkit_outfiles/W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.sorted_m.tsv") %>% dplyr::filter(score >= 5, Location != "N/A", Location != " ", Location != "intron" ) %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
W1118_r2_pseu_mRNA = read.delim("annotated_modkit_outfiles/W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.sorted_17802.tsv") %>% dplyr::filter(score >= 5, percent_modified > 15, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)

W1118_r3_a_mRNA = read.delim("annotated_modkit_outfiles/W1118_rep2_ref.sorted_a.tsv") %>% dplyr::filter(score >= 5, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
W1118_r3_m_mRNA = read.delim("annotated_modkit_outfiles/W1118_rep2_ref.sorted_m.tsv") %>% dplyr::filter(score >= 5, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
W1118_r3_pseu_mRNA = read.delim("annotated_modkit_outfiles/W1118_rep2_ref.sorted_17802.tsv") %>% dplyr::filter(score >= 5, percent_modified > 15, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)

colnames(W1118_r1_a_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                               "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(W1118_r1_m_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                               "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(W1118_r1_pseu_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                                  "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")

colnames(W1118_r2_a_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                               "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(W1118_r2_m_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                               "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(W1118_r2_pseu_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                                  "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")


colnames(W1118_r3_a_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                               "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(W1118_r3_m_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                               "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(W1118_r3_pseu_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                                  "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")


FR1113N_r1_a_mRNA = read.delim("annotated_modkit_outfiles/FR113N_05_11_2025_rep2_ref.sorted_a.tsv") %>% dplyr::filter(score >= 5, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
FR1113N_r1_m_mRNA = read.delim("annotated_modkit_outfiles/FR113N_05_11_2025_rep2_ref.sorted_m.tsv") %>% dplyr::filter(score >= 5, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
FR1113N_r1_pseu_mRNA = read.delim("annotated_modkit_outfiles/FR113N_05_11_2025_rep2_ref.sorted_17802.tsv") %>% dplyr::filter(score >= 5, percent_modified > 15, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)

FR1113N_r2_a_mRNA = read.delim("annotated_modkit_outfiles/FR113N_05_11_2025_rep3_ref.sorted_a.tsv") %>% dplyr::filter(score >= 5, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
FR1113N_r2_m_mRNA = read.delim("annotated_modkit_outfiles/FR113N_05_11_2025_rep3_ref.sorted_m.tsv") %>% dplyr::filter(score >= 5, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
FR1113N_r2_pseu_mRNA = read.delim("annotated_modkit_outfiles/FR113N_05_11_2025_rep3_ref.sorted_17802.tsv") %>% dplyr::filter(score >= 5, percent_modified > 15, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)

FR1113N_r3_a_mRNA = read.delim("annotated_modkit_outfiles/FR113N_Oct_4_mRNA_induro.sorted_a.tsv") %>% dplyr::filter(score >= 5, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
FR1113N_r3_m_mRNA = read.delim("annotated_modkit_outfiles/FR113N_Oct_4_mRNA_induro.sorted_m.tsv") %>% dplyr::filter(score >= 5, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
FR1113N_r3_pseu_mRNA = read.delim("annotated_modkit_outfiles/FR113N_Oct_4_mRNA_induro.sorted_17802.tsv") %>% dplyr::filter(score >= 5, percent_modified > 15, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)

colnames(FR1113N_r1_a_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                                 "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(FR1113N_r1_m_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                                 "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(FR1113N_r1_pseu_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                                    "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")

colnames(FR1113N_r2_a_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                                 "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(FR1113N_r2_m_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                                 "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(FR1113N_r2_pseu_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                                    "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")


colnames(FR1113N_r3_a_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                                 "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(FR1113N_r3_m_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                                 "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(FR1113N_r3_pseu_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand",
                                    "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")



Mettl3_r1_a_mRNA = read.delim("annotated_modkit_outfiles/Mettl3_KO_mRNA_RNA004.sorted_a.tsv") %>% dplyr::filter(score >= 10, percent_modified > 10, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
Mettl3_r1_m_mRNA = read.delim("annotated_modkit_outfiles/Mettl3_KO_mRNA_RNA004.sorted_m.tsv") %>% dplyr::filter(score >= 10, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
Mettl3_r1_pseu_mRNA = read.delim("annotated_modkit_outfiles/Mettl3_KO_mRNA_RNA004.sorted_17802.tsv") %>% dplyr::filter(score >= 10, percent_modified > 20, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)

Mettl3_r2_a_mRNA = read.delim("annotated_modkit_outfiles/Mettl3_KO_mRNA_rep2.sorted_a.tsv") %>% dplyr::filter(score >= 10, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
Mettl3_r2_m_mRNA = read.delim("annotated_modkit_outfiles/Mettl3_KO_mRNA_rep2.sorted_m.tsv") %>% dplyr::filter(score >= 10, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)
Mettl3_r2_pseu_mRNA = read.delim("annotated_modkit_outfiles/Mettl3_KO_mRNA_rep2.sorted_17802.tsv") %>% dplyr::filter(score >= 10, Location != "N/A", Location != " ", Location != "intron") %>% dplyr::select(-count_diff, -count_nocall, -valid_coverage)



colnames(Mettl3_r1_a_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand", "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(Mettl3_r1_m_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand", "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(Mettl3_r1_pseu_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand", "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")

colnames(Mettl3_r2_a_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand", "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(Mettl3_r2_m_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand", "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")
colnames(Mettl3_r2_pseu_mRNA) <- c("chrom", "start", "end", "mod", "score", "strand", "percent_modified", "Gene_ID", "Transcript_IDs", "Transcript_Biotype", "Location")

peak_files <- list(
  "w1118_Rep1_m6A" = W1118_r1_a_mRNA,
  "w1118_Rep1_m5C" = W1118_r1_m_mRNA,
  'w1118_Rep1_Ψ' = W1118_r1_pseu_mRNA,
  "w1118_Rep2_m6A" = W1118_r2_a_mRNA,
  'w1118_Rep2_m5C' = W1118_r2_m_mRNA,
  'w1118_Rep2_Ψ' = W1118_r2_pseu_mRNA,
  'w1118_Rep3_m6A' = W1118_r3_a_mRNA,
  'w1118_Rep3_m5C' = W1118_r3_m_mRNA,
  'w1118_Rep3_Ψ' = W1118_r3_pseu_mRNA,
  'FR113_Rep1_m6A' = FR1113N_r1_a_mRNA,
  'FR113_Rep1_m5C' = FR1113N_r1_m_mRNA,
  "FR113_Rep1_Ψ" = FR1113N_r1_pseu_mRNA,
  "FR113_Rep2_m6A" = FR1113N_r2_a_mRNA,
  "FR113_Rep2_m5C" = FR1113N_r2_m_mRNA, 
  "FR113_Rep2_Ψ" = FR1113N_r2_pseu_mRNA,
  "FR113_Rep3_m6A" = FR1113N_r3_a_mRNA,
  "FR113_Rep3_m5C" = FR1113N_r3_m_mRNA, 
  "FR113_Rep3_Ψ" = FR1113N_r3_pseu_mRNA,
  "Mettl3_Rep1_m6A" = Mettl3_r1_a_mRNA,
  "Mettl3_Rep1_m5C" = Mettl3_r1_m_mRNA, 
  "Mettl3_Rep1_Ψ" = Mettl3_r1_pseu_mRNA,
  "Mettl3_Rep2_m6A" = Mettl3_r2_a_mRNA,
  "Mettl3_Rep2_m5C" = Mettl3_r2_m_mRNA, 
  "Mettl3_Rep2_Ψ" = Mettl3_r2_pseu_mRNA)


peak_counts <- sapply(peak_files, nrow)
peak_counts_df <- data.frame(Sample = factor(names(peak_counts), levels = names(peak_counts)), Peaks = peak_counts)

# Plot the bar plot
mod.dist <- ggplot(peak_counts_df, aes(x = Sample, y = Peaks, fill = Sample)) +
  geom_bar(stat = "identity", show.legend = FALSE) +
  geom_text(aes(label = Peaks), vjust = -0.3, color = "black", size = 3.5) +
  scale_fill_manual(values = c('W1118_Rep1_m6A' = "#1b9e77", 'W1118_Rep1_m5C' = "#1b9e77", 'W1118_Rep1_Ψ' = "#1b9e77", 'W1118_Rep2_m6A' = "#1b9e77", 'W1118_Rep2_m5C' = "#1b9e77", 'W1118_Rep2_Ψ' = "#1b9e77", 'W1118_Rep3_m6A' = "#1b9e77", 'W1118_Rep3_m5C' = "#1b9e77", 'W1118_Rep3_Ψ' = "#1b9e77", "FR113_Rep1_m6A" = "#d95f02", "FR113_Rep1_m5C" = "#d95f02", "FR113_Rep1_Ψ" = "#d95f02", "FR113_Rep2_m6A" = "#d95f02", "FR113_Rep2_m5C" = "#d95f02", "FR113_Rep2_Ψ" = "#d95f02","FR113_Rep3_m6A" = "#d95f02","FR113_Rep3_m5C" = "#d95f02", "FR113_Rep3_Ψ" = "#d95f02","Mettl3_Rep1_m6A" = "#7570b3", "Mettl3_Rep1_m5C" = "#7570b3", "Mettl3_Rep1_Ψ" = "#7570b3", "Mettl3_Rep2_m6A" = "#7570b3","Mettl3_Rep2_m5C" = "#7570b3", "Mettl3_Rep2_Ψ" = "#7570b3")) +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  labs(title = "Number of Peaks per Sample(RNA004)", x = "Sample", y = "Number of mod identified") + theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + theme_classic(base_size = 22, base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "right", size = 12,
        legend.text=element_text(size=12, family = "Times New Roman")) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))


## 1. Define condition names (drop the replicate part)
condition_names <- sub("_Rep[0-9]+", "", names(peak_files))
# e.g. "W1118_Rep1_m6A" -> "W1118_m6A"

## 2. Split the list by condition
peaks_by_condition <- split(peak_files, condition_names)

## 3. For each condition, merge reps and remove duplicate sites

unique_counts <- map_int(peaks_by_condition, ~ {
  df_combined <- bind_rows(.x)
  df_unique <- distinct(df_combined, chrom, start, end, strand, .keep_all = TRUE)
  
  nrow(df_unique)
})

peak_counts_df <- data.frame(
  Condition = names(unique_counts),
  Peaks = as.integer(unique_counts))

# Split Condition into Group (genotype) and Mod (modification)
peak_counts_df <- peak_counts_df %>%
  mutate(Group = sub("_.*", "", Condition),          # W1118 / FR113 / Mettl3
         Mod = sub("^[^_]+_", "", Condition))       # m6A / m5C / Ψ)

# Order conditions if you want a specific x-axis order
peak_counts_df$Condition <- factor(
  peak_counts_df$Condition,
  levels = c("w1118_m6A", "w1118_m5C", "w1118_Ψ", "FR113_m6A", "FR113_m5C", "FR113_Ψ", "Mettl3_m6A","Mettl3_m5C","Mettl3_Ψ"))

## 5. Plot: one bar per condition (reps merged, duplicates collapsed)
mod.dist <- ggplot(peak_counts_df, aes(x = Condition, y = Peaks, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Peaks), vjust = -0.3, color = "black", size = 3.5) +
  scale_fill_manual(values = c("w1118"  = "#1b9e77", "FR113"  = "#d95f02", "Mettl3" = "#7570b3")) +
  labs(title = "Number of unique modified sites per condition", x = "Condition", y = "Number of unique sites") +
  theme_classic(base_size = 22, base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right",
        legend.text = element_text(size = 12, family = "Times New Roman"))

ggsave("fig2/percent_mod_ind_dis.pdf", plot = mod.dist, width = 9, height = 9, units = "in", device = cairo_pdf, bg = "white")

# Get number of genes per genotype and modification:

## 1. Define condition names (drop the replicate part)
condition_names <- sub("_Rep[0-9]+", "", names(peak_files))
# e.g. "W1118_Rep1_m6A" -> "W1118_m6A"

## 2. Split the list by condition (all reps for same condition grouped)
peaks_by_condition <- split(peak_files, condition_names)

## 3. For each condition, merge reps and count unique genes
##    >>> change `gene_id` to the actual column name for genes <<<
gene_counts <- map_int(peaks_by_condition, ~ {
  df_combined <- bind_rows(.x)
  n_distinct(df_combined$Gene_ID)
})

## 4. Build plotting data.frame
peak_gene_counts_df <- data.frame(
  Condition = names(gene_counts),
  Genes = as.integer(gene_counts))

# Split Condition into Group (genotype) and Mod (modification)
peak_gene_counts_df <- peak_gene_counts_df %>%
  mutate(Group = sub("_.*", "", Condition),          # W1118 / FR113 / Mettl3
         Mod = sub("^[^_]+_", "", Condition))    # m6A / m5C / Ψ


# Optional: order conditions on x-axis
peak_gene_counts_df$Condition <- factor(
  peak_gene_counts_df$Condition,
  levels = c("w1118_m6A", "w1118_m5C", "w1118_Ψ", "FR113_m6A", "FR113_m5C", "FR113_Ψ","Mettl3_m6A","Mettl3_m5C","Mettl3_Ψ"))

## 5. Plot: number of genes with ≥1 modified site per condition
mod.genes <- ggplot(peak_gene_counts_df, aes(x = Condition, y = Genes, fill = Group)) +
  geom_bar(stat = "identity") +
  geom_text(aes(label = Genes), vjust = -0.3, color = "black", size = 3.5) +
  scale_fill_manual(values = c("w1118"  = "#1b9e77","FR113"  = "#d95f02", "Mettl3" = "#7570b3")) +
  labs(title = "Number of genes with modifications per condition (RNA004)",
       x = "Condition", y = "Number of genes") +
  theme_classic(base_size = 22, base_family = "Times New Roman") +
  theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = "right",
        legend.text = element_text(size = 12, family = "Times New Roman"))
ggsave("fig2/num_gene_ind_dis.pdf", plot = mod.genes, width = 9, height = 9, units = "in", device = cairo_pdf, bg = "white")
make peak_counts

# Packages
library(dplyr)
library(purrr)
library(tidyr)
library(ggplot2)
library(ggrepel)
library(tibble)

## Helper: resolve object names robustly (handles FR113N vs FR1113N)
resolve_df <- function(candidates) {
  for (nm in candidates) {
    if (exists(nm, inherits = TRUE)) return(get(nm, inherits = TRUE))
  }
  stop("None of these objects exist: ", paste(candidates, collapse = ", "))
}

sample_list <- list(
  # W1118
  w1118_rep1_m6A = resolve_df(c("W1118_r1_a_mRNA")),
  w1118_rep1_m5C = resolve_df(c("W1118_r1_m_mRNA")),
  w1118_rep1_pseu = resolve_df(c("W1118_r1_pseu_mRNA")),
  w1118_rep2_m6A = resolve_df(c("W1118_r2_a_mRNA")),
  w1118_rep2_m5C = resolve_df(c("W1118_r2_m_mRNA")),
  w1118_rep2_pseu = resolve_df(c("W1118_r2_pseu_mRNA")),
  w1118_rep3_m6A = resolve_df(c("W1118_r3_a_mRNA")),
  w1118_rep3_m5C = resolve_df(c("W1118_r3_m_mRNA")),
  w1118_rep3_pseu = resolve_df(c("W1118_r3_pseu_mRNA")),
  
  # FR113N  (also accept FR1113N_* if that's how objects are named)
  FR113N_rep1_m6A = resolve_df(c("FR113N_r1_a_mRNA", "FR1113N_r1_a_mRNA")),
  FR113N_rep1_m5C = resolve_df(c("FR113N_r1_m_mRNA", "FR1113N_r1_m_mRNA")),
  FR113N_rep1_pseu = resolve_df(c("FR113N_r1_pseu_mRNA", "FR1113N_r1_pseu_mRNA")),
  FR113N_rep2_m6A = resolve_df(c("FR113N_r2_a_mRNA", "FR1113N_r2_a_mRNA")),
  FR113N_rep2_m5C = resolve_df(c("FR113N_r2_m_mRNA", "FR1113N_r2_m_mRNA")),
  FR113N_rep2_pseu = resolve_df(c("FR113N_r2_pseu_mRNA", "FR1113N_r2_pseu_mRNA")),
  FR113N_rep3_m6A = resolve_df(c("FR113N_r3_a_mRNA", "FR1113N_r3_a_mRNA")),
  FR113N_rep3_m5C = resolve_df(c("FR113N_r3_m_mRNA", "FR1113N_r3_m_mRNA")),
  FR113N_rep3_pseu = resolve_df(c("FR113N_r3_pseu_mRNA", "FR1113N_r3_pseu_mRNA")),
  
  # Mettl3-KO  (objects commonly named Mettl3_*; keep that but label as Mettl3-KO)
  `Mettl3-KO_rep1_m6A` = resolve_df(c("Mettl3_r1_a_mRNA")),
  `Mettl3-KO_rep1_m5C` = resolve_df(c("Mettl3_r1_m_mRNA")),
  `Mettl3-KO_rep1_pseu` = resolve_df(c("Mettl3_r1_pseu_mRNA")),
  `Mettl3-KO_rep2_m6A` = resolve_df(c("Mettl3_r2_a_mRNA")),
  `Mettl3-KO_rep2_m5C` = resolve_df(c("Mettl3_r2_m_mRNA")),
  `Mettl3-KO_rep2_pseu` = resolve_df(c("Mettl3_r2_pseu_mRNA")))

## 2) Build a site key and extract the score column
# Assumes input dfs have columns: Gene_ID, percent_modified
prep_one <- function(df, sample_name) {
  df %>%
    mutate(site = as.character(Gene_ID)) %>%
    distinct(site, .keep_all = TRUE) %>%
    dplyr::select(site, percent_modified) %>%
    rename(!!sample_name := percent_modified)
}

# Intersection of sites across all samples (change to full_join for union)
wide_scores <- imap(sample_list, prep_one) %>%
  reduce(inner_join, by = "site")

## 3) Matrix: rows = sites, cols = samples
mat_sites_x_samples <- wide_scores %>%
  arrange(site) %>%
  column_to_rownames("site") %>%
  as.matrix()

## 4) PCA on samples (transpose so samples are observations)
pca <- prcomp(t(mat_sites_x_samples), center = TRUE, scale. = TRUE)
var_expl <- (pca$sdev^2) / sum(pca$sdev^2)
percentVar <- round(100 * var_expl, 1)
pc1_lab <- paste0("PC1 (", percentVar[1], "%)")
pc2_lab <- paste0("PC2 (", percentVar[2], "%)")

## 5) Sample metadata parsed from names
pca_data <- as.data.frame(pca$x) %>%
  rownames_to_column("sample") %>%
  separate(sample, into = c("Strain", "Rep", "Mod"), sep = "_", remove = FALSE, fill = "right") %>%
  mutate(
    Rep = factor(Rep),
    Mod = factor(Mod, levels = c("m6A", "m5C", "pseu")),
    condition = factor(Strain, levels = c("w1118", "FR113N", "Mettl3-KO")))

custom_colors <- c("w1118" = "#1b9e77","FR113N" = "#d95f02", "Mettl3-KO" = "#7570b3")

percent_mod <- ggplot(pca_data, aes(PC1, PC2, color = condition, fill = condition)) +
  stat_ellipse(type = "euclid", alpha = 0.9, colour = "black") +
  geom_point(size = 4, shape = 21, stroke = 0.6) +
  labs(title = "PCA of samples (percent_modified)", x = paste0("PC1: ", percentVar[1], "% variance"), y = paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  scale_fill_manual(values = custom_colors, drop = FALSE) +
  scale_color_manual(values = custom_colors, drop = FALSE) +
  theme_classic(base_size = 12) +
  theme(panel.border = element_rect(fill = NA, color = "grey70"), plot.title = element_text(hjust = 0.5, face = "bold")) + theme_classic() +
  guides(color = guide_legend(title = "Condition"),
         fill  = guide_legend(title = "Condition")) +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "right",
        legend.text = element_text(size = 12, family = "Times New Roman"))
ggsave("fig2/percent_mod_PCA.pdf", plot = percent_mod, width = 9, height = 9, units = "in", device = cairo_pdf, bg = "white")


#Upset plot of shared genes across genotypes:

library(dplyr)
library(purrr)
library(tibble)
library(ComplexUpset)
library(ggplot2)

custom_colors <- c("w1118" = "#1b9e77", "FR113N" = "#d95f02", "Mettl3-KO" = "#7570b3")

# helper: genes per modification × genotype
get_genes_mod_genotype <- function(mod_pattern, geno_pattern) {
  peak_files[grepl(geno_pattern, names(peak_files)) &
               grepl(mod_pattern,  names(peak_files))] %>%
    map(~ .x$Gene_ID) %>%
    unlist() %>%
    unique()
}

make_upset_for_mod <- function(mod_pattern, mod_label = mod_pattern) {
  # gene sets for this modification
  genes_FR113 <- get_genes_mod_genotype(mod_pattern, "^FR113")
  genes_w1118 <- get_genes_mod_genotype(mod_pattern, "^w1118")
  genes_Mettl3 <- get_genes_mod_genotype(mod_pattern, "^Mettl3")
  
  all_genes <- Reduce(union, list(genes_FR113, genes_w1118, genes_Mettl3))
  
  df <- tibble(Gene_ID = all_genes,
               FR113N = Gene_ID %in% genes_FR113,
               w1118 = Gene_ID %in% genes_w1118,
               `Mettl3-KO` = Gene_ID %in% genes_Mettl3)
  
  message("Counts for ", mod_label, " genes:")
  print(
    df %>%
      summarise(
        FR113N = sum(FR113N),
        w1118 = sum(w1118),
        Mettl3_KO = sum(`Mettl3-KO`)))
  
  upset(df, intersect = c("FR113N", "w1118", "Mettl3-KO"), name = paste(mod_label, "genes"),
        base_annotations = list("Intersection size" = intersection_size(text = list(size = 3)))) +
    # Try to map your colors onto set-related geoms
    scale_fill_manual(values = custom_colors, breaks = names(custom_colors)) +
    scale_color_manual(values = custom_colors, breaks = names(custom_colors)) +
    labs(title = paste("Shared", mod_label, "genes across genotypes")) +
    theme(text = element_text(family = "Times New Roman", size = 12),
          plot.title = element_text(face = "bold", hjust = 0.5),
          axis.text.x = element_text(size = 10),
          axis.text.y = element_text(size = 10))
}

# Make the three “nice” plots
p_m6A_genes <- make_upset_for_mod("m6A", "m6A")
p_m5C_genes <- make_upset_for_mod("m5C", "m5C")
p_psi_genes <- make_upset_for_mod("Ψ",   "Ψ")

ggsave("fig2/upset_m6A_genes_FR113N_W1118_Mettl3KO.pdf", p_m6A_genes,
       width = 6, height = 4, units = "in", device = cairo_pdf, bg = "white")
ggsave("fig2/upset_m5C_genes_FR113N_W1118_Mettl3KO.pdf", p_m5C_genes,
       width = 6, height = 4, units = "in", device = cairo_pdf, bg = "white")
ggsave("fig2/upset_psi_genes_FR113N_W1118_Mettl3KO.pdf", p_psi_genes,
       width = 6, height = 4, units = "in", device = cairo_pdf, bg = "white")




Fig2

library(dplyr)
library(readr)

cols_to_drop <- c(
  "read_id","forward_read_position","ref_position","chrom","mod_strand",
  "ref_strand","ref_mod_strand","fw_soft_clipped_start","fw_soft_clipped_end",
  "alignment_start","alignment_end","inferred")

load_modkit <- function(path) {
  readr::read_tsv(path, show_col_types = FALSE) %>%
    dplyr::select(-dplyr::any_of(cols_to_drop)) %>%
    dplyr::mutate(call_prob = suppressWarnings(as.numeric(call_prob))) %>%
    dplyr::filter(!is.na(call_prob) & call_prob > 0)
}

W1118_r1_a_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/W1118_induro_mRNA_2025_102.modkit.mod.a.tsv") 
W1118_r1_m_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/W1118_induro_mRNA_2025_102.modkit.mod.m.tsv")
W1118_r1_pseu_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/W1118_induro_mRNA_2025_102.modkit.mod.17802.tsv")

W1118_r2_a_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.modkit.mod.a.tsv") 
W1118_r2_m_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.modkit.mod.m.tsv") 
W1118_r2_pseu_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.modkit.mod.17802.tsv") 

W1118_r3_a_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/w1118_ind_rev_rep1_mRNA_m6A_m5C_psi.modkit.mod.a.tsv") 
W1118_r3_m_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/w1118_ind_rev_rep1_mRNA_m6A_m5C_psi.modkit.mod.m.tsv") 
W1118_r3_pseu_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/w1118_ind_rev_rep1_mRNA_m6A_m5C_psi.modkit.mod.17802.tsv")


FR1113N_r1_a_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/FR113N_05_11_2025_rep2_ref.modkit.mod.a.tsv") 
FR1113N_r1_m_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/FR113N_05_11_2025_rep2_ref.modkit.mod.m.tsv") 
FR1113N_r1_pseu_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/FR113N_05_11_2025_rep2_ref.modkit.mod.17802.tsv") 

FR1113N_r2_a_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/FR113N_05_11_2025_rep3_ref.modkit.mod.a.tsv")
FR1113N_r2_m_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/FR113N_05_11_2025_rep3_ref.modkit.mod.m.tsv") 
FR1113N_r2_pseu_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/FR113N_05_11_2025_rep3_ref.modkit.mod.17802.tsv")

FR1113N_r3_a_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/FR113N_INDURO_mRNA_09_23.modkit.mod.a.tsv") 
FR1113N_r3_m_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/FR113N_INDURO_mRNA_09_23.modkit.mod.m.tsv") 
FR1113N_r3_pseu_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/FR113N_INDURO_mRNA_09_23.modkit.mod.17802.tsv")


Mettl3_r1_a_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/Mettl3_KO_mRNA_RNA004.modkit.mod.a.tsv") 
Mettl3_r1_m_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/Mettl3_KO_mRNA_RNA004.modkit.mod.m.tsv") 
Mettl3_r1_pseu_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/Mettl3_KO_mRNA_RNA004.modkit.mod.17802.tsv")

Mettl3_r2_a_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/Mettl3_KO_mRNA_rep2.modkit.mod.a.tsv") 
Mettl3_r2_m_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/Mettl3_KO_mRNA_rep2.modkit.mod.m.tsv") 
Mettl3_r2_pseu_mRNA_kmer = load_modkit("/mnt/raid1/george/direct_RNA/modkit_pileups/modkit_extra_calls/Mettl3_KO_mRNA_rep2.modkit.mod.17802.tsv")


library(dplyr)
library(ggplot2)
library(scales)

W1118_r1_a_mRNA_kmer  <- W1118_r1_a_mRNA_kmer  %>% mutate(call_code = as.character(call_code))
W1118_r2_a_mRNA_kmer  <- W1118_r2_a_mRNA_kmer  %>% mutate(call_code = as.character(call_code))
W1118_r3_a_mRNA_kmer <- W1118_r3_a_mRNA_kmer %>% mutate(call_code = as.character(call_code))

FR1113N_r1_a_mRNA_kmer <- FR1113N_r1_a_mRNA_kmer %>% mutate(call_code = as.character(call_code))
FR1113N_r2_a_mRNA_kmer <- FR1113N_r2_a_mRNA_kmer %>% mutate(call_code = as.character(call_code))
FR1113N_r3_a_mRNA_kmer <- FR1113N_r3_a_mRNA_kmer %>% mutate(call_code = as.character(call_code))

Mettl3_r1_a_mRNA_kmer <- Mettl3_r1_a_mRNA_kmer %>% mutate(call_code = as.character(call_code))
Mettl3_r2_a_mRNA_kmer <- Mettl3_r2_a_mRNA_kmer %>% mutate(call_code = as.character(call_code))


# Build one table with an explicit 'rep' column
all_reads <- bind_rows(
  W1118_r1_a_mRNA_kmer %>% mutate(sample="w1118", rep="Rep1"),
  W1118_r2_a_mRNA_kmer %>% mutate(sample="w1118", rep="Rep2"),
  W1118_r3_a_mRNA_kmer %>% mutate(sample="w1118", rep="Rep3"),
  FR1113N_r1_a_mRNA_kmer %>% mutate(sample="FR113N", rep="Rep1"),
  FR1113N_r2_a_mRNA_kmer %>% mutate(sample="FR113N", rep="Rep2"),
  FR1113N_r3_a_mRNA_kmer %>% mutate(sample="FR113N", rep="Rep3"),
  Mettl3_r1_a_mRNA_kmer %>% mutate(sample="Mettl3-KO", rep="Rep1"),
  Mettl3_r2_a_mRNA_kmer %>% mutate(sample="Mettl3-KO", rep="Rep2")) %>%
  mutate(call_prob = as.numeric(call_prob),
         read_length = as.numeric(read_length),
         sample = factor(sample, levels = c("w1118","FR113N","Mettl3-KO")),
         rep = factor(rep, levels = c("Rep1","Rep2","Rep3")))

# Keep PASS reads and sensible lengths
all_reads_pass <- all_reads %>%
  filter(!is.na(read_length), read_length > 0, !fail)

# Quick sanity check
dplyr::count(all_reads_pass, sample, rep)

# Colors per sample
cols <- c("w1118"="#1b9e77","FR113N"="#d95f02","Mettl3-KO"="#7570b3")

# Plot: x = replicate, facet by sample
read_length_by_rep <- ggplot(all_reads_pass, aes(x = rep, y = read_length, fill = sample)) +
  geom_violin(scale = "width", trim = TRUE, linewidth = 0.3, width = 0.9) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.3,
               fill = "white", alpha = 0.6) +
  facet_wrap(~ sample, nrow = 1) +
  scale_y_log10(labels = scales::label_comma()) +
  scale_fill_manual(values = cols, guide = "none") +
  labs(title = "Read length by replicate (PASS reads)", x = "Replicate", y = "Read length (bp, log10)") + 
  theme_classic() +
  theme(panel.grid.minor = element_blank(), strip.background = element_rect(fill = "grey95", colour = NA), strip.text = element_text(face = "bold"), plot.title = element_text(face = "bold", hjust = 0)) +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(plot.title= element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title= element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(hjust = 1),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text  = element_text(size = 12, family = "Times New Roman"))
# Save
dir.create("fig2", showWarnings = FALSE)
#ggsave("fig2/read_length_by_rep.png", read_length_by_rep,width = 8, height = 3, units = "in", dpi = 300, bg = "white")
ggsave("fig2/read_length_by_rep.pdf", plot = read_length_by_rep, width = 5, height = 5, units = "in", device = cairo_pdf, bg = "white")


# kmer distribution
library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(purrr)

# 1) Normalize each table so bind_rows doesn't choke
normalize_modkit <- function(df) {
  df %>%
    mutate(
      call_code = as.character(call_code),
      query_kmer = as.character(query_kmer),
      ref_kmer = as.character(ref_kmer),
      canonical_base = as.character(canonical_base),
      modified_primary_base = as.character(modified_primary_base),
      base_qual = suppressWarnings(as.numeric(base_qual)),
      read_length = suppressWarnings(as.numeric(read_length)),
      call_prob = suppressWarnings(as.numeric(call_prob)),
      flag = suppressWarnings(as.integer(flag)))
}

inputs <- list(
  # w1118
  W1118_r1_a_mRNA_kmer %>% mutate(sample="w1118", rep="Rep1", mod="m6A"),
  W1118_r2_a_mRNA_kmer %>% mutate(sample="w1118", rep="Rep2", mod="m6A"),
  W1118_r3_a_mRNA_kmer %>% mutate(sample="w1118", rep="Rep3", mod="m6A"),
  W1118_r1_m_mRNA_kmer %>% mutate(sample="w1118", rep="Rep1", mod="m5C"),
  W1118_r2_m_mRNA_kmer %>% mutate(sample="w1118", rep="Rep2", mod="m5C"),
  W1118_r3_m_mRNA_kmer %>% mutate(sample="w1118", rep="Rep3", mod="m5C"),
  W1118_r1_pseu_mRNA_kmer %>% mutate(sample="w1118", rep="Rep1", mod="psi"),
  W1118_r2_pseu_mRNA_kmer %>% mutate(sample="w1118", rep="Rep2", mod="psi"),
  W1118_r3_pseu_mRNA_kmer %>% mutate(sample="w1118", rep="Rep3", mod="psi"),
  # FR113N
  FR1113N_r1_a_mRNA_kmer %>% mutate(sample="FR113N", rep="Rep1", mod="m6A"),
  FR1113N_r2_a_mRNA_kmer %>% mutate(sample="FR113N", rep="Rep2", mod="m6A"),
  FR1113N_r3_a_mRNA_kmer %>% mutate(sample="FR113N", rep="Rep3", mod="m6A"),
  FR1113N_r1_m_mRNA_kmer %>% mutate(sample="FR113N", rep="Rep1", mod="m5C"),
  FR1113N_r2_m_mRNA_kmer %>% mutate(sample="FR113N", rep="Rep2", mod="m5C"),
  FR1113N_r3_m_mRNA_kmer %>% mutate(sample="FR113N", rep="Rep3", mod="m5C"),
  FR1113N_r1_pseu_mRNA_kmer %>% mutate(sample="FR113N", rep="Rep1", mod="psi"),
  FR1113N_r2_pseu_mRNA_kmer %>% mutate(sample="FR113N", rep="Rep2", mod="psi"),
  FR1113N_r3_pseu_mRNA_kmer %>% mutate(sample="FR113N", rep="Rep3", mod="psi"),
  # Mettl3-KO
  Mettl3_r1_a_mRNA_kmer %>% mutate(sample="Mettl3-KO", rep="Rep1", mod="m6A"),
  Mettl3_r2_a_mRNA_kmer %>% mutate(sample="Mettl3-KO", rep="Rep2", mod="m6A"),
  Mettl3_r1_m_mRNA_kmer %>% mutate(sample="Mettl3-KO", rep="Rep1", mod="m5C"),
  Mettl3_r2_m_mRNA_kmer %>% mutate(sample="Mettl3-KO", rep="Rep2", mod="m5C"),
  Mettl3_r1_pseu_mRNA_kmer%>% mutate(sample="Mettl3-KO", rep="Rep1", mod="psi"),
  Mettl3_r2_pseu_mRNA_kmer%>% mutate(sample="Mettl3-KO", rep="Rep2", mod="psi"))

all_reads <- inputs %>% map(normalize_modkit) %>% bind_rows() %>%
  mutate(
    sample = factor(sample, levels = c("w1118","FR113N","Mettl3-KO")),
    rep = factor(rep, levels = c("Rep1","Rep2","Rep3")),
    mod = factor(mod, levels = c("m6A","m5C","psi")))

# 2) Build a robust is_fail flag
all_reads <- all_reads %>%
  mutate(
    # if 'fail' column is missing, create NA
    fail = if ("fail" %in% names(.)) fail else NA,
    is_fail = case_when(
      is.logical(fail) ~ fail,
      is.numeric(fail) ~ fail != 0,
      TRUE ~ tolower(as.character(fail)) %in% c("true","t","yes","y","1")))

# Optional: inspect how many PASS/FAIL/NA
print(table(all_reads$is_fail, useNA = "ifany"))

# 3) Keep PASS (treat NA as PASS; change to 'FALSE' if you want to drop unknowns)
all_reads_pass <- all_reads %>%
  filter(is.na(is_fail) | is_fail == FALSE) %>%
  filter(!is.na(read_length), read_length > 0)

# Sanity check
nrow(all_reads_pass)
dplyr::count(all_reads_pass, sample, mod)

# 4) Build 5-mer table and plot (faceted mod ~ sample)
TOP_N <- 18

kmers_all <- all_reads_pass %>%
  filter(!is.na(query_kmer), query_kmer != ".") %>%
  mutate(
    kmer = toupper(query_kmer),
    kmer = if_else(str_detect(kmer, "^[ACGT]{5}$"), kmer, NA_character_)) %>%
  filter(!is.na(kmer)) %>%
  mutate(
    is_DRACH = (mod == "m6A") & str_detect(kmer, "^[AGT][AG]A[CT][ACT]$"),
    motif = if_else(is_DRACH, "DRACH", "Other"))

topk <- kmers_all %>%
  count(sample, mod, kmer, motif, sort = TRUE) %>%
  group_by(sample, mod) %>%
  slice_max(n, n = TOP_N, with_ties = FALSE) %>%
  mutate(kmer = fct_reorder(kmer, n)) %>%
  ungroup()

stopifnot(nrow(topk) > 0)  # will error early if still empty

kmer_lollipop_by_mod_sample <- ggplot(topk, aes(x = kmer, y = n, color = motif)) +
  geom_segment(aes(xend = kmer, y = 0, yend = n), linewidth = 0.35, alpha = 0.85) +
  geom_point(size = 1.8, alpha = 0.98) +
  coord_flip() +
  facet_grid(mod ~ sample, scales = "free_y") +
  scale_color_manual(values = c("DRACH" = "#e41a1c", "Other" = "grey50")) +
  labs(title = paste0("Top ", TOP_N, " query 5-mers per modification (DRACH highlighted for m6A)"),
       x = "k-mer (5-mer)", y = "Count", color = "Motif", family = "Times New Roman") +
  theme_classic(base_size = 12) +
  theme(
    panel.grid.major.y = element_blank(),
    axis.text.y = element_text(size = 7),
    strip.text = element_text(face = "bold"),
    plot.title = element_text(face = "bold", hjust = 0.5))


split kmer by mod:
  
  library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(tidytext)

dir.create("fig2", showWarnings = FALSE)
TOP_N <- 18


make_kmer_lollipop_by_mod <- function(dat, mod_name, top_n = TOP_N, out_dir = "fig2") {
  dat_mod <- dat %>%
    filter(mod == mod_name, !is.na(query_kmer), query_kmer != ".") %>%
    mutate(kmer = toupper(query_kmer),
           kmer = if_else(str_detect(kmer, "^[ACGT]{5}$"), kmer, NA_character_)) %>%
    filter(!is.na(kmer)) %>%
    mutate(is_DRACH = (mod_name == "m6A") & str_detect(kmer, "^[AGT][AG]A[CT][ACT]$"),
           motif = if_else(is_DRACH, "DRACH", "Other"))
  
  topk <- dat_mod %>%
    count(sample, kmer, motif, sort = TRUE) %>%
    group_by(sample) %>%
    slice_max(n, n = top_n, with_ties = FALSE) %>%
    ungroup() %>%
    # <- THIS is the key: order kmers within each sample
    mutate(kmer_reordered = reorder_within(kmer, n, sample))
  
  if (nrow(topk) == 0) {
    stop(paste0("No k-mers to plot for ", mod_name, ". Check PASS filtering or query_kmer content."))
  }
  
  p <- ggplot(topk, aes(x = kmer_reordered, y = n, color = motif)) +
    geom_segment(aes(xend = kmer_reordered, y = 0, yend = n), linewidth = 0.35, alpha = 0.85, color = "black") +
    geom_point(size = 3, alpha = 0.99) +
    coord_flip() +
    facet_wrap(~ sample, nrow = 1, scales = "free_y") +
    scale_x_reordered() +
    scale_color_manual(values = c("DRACH" = "darkorange", "Other" = "darkorange")) +
    labs(title = paste0("Top ", top_n, " query 5-mers — ", mod_name), x = "k-mer (5-mer)", y = "Count", color = "Motif", family = "Times New Roman") +
    theme_classic(base_size = 12, base_family = "Times New Roman") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text  = element_text(size = 12),
          axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "right",
          legend.title = element_text(face = "bold"),
          legend.text  = element_text(size = 12, family = "Times New Roman"))
  
  pdf_file <- file.path(out_dir, paste0("kmer_lollipop_", mod_name, "_by_sample.pdf"))
  ggsave(pdf_file, plot = p, width = 9, height = 9, units = "in",
         device = cairo_pdf, bg = "white")
  
  p
}


# If all_reads_pass isn’t already PASS-only, uncomment:
# all_reads_pass <- all_reads_pass %>% filter(is.na(is_fail) | is_fail == FALSE)

p_m6A <- make_kmer_lollipop_by_mod(all_reads_pass, "m6A")
p_m5C <- make_kmer_lollipop_by_mod(all_reads_pass, "m5C")
p_psi <- make_kmer_lollipop_by_mod(all_reads_pass, "psi")



Differential modification analysis Fig3:
  
  library(readr)
library(dplyr)
library(stringr)

# choose the first existing column from a preference list
.pick_col <- function(df, prefs) {
  hit <- intersect(prefs, names(df))
  if (length(hit) == 0) NULL else hit[[1]]
}

load_dmr <- function(path, base, effect_thresh = 0.05, pval_thresh = 0.05, b_max_pct_if_CT = 10) {
  df <- read_tsv(path, col_names = TRUE, show_col_types = FALSE)
  
  # column detection with fallbacks
  chrom_col <- .pick_col(df, c("#chrom","X.chrom","chrom"))
  start_col <- .pick_col(df, c("start","chromStart"))
  end_col <- .pick_col(df, c("end","chromEnd"))
  strand_col <- .pick_col(df, c("strand"))
  
  pval_col <- .pick_col(df, c("balanced_map_pvalue","map_pvalue","pvalue","p_value"))
  eff_col <- .pick_col(df, c("cohen_h","balanced_effect_size","effect_size"))
  a_pct_col <- .pick_col(df, c("a_pct_modified"))
  b_pct_col <- .pick_col(df, c("b_pct_modified"))
  
  # basic presence check
  need <- c(chrom_col, start_col, end_col, strand_col, eff_col, a_pct_col, b_pct_col)
  missing <- setdiff(need, names(df))
  if (length(missing) > 0) {
    stop("Missing required columns in ", basename(path), ": ", paste(missing, collapse = ", "))
  }
  
  # numeric coercion (guard against character parsing)
  num_cols <- c(eff_col, a_pct_col, b_pct_col, pval_col)
  num_cols <- num_cols[!is.null(num_cols)]
  df <- df %>% mutate(across(all_of(num_cols), suppressWarnings(as.numeric)))
  
  # effect-size filter (always, if available)
  if (!is.null(eff_col)) {
    df <- df %>% filter(abs(.data[[eff_col]]) > effect_thresh)
  } else {
    warning(basename(path), ": effect size column not found; skipping effect-size filter.")
  }
  
  # p-value filter (only if available)
  if (!is.null(pval_col)) {
    df <- df %>% filter(.data[[pval_col]] < pval_thresh)
  } else {
    warning(basename(path), ": p-value column not found; skipping p-value filter.")
  }
  
  # A: both > 0 ;  C/T: A > 0 and 0 ≤ B ≤ 10
  if (base == "A") {
    df <- df %>% filter(.data[[a_pct_col]] > 0, .data[[b_pct_col]] > 0)
  } else if (base %in% c("C","T")) {
    df <- df %>% filter(.data[[a_pct_col]] > 0, .data[[b_pct_col]] >= 0, .data[[b_pct_col]] <= b_max_pct_if_CT)
  } else {
    warning("Unknown base '", base, "'. No base-specific % filters applied.")
  }
  
  # Location filter (case-insensitive) if present
  if ("Location" %in% names(df)) {
    bad_locs <- c("intron","null","N/A","")
    df <- df %>%
      mutate(Location_lc = str_to_lower(str_trim(Location))) %>%
      dplyr::filter(!Location_lc %in% bad_locs) %>%
      dplyr::select(-Location_lc)
  }
  
  df
}

filter(a_total >= 5, b_total >= 5, a_pct_modified >= 10, b_pct_modified >= 10)


W1118_vs_Mettl3_KO_a_dmr <- read_tsv("/mnt/raid1/george/direct_RNA/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/pair_joint/W1118_vs_Mettl3KO.single_site_joint_A.tsv") %>% 
  dplyr::filter(abs(cohen_h) > 0.05, a_total >= 1, b_total >= 1) %>% 
  dplyr::filter(Location != "intron", Location != "null", Location != "N/A", Location != " ")

W1118_vs_Mettl3_m_dmr <- read_tsv("/mnt/raid1/george/direct_RNA/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/pair_joint/W1118_vs_Mettl3KO.single_site_joint_C.tsv") %>% 
  dplyr::filter(abs(cohen_h) > 0.05, a_total >= 1200, b_total >= 1200) %>% 
  dplyr::filter(Location != "Intron", Location != "null", Location != "N/A", Location != " ")

W1118_vs_Mettl3_psi_dmr <- read_tsv("/mnt/raid1/george/direct_RNA/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/pair_joint/W1118_vs_Mettl3KO.single_site_joint_T.tsv") %>% 
  dplyr::filter(abs(cohen_h) > 0.05, a_total >= 250, b_total >= 250) %>% 
  dplyr::filter(Location != "intron", Location != "null", Location != "N/A", Location != " ")

FR1113_vs_Mettl3_KO_a_dmr <- read_tsv("/mnt/raid1/george/direct_RNA/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/pair_joint/FR113N_vs_Mettl3KO.single_site_joint_A.tsv") %>% 
  dplyr::filter(abs(cohen_h) > 0.05, a_total >= 1, b_total >= 1)  %>% 
  dplyr::filter(Location != "intron", Location != "null", Location != "N/A", Location != " ")

FR1113_vs_Mettl3_KO_m_dmr <- read_tsv("/mnt/raid1/george/direct_RNA/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/pair_joint/FR113N_vs_Mettl3KO.single_site_joint_C.tsv") %>% 
  dplyr::filter(abs(cohen_h) > 0.05, a_total >= 1200, b_total >= 1200) %>% 
  dplyr::filter(Location != "intron", Location != "null", Location != "N/A", Location != " ")

FR1113_vs_Mettl3_KO_psi_dmr <- read_tsv("/mnt/raid1/george/direct_RNA/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/pair_joint/FR113N_vs_Mettl3KO.single_site_joint_T.tsv") %>% 
  dplyr::filter(abs(cohen_h) > 0.05, a_total >= 350, b_total >= 350) %>% 
  dplyr::filter(Location != "intron", Location != "null", Location != "N/A", Location != " ")


FR113_vs_w1118_a_dmr <- read_tsv("/mnt/raid1/george/direct_RNA/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/pair_joint/FR113N_vs_W1118.single_site_joint_A.tsv") %>% 
  dplyr::filter(abs(cohen_h) > 0.05, a_total >= 1, b_total >= 1)  %>% 
  dplyr::filter(Location != "intron", Location != "null", Location != "N/A", Location != " ")

FR113_vs_w1118_m_dmr <- read_tsv("/mnt/raid1/george/direct_RNA/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/pair_joint/FR113N_vs_W1118.single_site_joint_C.tsv") %>% 
  dplyr::filter(abs(cohen_h) > 0.05, a_total >= 800, b_total >= 800) %>% 
  dplyr::filter(Location != "intron", Location != "null", Location != "N/A", Location != " ")

FR113_vs_w1118_psi_dmr <- read_tsv("/mnt/raid1/george/direct_RNA/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/pair_joint/FR113N_vs_W1118.single_site_joint_T.tsv") %>% 
  dplyr::filter(abs(cohen_h) > 0.05, a_total >=350, b_total >= 350) %>% 
  dplyr::filter(Location != "intron", Location != "null", Location != "N/A", Location != " ")

library(dplyr)
library(ggplot2)
library(stringr)
library(purrr)
library(tibble)

list_dmr <- list(
  `w1118 vs Mettl3-KO m6A` = W1118_vs_Mettl3_KO_a_dmr,
  `w1118 vs Mettl3-KO m5C` = W1118_vs_Mettl3_m_dmr,
  `w1118 vs Mettl3-KO psi` = W1118_vs_Mettl3_psi_dmr,
  `FR1113 vs Mettl3-KO m6A` = FR1113_vs_Mettl3_KO_a_dmr,
  `FR1113 vs Mettl3-KO m5C` = FR1113_vs_Mettl3_KO_m_dmr,
  `FR1113 vs Mettl3-KO psi` = FR1113_vs_Mettl3_KO_psi_dmr,
  `FR113 vs w1118 m6A` = FR113_vs_w1118_a_dmr,
  `FR113 vs w1118 m5C` = FR113_vs_w1118_m_dmr,
  `FR113 vs w1118 psi` = FR113_vs_w1118_psi_dmr)


dmr_combined <- imap_dfr(list_dmr, ~ {
  mod <- case_when(
    str_detect(.y, regex("\\bm6A\\b", ignore_case = TRUE)) ~ "m6A",
    str_detect(.y, regex("\\bm5C\\b", ignore_case = TRUE)) ~ "m5C",
    str_detect(.y, regex("psi|Ψ", ignore_case = TRUE)) ~ "psi",
    TRUE ~ "unknown")
  .x %>% mutate(Dataset = .y, Modification = mod)
})

#Count DMRs per dataset & modification
dmr_counts <- dmr_combined %>%
  group_by(Dataset, Modification) %>%
  summarise(DMR_Count = n(), .groups = "drop") %>%
  # order bars by total DMRs (optional)
  group_by(Dataset) %>%
  mutate(Total = sum(DMR_Count)) %>%
  ungroup() %>%
  mutate(Dataset = fct_reorder(Dataset, Total)) %>%
  arrange(Dataset)

# Lollipop plot (grouped by modification)
colors <- c("m6A" = "#95CC5E", "m5C" = "#748AA6", "psi" = "#EAB8D1")

pd <- position_dodge(width = 0.7)

p_lollipop_single <- ggplot(dmr_counts, aes(x = Dataset)) +
  geom_linerange(aes(ymin = 0, ymax = DMR_Count, color = Modification), position = pd, linewidth = 0.9, alpha = 0.9, color = "black") +
  geom_point(aes(y = DMR_Count, fill = Modification, color = Modification), shape = 21, stroke = 1.2, size = 5, position = pd) +
  geom_text(aes(y = DMR_Count, label = DMR_Count), position = pd, vjust = -1.0, size = 3.8, family = "Times New Roman") +
  scale_color_manual(values = colors) +
  scale_fill_manual(values = colors) +
  labs(title = "Number of differentially modified sites", x = NULL, y = "DMR count", color = "Modification", fill = "Modification",  family = "Times New Roman") +
  theme_classic(base_size = 13) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.text.y = element_text(size = 10), plot.title = element_text(face = "bold", hjust = 0.5)) +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "right",
        legend.text = element_text(size = 12, family = "Times New Roman"))

ggsave("fig2/num_diff_modified.pdf", plot = p_lollipop_single, width = 9, height = 9, units = "in", device = cairo_pdf, bg = "white")


effect size distribution Fig2:
  
  library(dplyr)
library(ggplot2)
library(stringr)
library(forcats)
library(scales)
library(gridExtra)

# --- 1) Function to normalize regions ---
normalize_region <- function(x) {
  y <- tolower(trimws(x))
  y <- gsub("[ _-]", "", y)
  dplyr::case_when(
    str_detect(y, "^(5'utr|utr5'|5utr|utr5|fiveprimeutr|fiveprime)$") ~ "5'UTR",
    str_detect(y, "^(cds|codingsequence|cdsregion)$") ~ "CDS",
    str_detect(y, "^(3'utr|utr3'|3utr|utr3|threeprimeutr|threeprime)$") ~ "3'UTR",
    TRUE ~ NA_character_)
}

# --- 2) Combine DMR datasets ---
list_dmr <- list(
  W1118_vs_Mettl3_KO_a_dmr,
  W1118_vs_Mettl3_m_dmr,
  W1118_vs_Mettl3_psi_dmr,
  FR1113_vs_Mettl3_KO_a_dmr,
  FR1113_vs_Mettl3_KO_m_dmr,
  FR1113_vs_Mettl3_KO_psi_dmr,
  FR113_vs_w1118_a_dmr,
  FR113_vs_w1118_m_dmr,
  FR113_vs_w1118_psi_dmr)

names(list_dmr) <- c(
  "w1118 vs Mettl3-KO m6A", "w1118 vs Mettl3-KO m5C", "w1118 vs Mettl3-KO psi",
  "FR1113 vs Mettl3-KO m6A", "FR1113 vs Mettl3-KO m5C", "FR1113 vs Mettl3-KO psi",
  "FR113 vs w1118 m6A", "FR113 vs w1118 m5C", "FR113 vs w1118 psi")

# Add dataset & modification type
dmr_combined <- bind_rows(lapply(names(list_dmr), function(n) {
  df <- list_dmr[[n]]
  mod <- ifelse(grepl("_a_", n), "m6A", ifelse(grepl("_m_", n), "m5C", "psi"))
  df %>% mutate(Dataset = n, Modification = mod)
}))


dmr_combined <- dmr_combined %>%
  mutate(Region = normalize_region(Location)) %>%
  dplyr::filter(!is.na(Region))

# --- 4) Calculate percent DMRs per region per dataset ---
region_counts <- dmr_combined %>%
  group_by(Dataset, Modification, Region) %>%
  summarise(Count = n(), .groups = "drop")

totals <- region_counts %>%
  group_by(Dataset, Modification) %>%
  summarise(Total = sum(Count), .groups = "drop")

region_pct <- region_counts %>%
  left_join(totals, by = c("Dataset", "Modification")) %>%
  mutate(percent = Count / Total)

# --- 5) Barplot: percent of DMRs per region ---
region_colors <- c("5'UTR" = "#268BD2", "CDS" = "#DC322F", "3'UTR" = "#2AA198")

p_region_bar <- ggplot(region_pct, aes(x = Dataset, y = percent, fill = Region)) +
  geom_col(position = "stack", color = "black") +
  facet_wrap(~ Modification, nrow = 1) +
  scale_fill_manual(values = region_colors) +
  scale_y_continuous(labels = percent_format(accuracy = 1)) +
  labs(title = "Percent of DMRs in 5'UTR, CDS, 3'UTR by dataset and modification",
       x = "Dataset", y = "Percent of DMRs") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# --- 6) Boxplot: effect size by region ---
dmr_combined$Region <- factor(dmr_combined$Region, levels = c("5'UTR", "CDS", "3'UTR"))
p_region_box <- ggplot(dmr_combined, aes(x = Region, y = effect_size, fill = Region)) +
  geom_boxplot(outlier.shape = NA, notch = T) + 
  #geom_violin(trim = FALSE, alpha = 0.9) +
  #geom_jitter(aes(color = Region), width = 0.2, alpha = 0.5, size = 1.5) +
  #geom_jitter(color = "black", width = 0.05, alpha = 0.05, size = 1.5) +
  facet_wrap(~ Modification + Dataset, scales = "free") +
  scale_fill_manual(values = region_colors) +
  scale_color_manual(values = region_colors) +
  labs(title = "Effect size distribution by region", x = "Region", y = "Effect size") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_text(size = 10))

# --- Display plots ---
grid.arrange(p_region_bar, nrow = 1)
print(p_region_box)

ML 
Lollipop plot of DMR counts per dataset and modification:
  
  
  
  # summarize one DMR table to gene level --------
summarize_gene <- function(df, contrast_label, base_label, fdr_cut = 0.05, min_abs_es = 0.0) {
  df %>%
    #mutate(FDR = p.adjust(balanced_map_pvalue, method = "BH")) %>%
    #filter(!is.na(Gene_ID), FDR < fdr_cut, abs(effect_size) >= min_abs_es) %>%
    filter(!is.na(Gene_ID), abs(effect_size) >= min_abs_es) %>%
    group_by(Gene_ID) %>%
    summarise(n_sites = n(), median_es = median(effect_size, na.rm = TRUE), top_loc = names(sort(table(Location), decreasing = TRUE))[1], 
              .groups = "drop") %>%
    mutate(contrast = contrast_label, base = base_label)
}

gene_all <- bind_rows(
  summarize_gene(W1118_vs_Mettl3_KO_a_dmr, "w1118 vs Mettl3-KO", "m6A"),
  summarize_gene(W1118_vs_Mettl3_m_dmr, "w1118 vs Mettl3-KO", "m5C"),
  summarize_gene(W1118_vs_Mettl3_psi_dmr, "w1118 vs Mettl3-KO", "psi"),
  summarize_gene(FR1113_vs_Mettl3_KO_a_dmr, "FR113N vs Mettl3-KO", "m6A"),
  summarize_gene(FR1113_vs_Mettl3_KO_m_dmr, "FR113N vs Mettl3-KO", "m5C"),
  summarize_gene(FR1113_vs_Mettl3_KO_psi_dmr, "FR113N vs Mettl3-KO", "psi"),
  summarize_gene(FR113_vs_w1118_a_dmr, "FR113N vs w1118", "m6A"),
  summarize_gene(FR113_vs_w1118_m_dmr, "FR113N vs w1118", "m5C"),
  summarize_gene(FR113_vs_w1118_psi_dmr, "FR113N vs w1118", "psi")) %>%
  mutate(contrast = factor(contrast, levels = c("w1118 vs Mettl3-KO","FR113N vs Mettl3-KO","FR113N vs w1118")),
         base = factor(base, levels = c("m6A","m5C", "psi")),
         direction = ifelse(median_es >= 0, "Higher in A", "Higher in B"))

# pick top genes per facet by |median_es| --------
topN <- 10
plot_df <- gene_all %>%
  group_by(contrast, base) %>%
  slice_max(order_by = abs(median_es), n = topN, with_ties = FALSE) %>%
  ungroup() %>%
  mutate(Gene_ID_fac = tidytext::reorder_within(Gene_ID, median_es, interaction(contrast, base)))


col_dir <- c("Higher in A" = "#1b9e77",   # A group higher (left side of label)
             "Higher in B" = "#d95f02")   # B group higher


topdiff <- ggplot(plot_df,aes(x = Gene_ID_fac, y = median_es, color = direction)) +
  geom_segment(aes(xend = Gene_ID_fac, y = 0, yend = median_es), linewidth = 0.6, alpha = 0.8, color = "black") +
  geom_point(aes(size = n_sites), stroke = 1.2, alpha = 0.95, size = 3) +
  scale_color_manual(values = col_dir, name = "Direction") +
  scale_size_continuous(name = "# sites", range = c(1.2, 4.0)) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  coord_flip() +
  facet_grid(base ~ contrast, scales = "free_y", space = "free_y") +
  tidytext::scale_x_reordered() +
  labs(x = NULL, y = expression("Median " * Delta * " modified fraction (A − B) per gene"), title = "Top genes by net change across contrasts and base types", subtitle = "Points sized by # significant sites; color indicates which group has higher modification (A vs B)", vjust = -1.2, size = 6, fontface = "bold", family = "Times New Roman") +
  theme_classic(base_size = 11) +
  theme(legend.position = "right", strip.background = element_blank(), strip.text = element_text(face = "bold"), panel.grid = element_blank()) +
  theme_classic(base_size = 11) +
  theme(legend.position = "right", panel.spacing.x = unit(8, "pt"), strip.background = element_blank(), strip.text = element_text(face = "bold")) +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(axis.text.x = element_text(hjust = 1),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 12),
        legend.position = "right",
        legend.text = element_text(size = 12, family = "Times New Roman"))


make box plot of the effect size:
  
  library(dplyr)
library(stringr)
library(ggplot2)
library(forcats)

## 1. Prepare data for plotting

dmr_for_plot <- dmr_combined %>%
  # effect size column (rename for clarity)
  rename(EffectSize = cohen_h) %>%
  filter(!is.na(EffectSize)) %>%
  
  # clean up comparison labels (Dataset -> Comparison)
  mutate(Comparison = case_when(
    str_detect(Dataset, "w1118 vs Mettl3-KO") ~ "w1118 vs Mettl3-KO",
    str_detect(Dataset, "FR1113 vs Mettl3-KO") |
      str_detect(Dataset, "FR113 vs Mettl3-KO") ~ "FR113N vs Mettl3-KO",
    str_detect(Dataset, "FR113 vs w1118") ~ "FR113N vs w1118",
    TRUE ~ Dataset),
    # order comparisons the way you want on the x-axis
    Comparison = factor(Comparison, levels = c("FR113N vs w1118", "w1118 vs Mettl3-KO", "FR113N vs Mettl3-KO")),
    # nice labels for modifications
    Modification = factor(Modification, levels = c("m6A","m5C","psi"), labels = c("m⁶A","m⁵C","Ψ")))

## 2. Colors for each modification 
mod_cols <- c("m⁶A" = "#95CC5E", "m⁵C" = "#748AA6", "Ψ" = "#EAB8D1")

## 3. Boxplot of effect sizes

box_eff <- ggplot(dmr_for_plot, aes(x = Comparison, y = EffectSize, fill = Modification)) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.4) +
  geom_boxplot(width = 0.7, outlier.alpha = 0.2, position = position_dodge(width = 0.8),linewidth = 0.4, notch=T) +
  scale_fill_manual(values = mod_cols, name = "Modification") +
  labs(x = NULL, y = "Effect size (Δ modification fraction, A − B)") +
  coord_cartesian(ylim = c(-1, 1)) +  # tweak if you want a different range
  theme_classic(base_size = 12) +
  theme(axis.text.x  = element_text(angle = 20, hjust = 1), legend.position = "top", legend.title = element_text(face = "bold"),
        panel.border = element_rect(colour = "black", fill = NA, linewidth = 0.6))

ggsave("fig2/box_eff_forall.pdf", plot = box_eff, width = 5, height = 5, units = "in", device = cairo_pdf, bg = "white")


Split the lollipop plot :
  
  library(dplyr)
library(ggplot2)
library(tidytext)
library(grid)   # for unit()

# keep your original colors
col_dir <- c("Higher in A" = "#1b9e77",
             "Higher in B" = "#d95f02")

make_topdiff_plot_one_contrast <- function(gene_all, contrast_name, topN = 15) {
  
  plot_df <- gene_all %>%
    dplyr::filter(contrast == contrast_name) %>%
    group_by(contrast, base) %>%
    slice_max(order_by = abs(median_es), n = topN, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(Gene_ID_fac = tidytext::reorder_within(Gene_ID, median_es, interaction(contrast, base)))
  
  ggplot(plot_df, aes(x = Gene_ID_fac, y = median_es, color = direction)) +
    geom_segment(aes(xend = Gene_ID_fac, y = 0, yend = median_es),linewidth = 0.6, alpha = 0.8, color = "black") +
    geom_point(aes(size = n_sites), stroke = 1.2, alpha = 0.95, size = 3) +
    scale_color_manual(values = col_dir, name = "Direction") +
    scale_size_continuous(name = "# sites", range = c(1.2, 4.0)) +
    geom_hline(yintercept = 0, linewidth = 0.3) +
    coord_flip() +
    facet_grid(base ~ contrast, scales = "free_y", space = "free_y") +
    tidytext::scale_x_reordered() +
    labs(x = NULL, y = expression("Median " * Delta * " modified fraction (A − B) per gene"),
         title = paste0("Top genes by net change across base types: ", contrast_name),
         subtitle = "Points sized by # significant sites; color indicates which group has higher modification (A vs B)") +
    theme_classic(base_size = 11) +
    theme(legend.position = "right",
          strip.background = element_blank(),
          strip.text = element_text(face = "bold"),
          panel.grid = element_blank()) +
    theme_classic(base_size = 11) +
    theme(legend.position = "right",
          panel.spacing.x = unit(8, "pt"),
          strip.background = element_blank(),
          strip.text = element_text(face = "bold")) +
    theme_classic(base_size = 12, base_family = "Times New Roman") +
    theme(axis.text.x = element_text(hjust = 1),
          plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.position = "right",
          legend.text = element_text(size = 12, family = "Times New Roman"))
}

# make the 3 plots
p_w1118_mettl3 <- make_topdiff_plot_one_contrast(gene_all, "w1118 vs Mettl3-KO", topN = 15)
p_fr113_mettl3 <- make_topdiff_plot_one_contrast(gene_all, "FR113N vs Mettl3-KO", topN = 15)
p_fr113_w1118 <- make_topdiff_plot_one_contrast(gene_all, "FR113N vs w1118", topN = 15)

p_all <- (p_w1118_mettl3 | p_fr113_mettl3 | p_fr113_w1118) + plot_layout(guides = "collect") + plot_annotation(tag_levels = "A") &
  theme(legend.position = "right")

ggsave("fig2/topdiff_all_contrasts_patchwork_wide.pdf", p_all, width = 20, height = 9, units = "in", device = cairo_pdf, bg = "white")




Venn diagram of DMR-overlapping genes:
  
  sig_gene_set <- function(df, fdr_cut = 0.05, min_abs_es = 0.1, min_cov = 20) {
    df %>%
      #mutate(FDR = p.adjust(balanced_map_pvalue, "BH"),
      mutate(min_cov = pmin(a_total, b_total)) %>%
      #filter(!is.na(Gene_ID), FDR < fdr_cut, abs(effect_size) >= min_abs_es, min_cov >= min_cov) %>%
      filter(!is.na(Gene_ID), abs(effect_size) >= min_abs_es, min_cov >= min_cov) %>%
      distinct(Gene_ID, .keep_all = TRUE)
  }


# m6A (A-base)
A_WT_KO <- sig_gene_set(W1118_vs_Mettl3_KO_a_dmr)
A_FR_KO <- sig_gene_set(FR1113_vs_Mettl3_KO_a_dmr)
A_FR_WT <- sig_gene_set(FR113_vs_w1118_a_dmr)

sets_A <- list(`W1118 vs KO`= A_WT_KO$Gene_ID, `FR113N vs KO` = A_FR_KO$Gene_ID, `FR113N vs W1118` = A_FR_WT$Gene_ID)

sharem6A <- ggvenn(sets_A, fill_color = c("#1b9e77", "#d95f02", "#7570b3"), stroke_size = 0.6, set_name_size = 4, text_size = 4, show_percentage = FALSE)+ theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(plot.title= element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title= element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(hjust = 1),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text  = element_text(size = 12, family = "Times New Roman"))
# ggsave("venn_m6A_genes.pdf", width = 5, height = 5, useDingbats = FALSE)
ggsave("fig2/Fig2_sharem6A_diff.pdf", plot = sharem6A, width = 9, height = 9, units = "in", device = cairo_pdf, bg = "white")

# m5C (C-base)
C_WT_KO <- sig_gene_set(W1118_vs_Mettl3_m_dmr)
C_FR_KO <- sig_gene_set(FR1113_vs_Mettl3_KO_m_dmr)
C_FR_WT <- sig_gene_set(FR113_vs_w1118_m_dmr)

sets_C <- list(`W1118 vs KO` = C_WT_KO$Gene_ID, `FR113N vs KO` = C_FR_KO$Gene_ID, `FR113N vs W1118` = C_FR_WT$Gene_ID)

sharem5C <- ggvenn(sets_C, fill_color = c("#1b9e77", "#d95f02", "#7570b3"), stroke_size = 0.6, set_name_size = 4, text_size = 4, show_percentage = FALSE)+ theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(plot.title= element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title= element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(hjust = 1),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text  = element_text(size = 12, family = "Times New Roman"))

# psi (T-base)
P_WT_KO <- sig_gene_set(W1118_vs_Mettl3_psi_dmr)
P_FR_KO <- sig_gene_set(FR1113_vs_Mettl3_KO_psi_dmr)
P_FR_WT <- sig_gene_set(FR113_vs_w1118_psi_dmr)

sets_P <- list(`W1118 vs KO` = P_WT_KO$Gene_ID, `FR113N vs KO` = P_FR_KO$Gene_ID, `FR113N vs W1118` = P_FR_WT$Gene_ID)

sharepsi <- ggvenn(sets_P, fill_color = c("#1b9e77", "#d95f02", "#7570b3"), stroke_size = 0.6, set_name_size = 4, text_size = 4, show_percentage = FALSE)+ theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(plot.title= element_text(size = 12, face = "bold", hjust = 0.5),
        axis.title= element_text(size = 12),
        axis.text = element_text(size = 12),
        axis.text.x = element_text(hjust = 1),
        legend.position = "right",
        legend.title = element_text(face = "bold"),
        legend.text  = element_text(size = 12, family = "Times New Roman"))

#modkit dmr pair -a W1118_induro_mRNA_2025_102.sorted.bed.gz -a W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.sorted.bed.gz -a W1118_rep2_ref.sorted.bed.gz -b Mettl3_KO_mRNA_RNA004.sorted.bed.gz -b Mettl3_KO_mRNA_rep2.sorted.bed.gz -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/w1118_Mettl3_KO_single_base_replicates_T.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base T --threads 32 --header
#modkit dmr pair -a W1118_induro_mRNA_2025_102.sorted.bed.gz -a W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.sorted.bed.gz -a W1118_rep2_ref.sorted.bed.gz -b Mettl3_KO_mRNA_RNA004.sorted.bed.gz -b Mettl3_KO_mRNA_rep2.sorted.bed.gz -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/w1118_Mettl3_KO_single_base_replicates_A.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base A --threads 32 --header
#modkit dmr pair -a W1118_induro_mRNA_2025_102.sorted.bed.gz -a W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.sorted.bed.gz -a W1118_rep2_ref.sorted.bed.gz -b Mettl3_KO_mRNA_RNA004.sorted.bed.gz -b Mettl3_KO_mRNA_rep2.sorted.bed.gz -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/w1118_Mettl3_KO_single_base_replicates_C.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base C --threads 32 --header 


#modkit dmr pair -a FR113N_05_11_2025_rep2_ref.sorted.bed.gz -a FR113N_05_11_2025_rep3_ref.sorted.bed.gz -a FR113N_Oct_4_mRNA_induro.sorted.bed.gz -b W1118_induro_mRNA_2025_102.sorted.bed.gz -b W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.sorted.bed.gz -b W1118_rep2_ref.sorted.bed.gz -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/FR113N_W1118_single_base_replicates_T.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base T --threads 32 --header
#modkit dmr pair -a FR113N_05_11_2025_rep2_ref.sorted.bed.gz -a FR113N_05_11_2025_rep3_ref.sorted.bed.gz -a FR113N_Oct_4_mRNA_induro.sorted.bed.gz -b W1118_induro_mRNA_2025_102.sorted.bed.gz -b W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.sorted.bed.gz -b W1118_rep2_ref.sorted.bed.gz -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/FR113N_W1118_single_base_replicates_A.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base A --threads 32 --header
#modkit dmr pair -a FR113N_05_11_2025_rep2_ref.sorted.bed.gz -a FR113N_05_11_2025_rep3_ref.sorted.bed.gz -a FR113N_Oct_4_mRNA_induro.sorted.bed.gz -b W1118_induro_mRNA_2025_102.sorted.bed.gz -b W1118_Rep1_mRNA_m6A_m5C_pseu.aligned.sorted.bed.gz -b W1118_rep2_ref.sorted.bed.gz -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/FR113N_W1118_single_base_replicates_C.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base C --threads 32 --header

#modkit dmr pair -a FR113N_05_11_2025_rep2_ref.sorted.bed.gz -a FR113N_05_11_2025_rep3_ref.sorted.bed.gz -a FR113N_Oct_4_mRNA_induro.sorted.bed.gz -b Mettl3_KO_mRNA_RNA004.sorted.bed.gz -b Mettl3_KO_mRNA_rep2.sorted.bed.gz -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/FR113N_Mettl3_KO_single_base_replicates_T.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base T --threads 32 --header
#modkit dmr pair -a FR113N_05_11_2025_rep2_ref.sorted.bed.gz -a FR113N_05_11_2025_rep3_ref.sorted.bed.gz -a FR113N_Oct_4_mRNA_induro.sorted.bed.gz -b Mettl3_KO_mRNA_RNA004.sorted.bed.gz -b Mettl3_KO_mRNA_rep2.sorted.bed.gz -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/FR113N_Mettl3_KO_single_base_replicates_A.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base A --threads 32 --header
#modkit dmr pair -a FR113N_05_11_2025_rep2_ref.sorted.bed.gz -a FR113N_05_11_2025_rep3_ref.sorted.bed.gz -a FR113N_Oct_4_mRNA_induro.sorted.bed.gz -b Mettl3_KO_mRNA_RNA004.sorted.bed.gz -b Mettl3_KO_mRNA_rep2.sorted.bed.gz -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/FR113N_Mettl3_KO_single_base_replicates_C.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base C --threads 32 --header

modkit dmr multi -s FR113N_05_11_2025_rep2_ref.sorted.bed.gz FR113N_Rep1 -s FR113N_05_11_2025_rep3_ref.sorted.bed.gz FR113N_Rep2 -s FR113N_Oct_4_mRNA_induro.sorted.bed.gz FR113N_Rep3 -s Mettl3_KO_mRNA_RNA004.sorted.bed.gz Mettl3_Rep1 -s Mettl3_KO_mRNA_rep2.sorted.bed.gz Mettl3_Rep2 -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/FR113N_Mettl3_KO_single_base_replicates_T.bed -r /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/regions/FR113N_vs_Mettl3KO.filtered.exonic.regions.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base T --threads 32 --header
modkit dmr multi -s FR113N_05_11_2025_rep2_ref.sorted.bed.gz FR113N_Rep1 -s FR113N_05_11_2025_rep3_ref.sorted.bed.gz FR113N_Rep2 -s FR113N_Oct_4_mRNA_induro.sorted.bed.gz FR113N_Rep3 -s Mettl3_KO_mRNA_RNA004.sorted.bed.gz Mettl3_Rep1 -s Mettl3_KO_mRNA_rep2.sorted.bed.gz Mettl3_Rep2 -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/FR113N_Mettl3_KO_single_base_replicates_A.bed -r /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/regions/FR113N_vs_Mettl3KO.filtered.exonic.regions.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base A --threads 32 --header
modkit dmr multi -s FR113N_05_11_2025_rep2_ref.sorted.bed.gz FR113N_Rep1 -s FR113N_05_11_2025_rep3_ref.sorted.bed.gz FR113N_Rep2 -s FR113N_Oct_4_mRNA_induro.sorted.bed.gz FR113N_Rep3 -s Mettl3_KO_mRNA_RNA004.sorted.bed.gz Mettl3_Rep1 -s Mettl3_KO_mRNA_rep2.sorted.bed.gz Mettl3_Rep2 -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/FR113N_Mettl3_KO_single_base_replicates_C.bed -r /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/regions/FR113N_vs_Mettl3KO.filtered.exonic.regions.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base C --threads 32 --header

modkit dmr multi -s W1118_induro_mRNA_2025_102.sorted.bed.gz W1118_Rep1 -s W1118_rep2_ref.sorted.bed.gz W1118_Rep2 -s Mettl3_KO_mRNA_RNA004.sorted.bed.gz Mettl3_Rep1 -s Mettl3_KO_mRNA_rep2.sorted.bed.gz Mettl3_Rep2 -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/w1118_Mettl3_KO_single_base_replicates_T.bed -r /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/regions/W1118_vs_Mettl3KO.filtered.exonic.regions.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base T --threads 32 --header
modkit dmr multi -s W1118_induro_mRNA_2025_102.sorted.bed.gz W1118_Rep1 -s W1118_rep2_ref.sorted.bed.gz W1118_Rep2 -s Mettl3_KO_mRNA_RNA004.sorted.bed.gz Mettl3_Rep1 -s Mettl3_KO_mRNA_rep2.sorted.bed.gz Mettl3_Rep2 -o  /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/w1118_Mettl3_KO_single_base_replicates_A.bed -r /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/regions/W1118_vs_Mettl3KO.filtered.exonic.regions.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base A --threads 32 --header
modkit dmr multi -s W1118_induro_mRNA_2025_102.sorted.bed.gz W1118_Rep1 -s W1118_rep2_ref.sorted.bed.gz W1118_Rep2 -s Mettl3_KO_mRNA_RNA004.sorted.bed.gz Mettl3_Rep1 -s Mettl3_KO_mRNA_rep2.sorted.bed.gz Mettl3_Rep2 -o  /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/w1118_Mettl3_KO_single_base_replicates_C.bed -r /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/regions/W1118_vs_Mettl3KO.filtered.exonic.regions.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base C --threads 32 --header


modkit dmr multi -s FR113N_05_11_2025_rep2_ref.sorted.bed.gz FR113N_Rep1 -s FR113N_05_11_2025_rep3_ref.sorted.bed.gz FR113N_Rep2 -s FR113N_Oct_4_mRNA_induro.sorted.bed.gz FR113N_Rep3 -s W1118_induro_mRNA_2025_102.sorted.bed.gz W1118_Rep1 -s W1118_Rep2_mRNA_m6A_m5C_pseu.aligned.sorted.bed.gz W1118_Rep2 -s W1118_rep2_ref.sorted.bed.gz  W1118_Rep3 -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/FR113N_W1118_single_base_replicates_T.bed -r /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/regions/FR113N_vs_W1118.filtered.exonic.regions.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base T --threads 32 --header
modkit dmr multi -s FR113N_05_11_2025_rep2_ref.sorted.bed.gz FR113N_Rep1 -s FR113N_05_11_2025_rep3_ref.sorted.bed.gz FR113N_Rep2 -s FR113N_Oct_4_mRNA_induro.sorted.bed.gz FR113N_Rep3 -s W1118_induro_mRNA_2025_102.sorted.bed.gz W1118_Rep1 -s W1118_Rep2_mRNA_m6A_m5C_pseu.aligned.sorted.bed.gz W1118_Rep2 -s W1118_rep2_ref.sorted.bed.gz  W1118_Rep3 -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/FR113N_W1118_single_base_replicates_A.bed -r /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/regions/FR113N_vs_W1118.filtered.exonic.regions.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base A --threads 32 --header
modkit dmr multi -s FR113N_05_11_2025_rep2_ref.sorted.bed.gz FR113N_Rep1 -s FR113N_05_11_2025_rep3_ref.sorted.bed.gz FR113N_Rep2 -s FR113N_Oct_4_mRNA_induro.sorted.bed.gz FR113N_Rep3 -s W1118_induro_mRNA_2025_102.sorted.bed.gz W1118_Rep1 -s W1118_Rep2_mRNA_m6A_m5C_pseu.aligned.sorted.bed.gz W1118_Rep2 -s W1118_rep2_ref.sorted.bed.gz  W1118_Rep3 -o /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/FR113N_W1118_single_base_replicates_C.bed -r /mnt/george_drive/george/direct_RNAseq/modkit_pileups/FR113N_Oct_4_mRNA_Mettl3_KO_mRNA_dmr/regions/FR113N_vs_W1118.filtered.exonic.regions.bed --ref /mnt/george_drive/george/m6a/dm6.fasta --base C --threads 32 --header

# Convert a GTF to BED of exons (or genes/UTRs if you filter differently)
#gffread -E dm6.gtf -T -o- | awk 'BEGIN{FS=OFS="\t"} $3=="exon" {print $1,$4-1,$5,$9}' | bedtools sort -i - > dm6_exons.bed

# Optional: merge overlapping exons per gene
#bedtools merge -i dm6_exons.bed -c 4 -o distinct > dm6_exons_merged.bed
W1118_Rep2_mRNA_m6A_m5C_pseu.aligned.sort.bed

PAS for Mettl3-KO vs FR113N = Mettl3_FR113_PAS
PAS for Mettl3-KO vs W1118 = Mettl3_w1118_PAS
PAS for FR113N vs W1118 = FR113_w1118_PAS



library(dplyr)
library(stringr)
library(forcats)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

clean_filter <- function(df) {
  df %>% dplyr::filter(!Location %in% c("", "null", "N/A", "intron"))
}

peak_dfs <- list(
  W1118_r1_a_mRNA, W1118_r1_m_mRNA, W1118_r1_pseu_mRNA,
  W1118_r2_a_mRNA, W1118_r2_m_mRNA, W1118_r2_pseu_mRNA,
  W1118_r3_a_mRNA, W1118_r3_m_mRNA, W1118_r3_pseu_mRNA,
  FR1113N_r1_a_mRNA, FR1113N_r1_m_mRNA, FR1113N_r1_pseu_mRNA,
  FR1113N_r2_a_mRNA, FR1113N_r2_m_mRNA, FR1113N_r2_pseu_mRNA,
  FR1113N_r3_a_mRNA, FR1113N_r3_m_mRNA, FR1113N_r3_pseu_mRNA,
  Mettl3_r1_a_mRNA, Mettl3_r1_m_mRNA, Mettl3_r1_pseu_mRNA,
  Mettl3_r2_a_mRNA, Mettl3_r2_m_mRNA, Mettl3_r2_pseu_mRNA)

sample_names <- c(
  "W1118_Rep1_m6A", "W1118_Rep1_m5C", "W1118_Rep1_Pseu",
  "W1118_Rep2_m6A", "W1118_Rep2_m5C", "W1118_Rep2_Pseu",
  "W1118_Rep3_m6A", "W1118_Rep3_m5C", "W1118_Rep3_Pseu",
  "FR113N_Rep1_m6A", "FR113N_Rep1_m5C", "FR113N_Rep1_Pseu",
  "FR113N_Rep2_m6A", "FR113N_Rep2_m5C", "FR113N_Rep2_Pseu",
  "FR113N_Rep3_m6A", "FR113N_Rep3_m5C", "FR113N_Rep3_Pseu",
  "Mettl3-KO_Rep1_m6A", "Mettl3-KO_Rep1_m5C", "Mettl3-KO_Rep1_Pseu",
  "Mettl3-KO_Rep2_m6A", "Mettl3-KO_Rep2_m5C", "Mettl3-KO_Rep2_Pseu")

mod_types <- ifelse(grepl("pseu", tolower(sample_names)), "Ψ",
                    ifelse(grepl("m5c", tolower(sample_names)), "m5C", "m6A"))

annotated_peaks <- Map(function(df, sample, mod) {
  if ("mod" %in% names(df)) df <- df[, !(names(df) %in% "mod")]
  
  clean_filter(df) %>%
    mutate(Sample = sample,
           Group = case_when(grepl("^W1118", sample) ~ "W1118", grepl("^FR113N", sample) ~ "FR113N", grepl("^Mettl3-KO", sample) ~ "Mettl3-KO", TRUE ~ "Other"),
           Modification = mod, Replicate = case_when(grepl("Rep1", sample) ~ "Rep1", grepl("Rep2", sample) ~ "Rep2", grepl("Rep3", sample) ~ "Rep3", TRUE ~ "Other"))
}, peak_dfs, sample_names, mod_types) %>%
  bind_rows()

annotated_peaks$chrom <- factor(annotated_peaks$chrom, levels = sort(unique(annotated_peaks$chrom)))
annotated_peaks$Modification <- factor(annotated_peaks$Modification, levels = c("m6A", "m5C", "Ψ"))
normalize_region <- function(x) {
  y <- tolower(trimws(x))
  y <- gsub("[ _-]", "", y)
  dplyr::case_when(
    str_detect(y, "^(3'utr|utr3'|3utr|utr3|threeprimeutr|threeprime)$") ~ "3'UTR",
    str_detect(y, "^(5'utr|utr5'|5utr|utr5|fiveprimeutr|fiveprime)$") ~ "5'UTR",
    str_detect(y, "^(cds|codingsequence|cdsregion)$") ~ "CDS",
    TRUE ~ NA_character_)
}

if (!"Dataset" %in% names(annotated_peaks)) {
  annotated_peaks <- annotated_peaks %>%
    mutate(Dataset = case_when(
      grepl("^W1118", Sample) ~ "W1118",
      grepl("^FR113N", Sample) ~ "FR113N",
      grepl("^Mettl3-KO", Sample) ~ "Mettl3-KO",
      TRUE ~ "Other"))
}

ap <- annotated_peaks %>%
  mutate(Region = normalize_region(Location)) %>%
  dplyr::filter(!is.na(Region), !is.na(Gene_ID), Gene_ID != "") %>%
  mutate(Region = factor(Region, levels = c("5'UTR", "CDS", "3'UTR")))

gene_counts <- ap %>%
  group_by(Dataset, Modification, Region, Gene_ID) %>%
  summarise(Peak_Count = n(), .groups = "drop")

K <- 10
top_genes_by_ds_mod_region <- gene_counts %>%
  group_by(Dataset, Modification, Region) %>%
  arrange(desc(Peak_Count), Gene_ID) %>%
  mutate(Rank = row_number()) %>%
  dplyr::filter(Rank <= K) %>%
  ungroup()

replicate_consistent_top <- ap %>%
  group_by(Dataset, Modification, Region, Gene_ID, Replicate) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Dataset, Modification, Region, Gene_ID) %>%
  summarise(Replicate_Presence = n_distinct(Replicate),
            Peak_Count = sum(n), .groups = "drop")

region_cols <- c("5'UTR" = "#268BD2", "CDS" = "#DC322F", "3'UTR" = "#2AA198")

plot_one_dataset <- function(df, dataset_name) {
  ds_df <- df %>% filter(Dataset == dataset_name)
  
  if (nrow(ds_df) == 0) {
    return(ggplot() +
             annotate("text", x = 0, y = 0, label = paste("No data for", dataset_name)) +
             theme_void())
  }
  
  ggplot(ds_df, aes(x = fct_reorder(Gene_ID, Peak_Count), y = Peak_Count, fill = Region)) +
    geom_col(color = "black", linewidth = 0.2, width = 0.8) +
    geom_text(aes(label = Peak_Count), hjust = -0.15, size = 5) +
    coord_flip(clip = "off") +
    facet_grid(Region ~ Modification, scales = "free_y") +
    scale_fill_manual(values = region_cols) +
    scale_y_continuous(expand = expansion(mult = c(0.02, 0.12))) +
    labs(title = paste0("Top ", K, " genes per region - ", dataset_name), x = "Gene", y = "Number of modified sites", fill = "Region",  family = "Times New Roman") +
    theme_classic(base_size = 12, base_family = "Times New Roman") +
    theme(plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
          axis.title = element_text(size = 12),
          axis.text = element_text(size = 12),
          legend.position = "right", size = 12,
          legend.text=element_text(size=12, family = "Times New Roman"), 
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}
datasets <- unique(top_genes_by_ds_mod_region$Dataset)
plot_list <- setNames(lapply(datasets, function(ds)
  plot_one_dataset(top_genes_by_ds_mod_region, ds)), datasets)

combined_fig <- wrap_plots(plot_list, ncol = 3) +
  plot_annotation(
    title = "Top genes per transcript region by dataset and modification",
    caption = paste0("Top ", K, " genes within each Region × Modification per dataset."),
    theme = theme_classic(base_size = 12) +
      theme(plot.title = element_text(face = "bold", hjust = 0.5)))


top_shared_genes <- (plot_list$W1118 | plot_list$FR113N| plot_list$`Mettl3-KO`) +
  plot_annotation(title = " ", theme = theme(plot.title = element_text(size = 18, family = "Times New Roman")), tag_levels = 'A') &
  theme(text = element_text(size = 18),  axis.title = element_text(size = 18), axis.text = element_text(size = 18), strip.text = element_text(size = 18), legend.text = element_text(size = 18), legend.title = element_text(size = 18))

Head_indepen_lia_lab = read.delim("WT_head_lai/wt_vs_ime4_head_minus_cond1.tsv") %>% dplyr::filter(Gene_ID != "null") %>% select(-chr, -start, -end, -strand, -PeakID, -fold_enrichment, -Transcript_IDs, -Transcript_Biotype, -X)

Head_depen_lia_lab = read.delim("WT_head_lai/wt_vs_ime4_head_minus_cond2.tsv") %>% dplyr::filter(Gene_ID != "null") %>% select(-chr, -start, -end, -strand, -PeakID, -fold_enrichment, -Transcript_IDs, -Transcript_Biotype)


neuron_biased_cleary_lab = read.delim("WT_head_lai/meRIP_neuron_biased_brain_2_minus_peaks.tsv") %>% dplyr::filter(Gene_ID != "null") %>% select(-chr, -start, -end, -strand, -name, -fold_enrichment, -Transcript_IDs, -Transcript_Biotype, -length, -abs_summit, -pileup, -X.log10.pvalue., -X.log10.qvalue., -Location)


neuroblast_cleary_lab = read.delim("WT_head_lai/meRIP_neuroblast_biased_brain_2_minus_peaks.tsv") %>% dplyr::filter(Gene_ID != "null") %>% select(-chr, -start, -end, -strand, -name, -fold_enrichment, -Transcript_IDs, -Transcript_Biotype, -length, -abs_summit, -pileup, -X.log10.pvalue., -X.log10.qvalue., -Location)

Mettl3null_brain_cleary_lab = read.delim("WT_head_lai/meRIP_Mettl3null_brain_1_minus_peaks.tsv") %>% dplyr::filter(Gene_ID != "null") %>% select(-chr, -start, -end, -strand, -name, -fold_enrichment, -Transcript_IDs, -Transcript_Biotype, -length, -abs_summit, -pileup, -X.log10.pvalue., -X.log10.qvalue., -Location)

boninie_brain = read.delim("WT_head_lai/Mettl3_dependent_genes.txt", header =T, sep ="\t") %>% dplyr::filter(Gene_ID != "null")

library(dplyr)
library(purrr)
library(tibble)

# External datasets: gene vectors

lai_indep_genes <- Head_indepen_lia_lab %>%
  pull(Gene_ID) %>%
  unique()

lai_dep_genes <- Head_depen_lia_lab %>%
  pull(Gene_ID) %>%
  unique()

cleary_neuron_genes <- neuron_biased_cleary_lab %>%
  pull(Gene_ID) %>%
  unique()

cleary_neuroblast_genes <- neuroblast_cleary_lab %>%
  pull(Gene_ID) %>%
  unique()

cleary_Mettl3null_genes <- Mettl3null_brain_cleary_lab %>%
  pull(Gene_ID) %>%
  unique()

bonini_genes <- boninie_brain %>%
  pull(Gene_ID) %>%
  unique()

# Put them in a named list for convenience
external_sets <- list(Lai_head_independent = lai_indep_genes,
                      Lai_head_dependent = lai_dep_genes,
                      Cleary_neuron_biased = cleary_neuron_genes,
                      Cleary_neuroblast_biased = cleary_neuroblast_genes,
                      Cleary_Mettl3null_brain = cleary_Mettl3null_genes,
                      Bonini_brain = bonini_genes)

# the m6A genes by genotype

get_genotype_m6A_genes <- function(geno_pattern) {
  peak_files[grepl(paste0("^", geno_pattern, ".*m6A$"), names(peak_files))] %>%
    map(~ .x$Gene_ID) %>%
    unlist() %>%
    unique()
}

geno_m6A <- list(
  w1118 = get_genotype_m6A_genes("w1118"),
  FR113N = get_genotype_m6A_genes("FR113"),
  Mettl3KO = get_genotype_m6A_genes("Mettl3"))

# Quick sanity check
sapply(geno_m6A, length)


#  Overlap by genotype & external dataset

overlap_genotype_full <- expand.grid(
  genotype = names(geno_m6A),
  dataset  = names(external_sets),
  stringsAsFactors = FALSE) %>%
  as_tibble() %>%
  rowwise() %>%
  mutate(shared_genes = list(intersect(geno_m6A[[genotype]], external_sets[[dataset]])),
         n_shared = length(shared_genes),
         geno_total = length(geno_m6A[[genotype]]),
         dataset_total = length(external_sets[[dataset]]),
         frac_of_genotype = ifelse(geno_total > 0, n_shared / geno_total, NA_real_),
         frac_of_dataset = ifelse(dataset_total > 0, n_shared / dataset_total, NA_real_)) %>%
  ungroup()

# Numeric summary (no gene lists column)
overlap_genotype_summary_full <- overlap_genotype_full %>%
  select(genotype, dataset, n_shared, geno_total, dataset_total,
         frac_of_genotype, frac_of_dataset)

overlap_genotype_summary_full %>% print(width = Inf)


library(dplyr)
library(tidyr)
library(tibble)

# starting from overlap_genotype_summary_full
# columns: genotype, dataset, n_shared, geno_total, dataset_total, frac_of_genotype, frac_of_dataset

heat_mat <- overlap_genotype_summary_full %>%
  select(dataset, genotype, frac_of_dataset) %>%
  pivot_wider(names_from  = genotype,
              values_from = frac_of_dataset) %>%
  column_to_rownames("dataset") %>%
  as.matrix()

heat_mat

library(pheatmap)
library(RColorBrewer)

# continuous palette from low (white) to high (blue)
heat_colors <- colorRampPalette(c("beige", "#1f78b4"))(100)

heatmap_overlap <- pheatmap(heat_mat, color = heat_colors, cluster_rows = FALSE, cluster_cols = FALSE, scale = "none", show_rownames = TRUE, show_colnames = TRUE, main = "Fraction of external m6A datasets recovered per genotype", legend = TRUE, legend_breaks = c(0, 0.5, 1), legend_labels = c("0", "0.5", "1.0"), border_color = "grey80", annotation_row = annotation_row, treeheight_row = 0, treeheight_col = 0, fontsize = 10)


##psi overlap
library(dplyr)
library(purrr)
library(tibble)

# Helper: get union of Ψ Gene_IDs per genotype
# Match any of "Ψ", "psi", "pseu" in the sample name, just to be robust
get_genotype_Psi_genes <- function(geno_pattern) {
  peak_files[
    grepl(paste0("^", geno_pattern), names(peak_files)) &
      grepl("Ψ|psi|pseu", names(peak_files))
  ] %>%
    map(~ .x$Gene_ID) %>%
    unlist() %>%
    unique()
}

geno_Psi <- list(w1118 = get_genotype_Psi_genes("w1118"), FR113N = get_genotype_Psi_genes("FR113"), Mettl3KO = get_genotype_Psi_genes("Mettl3"))

# Quick check: how many Ψ genes per genotype
sapply(geno_Psi, length)


Dan_trace_Psi <- read.delim("Dan_trace_Psi.txt", header = TRUE, sep = "\t") %>% dplyr::filter(Gene_ID != "null")

dan_psi_genes <- Dan_trace_Psi %>%
  pull(Gene_ID) %>%
  unique()

length(dan_psi_genes)  # how many Ψ genes in Dan et al. mRNA dataset

# genotype × Dan Psi (mRNA) overlap
overlap_Psi_genotype <- tibble(
  genotype = names(geno_Psi)) %>%
  rowwise() %>%
  mutate(shared_genes = list(intersect(geno_Psi[[genotype]], dan_psi_genes)),
         n_shared = length(shared_genes),
         geno_total = length(geno_Psi[[genotype]]),
         dataset_total = length(dan_psi_genes),
         frac_of_genotype = ifelse(geno_total > 0, n_shared / geno_total, NA_real_),
         frac_of_dataset = ifelse(dataset_total > 0, n_shared / dataset_total, NA_real_)) %>%
  ungroup()

overlap_Psi_genotype


# Go for exclusive genes 
install.packages("BiocManager")
BiocManager::install(c("clusterProfiler", "org.Dm.eg.db"))
library(clusterProfiler)
library(org.Dm.eg.db)


library(dplyr)
library(purrr)
library(tibble)

## m6A: exclusive genes per genotype

# union of all external m6A gene sets (Lai + Cleary + Bonini)
external_m6A_union <- external_sets %>%
  unlist() %>%
  unique()

exclusive_m6A <- lapply(names(geno_m6A), function(g) {
  setdiff(geno_m6A[[g]], external_m6A_union)
}) %>%
  setNames(names(geno_m6A))

exclusive_m6A

# Optional: counts
sapply(exclusive_m6A, length)


## Psi: exclusive genes per genotype

exclusive_Psi <- lapply(names(geno_Psi), function(g) {
  setdiff(geno_Psi[[g]], dan_psi_genes)
}) %>%
  setNames(names(geno_Psi))

exclusive_Psi

# Optional: counts
sapply(exclusive_Psi, length)

# all my m6A genes (union of all genotypes)
universe_m6A <- geno_m6A %>%
  unlist() %>%
  unique()

# all my Psi genes (union of all genotypes)
universe_Psi <- geno_Psi %>%
  unlist() %>%
  unique()

run_go_for_gene_set <- function(genes, universe, keytype = "FLYBASE", ont = "BP", p_cutoff = 0.05, q_cutoff = 0.2) {
  enrichGO(gene = genes, universe = universe, OrgDb = org.Dm.eg.db, keyType = keytype,
           ont = ont,             # "BP", "MF", or "CC"
           pAdjustMethod = "BH", pvalueCutoff = p_cutoff, qvalueCutoff = q_cutoff, readable = TRUE)
}


## GO for m6A-exclusive genes 

go_m6A_exclusive <- lapply(names(exclusive_m6A), function(g) {
  genes <- exclusive_m6A[[g]]
  
  message("Running GO for m6A-exclusive genes in ", g,
          " (n = ", length(genes), ")")
  
  ego <- run_go_for_gene_set(genes = genes, universe = universe_m6A, keytype = "SYMBOL")
  
  ego@result %>%
    as_tibble() %>%
    mutate(genotype = g, mark = "m6A")
}) %>%
  bind_rows()

go_m6A_exclusive %>% dplyr::glimpse()

#write.table(go_m6A_exclusive, file = "GO_m6A_exclusive_by_genotype.txt", sep = "\t", quote = FALSE, row.names = FALSE)



# GO for Ψ-exclusive genes per genotype

## GO for Psi-exclusive genes

go_Psi_exclusive <- lapply(names(exclusive_Psi), function(g) {
  genes <- exclusive_Psi[[g]]
  
  message("Running GO for Psi-exclusive genes in ", g,
          " (n = ", length(genes), ")")
  
  ego <- run_go_for_gene_set(genes = genes, universe = universe_Psi, keytype = "SYMBOL")
  
  ego@result %>%
    as_tibble() %>%
    mutate(genotype = g, mark = "Psi")
}) %>%
  bind_rows()

go_Psi_exclusive %>% dplyr::glimpse()

#write.table(go_Psi_exclusive, file = "GO_Psi_exclusive_by_genotype.txt", sep = "\t", quote = FALSE, row.names = FALSE)



library(dplyr)
library(ggplot2)
library(forcats)

TOP_N_GO <- 10   # change to 15, 20, etc if you prefer

# keep only significant-ish terms (optional filter)
go_m6A_top <- go_m6A_exclusive %>%
  filter(!is.na(p.adjust)) %>%
  group_by(genotype) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = TOP_N_GO) %>%
  ungroup() %>%
  mutate(Description = fct_reorder(Description, -log10(p.adjust)))


TOP_N_GO <- 10  # top GO terms per genotype

##Plot for m6A-exclusive

go_m6A_top <- go_m6A_exclusive_manual %>%
  filter(!is.na(p.adjust)) %>%
  group_by(genotype) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = TOP_N_GO) %>%
  ungroup() %>%
  mutate(Description = fct_reorder(Description, -log10(p.adjust)))

p_GO_m6A <- ggplot(go_m6A_top, aes(x = Description, y = -log10(p.adjust))) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ genotype, ncol = 1, scales = "free_y") +
  labs(title = "GO BP enrichment for m6A-exclusive genes", x = "GO Biological Process", y = expression(-log[10]("adj p"))) +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title.y = element_text(size = 9),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 12, family = "Times New Roman"))


ggsave(filename = "WT_head_lai/GO_m6A_exclusive_by_genotype.pdf", plot = p_GO_m6A, units = "in", width = 10, height = 9, device = cairo_pdf, bg = "white")


## Plot for Psi-exclusive

go_Psi_top <- go_Psi_exclusive_manual %>%
  filter(!is.na(p.adjust)) %>%
  group_by(genotype) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  slice_head(n = TOP_N_GO) %>%
  ungroup() %>%
  mutate(
    Description = fct_reorder(Description, -log10(p.adjust)))

p_GO_Psi <- ggplot(go_Psi_top, aes(x = Description, y = -log10(p.adjust))) +
  geom_col() +
  coord_flip() +
  facet_wrap(~ genotype, ncol = 1, scales = "free_y") +
  labs(title = "GO BP enrichment for Psi-exclusive genes", x = "GO Biological Process", y = expression(-log[10](adj~p))) +
  theme_classic(base_size = 12, base_family = "Times New Roman") +
  theme(
    strip.text = element_text(face = "bold"),
    axis.title.y = element_text(size = 9),
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 12),
    legend.position = "right",
    legend.text = element_text(size = 12, family = "Times New Roman"))


ggsave(filename = "WT_head_lai/GO_Psi_exclusive_by_genotype.pdf", plot = p_GO_Psi, units = "in", width = 10, height = 9, device = cairo_pdf, bg = "white")


# GO for m5C

library(dplyr)
library(purrr)
library(tibble)

# helper: union of m5C Gene_IDs per genotype
get_genotype_m5C_genes <- function(geno_pattern) {
  peak_files[
    grepl(paste0("^", geno_pattern), names(peak_files)) &
      grepl("m5C", names(peak_files))
  ] %>%
    map(~ .x$Gene_ID) %>%
    unlist() %>%
    unique()
}

geno_m5C <- list(
  w1118 = get_genotype_m5C_genes("w1118"),
  FR113N = get_genotype_m5C_genes("FR113"),
  Mettl3KO = get_genotype_m5C_genes("Mettl3"))

# quick check
sapply(geno_m5C, length)

universe_m5C <- geno_m5C %>%
  unlist() %>%
  unique()
length(universe_m5C)


library(AnnotationDbi)
library(org.Dm.eg.db)
library(GO.db)

KEYTYPE <- "SYMBOL"

build_go_mapping <- function(universe_genes, keytype = "SYMBOL") {
  annot <- AnnotationDbi::select(
    org.Dm.eg.db,
    keys = universe_genes,
    columns = c("GO", "ONTOLOGY"),
    keytype = keytype)
  
  annot_bp <- annot %>%
    filter(ONTOLOGY == "BP", !is.na(GO))
  
  go2genes <- split(annot_bp[[keytype]], annot_bp$GO) %>%
    lapply(unique)
  
  go2genes
}

go2genes_m5C <- build_go_mapping(universe_m5C, keytype = KEYTYPE)


go_enrich_manual <- function(gene_set, universe_genes, go2genes, keytype = "SYMBOL") {
  N <- length(universe_genes)
  M <- length(gene_set)
  
  res <- purrr::imap_dfr(go2genes, function(term_genes, go_id) {
    term_genes <- intersect(term_genes, universe_genes)
    n_term     <- length(term_genes)
    if (n_term == 0) return(NULL)
    
    k <- length(intersect(gene_set, term_genes))
    if (k == 0) return(NULL)
    
    a <- k
    b <- M - k
    c <- n_term - k
    d <- N - n_term - b
    
    mat <- matrix(c(a, b, c, d), nrow = 2)
    ft  <- fisher.test(mat, alternative = "greater")
    
    tibble(GO_ID = go_id, k_in_set = a, M_set_size = M, n_in_term = n_term, N_universe_size  = N,
           pvalue = ft$p.value)
  })
  
  if (nrow(res) == 0) return(NULL)
  
  res <- res %>%
    mutate(p.adjust = p.adjust(pvalue, method = "BH"))
  
  # use GOTERM + AnnotationDbi::Term() to get human-readable names
  term_vec <- vapply(res$GO_ID, function(id) {
    obj <- GO.db::GOTERM[[id]]
    if (is.null(obj)) {
      NA_character_
    } else {
      AnnotationDbi::Term(obj)
    }
  }, character(1))
  
  res$Description <- term_vec
  
  res
}

go_m5C <- purrr::imap_dfr(geno_m5C, function(gene_set, geno) {
  message("Running GO for m5C genes in ", geno,
          " (n = ", length(gene_set), ")")
  
  res <- go_enrich_manual(gene_set = gene_set, universe_genes = universe_m5C, go2genes = go2genes_m5C, keytype = KEYTYPE)
  
  if (is.null(res)) return(NULL)
  
  res %>%
    mutate(genotype = geno, mark = "m5C")
})


# save raw table if you want
write.table(go_m5C, file  = "WT_head_lai/GO_m5C_by_genotype_manual.txt", sep = "\t", quote = FALSE, row.names = FALSE)

