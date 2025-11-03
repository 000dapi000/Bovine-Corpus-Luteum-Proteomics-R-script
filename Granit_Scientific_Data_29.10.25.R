library(clValid)
library(openxlsx)  # for saving Excel files

# --- Set working directory ---
setwd("C:/Users/ga53hil/Desktop/Granit_proteomics/28.08.25_result")

# --- Load and preprocess data ---
data <- read.table("RAW_P337_02B_proteinGroups.txt", 
                   header = TRUE, sep = "\t", quote = "", check.names = FALSE)

# --- Filter out unwanted rows ---
data_filtered <- data[
  data$`Only identified by site` != "+" &
    data$Reverse != "+" &
    data$`Potential contaminant` != "+",
]

# --- Reset row names ---
row.names(data_filtered) <- NULL

# --- Save filtered data (TXT) ---
write.table(data_filtered,
            "P337_02B_proteinGroups_filtered.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save filtered data (Excel) ---
write.xlsx(data_filtered,
           "P337_02B_proteinGroups_filtered.xlsx",
           overwrite = TRUE)


################################################################################
################################################################################
################################################################################

# --- Keep proteins with Peptides > 1 ---
# Coerce safely in case it's read as character
data_filtered$Peptides_num <- suppressWarnings(as.numeric(data_filtered$Peptides))

data_filtered2 <- subset(
  data_filtered,
  !is.na(Peptides_num) & Peptides_num > 1
)

# Drop helper column
data_filtered2$Peptides_num <- NULL

# --- Save result (TXT) ---
write.table(data_filtered2,
            "P337_02B_proteinGroups_filtered_pepGT1.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save result (Excel) ---
library(openxlsx)
write.xlsx(data_filtered2,
           "P337_02B_proteinGroups_filtered_pepGT1.xlsx",
           overwrite = TRUE)


################################################################################
################################################################################
################################################################################

# --- Identify LFQ intensity columns (original names) ---
lfq_cols <- grep("^LFQ intensity\\s", colnames(data_filtered2), value = TRUE)

# --- Copy before transforming ---
log2_data_filtered2 <- data_filtered2

# --- Robust log2 transform of LFQ columns ---
log2_data_filtered2[lfq_cols] <- lapply(data_filtered2[lfq_cols], function(x) {
  # 1) force to character, strip thousands separators, then to numeric
  x_num <- suppressWarnings(as.numeric(gsub(",", "", as.character(x))))
  # 2) treat zeros/negatives as missing (will be imputed later)
  x_num[x_num <= 0] <- NA_real_
  # 3) log2 and clean non-finite
  y <- log2(x_num)
  y[!is.finite(y)] <- NA_real_
  y
})

# --- Save (TXT) ---
write.table(log2_data_filtered2,
            "P337_02B_proteinGroups_filtered_pepGT1_log2LFQ.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save (Excel) ---
library(openxlsx)
write.xlsx(log2_data_filtered2,
           "P337_02B_proteinGroups_filtered_pepGT1_log2LFQ.xlsx",
           overwrite = TRUE)

################################################################################
################################################################################
################################################################################

# Mapping from Roman numerals to numbers
roman_to_num <- c(
  "I" = 1, "II" = 2, "III" = 3, "IV" = 4, "V" = 5,
  "VI" = 6, "VII" = 7, "VIII" = 8, "IX" = 9, "X" = 10
)

# Loop through LFQ intensity columns and rename
new_colnames <- colnames(log2_data_filtered2)
for (i in seq_along(new_colnames)) {
  if (grepl("^LFQ intensity", new_colnames[i])) {
    # Extract the Roman numeral and replicate number
    parts <- unlist(strsplit(gsub("^LFQ intensity ", "", new_colnames[i]), "-"))
    roman <- trimws(parts[1])
    rep <- trimws(parts[2])
    # Map to new name format
    new_colnames[i] <- paste0("T", roman_to_num[roman], "_", rep)
  }
}

# Apply new column names
colnames(log2_data_filtered2) <- new_colnames

################################################################################
################################################################################
################################################################################

# --- Identify LFQ (T#_#) columns ---
lfq_cols <- grep("^T\\d+_\\d+$", colnames(log2_data_filtered2), value = TRUE)

# --- Clean NaN/Inf -> NA in LFQ columns ---
log2_data_filtered2[lfq_cols] <- lapply(log2_data_filtered2[lfq_cols], function(x) {
  x[!is.finite(x)] <- NA_real_
  x
})

# --- Build group list (T1..T10 detected from column names) ---
groups <- sort(unique(sub("^T(\\d+)_.*", "\\1", lfq_cols)))

# --- Keep rows where at least one group has >=70% valid values ---
passes_any_group <- Reduce(`|`, lapply(groups, function(g) {
  cols_g <- grep(paste0("^T", g, "_\\d+$"), colnames(log2_data_filtered2), value = TRUE)
  if (length(cols_g) == 0) return(rep(FALSE, nrow(log2_data_filtered2)))
  thr <- ceiling(0.70 * length(cols_g))   # e.g., 0.70 * 8 = 5.6 -> 6
  rowSums(!is.na(log2_data_filtered2[cols_g])) >= thr
}))

filtered_70 <- log2_data_filtered2[passes_any_group, , drop = FALSE]
row.names(filtered_70) <- NULL

# --- Save result (TXT) ---
write.table(filtered_70,
            "P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# --- Save result (Excel) ---
library(openxlsx)
write.xlsx(filtered_70,
           "P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70.xlsx",
           overwrite = TRUE)

################################################################################
################################################################################
################################################################################

library(dplyr)
library(ggplot2)

groups <- list(
  T1  = paste0("Peptides I-",     1:8),
  T2  = paste0("Peptides II-",    1:8),
  T3  = paste0("Peptides III-",   1:8),
  T4  = paste0("Peptides IV-",    1:8),
  T5  = paste0("Peptides V-",     1:8),
  T6  = paste0("Peptides VI-",    1:8),
  T7  = paste0("Peptides VII-",   1:8),
  T8  = paste0("Peptides VIII-",  1:8),
  T9  = paste0("Peptides IX-",    1:8),
  T10 = paste0("Peptides X-",     1:8)
)
all_cols <- unlist(groups, use.names = FALSE)
filtered_70[all_cols] <- lapply(filtered_70[all_cols], function(x) suppressWarnings(as.numeric(x)))

replicate_sums_mat <- sapply(groups, function(cols) {
  colSums(filtered_70[, cols, drop = FALSE], na.rm = TRUE)
})
df <- data.frame(
  Condition  = rep(names(groups), each = 8),
  Replicate  = rep(1:8, times = length(groups)),
  PeptideSum = as.vector(replicate_sums_mat)
)
df$Condition <- factor(df$Condition, levels = paste0("T", 1:10))

summary_df <- df %>%
  group_by(Condition) %>%
  summarise(
    mean = mean(PeptideSum, na.rm = TRUE),
    sd   = sd(PeptideSum,   na.rm = TRUE),
    .groups = "drop"
  )

y_max <- max(
  summary_df$mean + summary_df$sd,
  df$PeptideSum
) * 1.10

plot_peptides <- ggplot(summary_df, aes(x = Condition, y = mean)) +
  geom_bar(stat = "identity", fill = "#0065bd", width = 0.85) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.3, linewidth = 1, color = "black") +
  geom_point(
    data = df, aes(x = Condition, y = PeptideSum),
    position = position_jitter(width = 0.10, height = 0),
    size = 3.5, color = "black", fill = "white", shape = 21, stroke = 1.2, alpha = 1
  ) +
  labs(
    x = "Time points",
    y = "Number of peptides",
    title = "Total Peptides per Condition"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
  theme_classic(base_size = 20) +
  theme(
    axis.line        = element_line(color = "black", linewidth = 1.2),
    axis.ticks       = element_line(color = "black", linewidth = 1),
    axis.text.x      = element_text(size = 22, color = "black", face = "bold"),
    axis.text.y      = element_text(size = 22, color = "black", face = "bold"),
    axis.title.x     = element_text(size = 24, color = "black", face = "bold", margin = margin(t = 18)),
    axis.title.y     = element_text(size = 24, color = "black", face = "bold", margin = margin(r = 18)),
    plot.title       = element_text(size = 26, color = "black", face = "bold", hjust = 0.5, margin = margin(b = 22)),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(plot_peptides)

# Export PNG, white background, 300 dpi, size 10x7 inches
ggsave("peptide_totals_T1_T10_opencircle.png", plot_peptides, width = 10, height = 7, dpi = 300, bg = "white")

################################################################################
################################################################################
################################################################################

library(dplyr)
library(ggplot2)

# 1. Define sample columns in T1_1 to T10_8 order
sample_cols <- c(
  paste0("T1_", 1:8), paste0("T2_", 1:8), paste0("T3_", 1:8), paste0("T4_", 1:8),
  paste0("T5_", 1:8), paste0("T6_", 1:8), paste0("T7_", 1:8), paste0("T8_", 1:8),
  paste0("T9_", 1:8), paste0("T10_", 1:8)
)

# 2. Make protein presence matrix (TRUE = present, FALSE/NA = missing)
protein_matrix <- !is.na(as.matrix(filtered_70[sample_cols]))

# 3. Protein count per sample (number of TRUEs per column)
protein_counts_per_sample <- colSums(protein_matrix)
sample_names <- sample_cols

# 4. Calculate cumulative and shared trends
cumulative_protein <- sapply(1:length(sample_names), function(i) {
  sum(rowSums(protein_matrix[, 1:i, drop = FALSE]) > 0)
})

shared_protein <- sapply(1:length(sample_names), function(i) {
  sum(rowSums(protein_matrix[, 1:i, drop = FALSE]) == i)
})

# 5. Build dataframe for plotting
df <- data.frame(
  Sample = factor(sample_names, levels = sample_names),
  ProteinCount = protein_counts_per_sample,
  Cumulative = cumulative_protein,
  Shared = shared_protein
)

# 6. Set y-axis upper limit for space above trend lines/bars
y_max <- max(df$ProteinCount, df$Cumulative, df$Shared, na.rm = TRUE) * 1.10

# 7. Set plot width automatically: e.g., 0.25 inch per sample (20 samples = 5", 80 = 20")
n_samples <- length(sample_names)
plot_width <- max(12, n_samples * 0.25)  # You can adjust 0.25 for your monitor or journal

# 8. Plot with all sample names
p <- ggplot(df, aes(x = Sample)) +
  geom_bar(aes(y = ProteinCount), stat = "identity", fill = "#0065bd", width = 0.85) +
  geom_line(aes(y = Cumulative, group = 1, color = "Cumulative"), size = 1.2) +
  geom_point(aes(y = Cumulative, color = "Cumulative"), size = 2.2) +
  geom_line(aes(y = Shared, group = 1, color = "Shared"), size = 1.2) +
  geom_point(aes(y = Shared, color = "Shared"), size = 2.2) +
  scale_color_manual(values = c("Cumulative" = "green4", "Shared" = "orange2")) +
  labs(
    y = "Number of proteins (ID)",
    x = "Sample (T1_1 to T10_8)",
    title = "Protein IDs per Sample with Shared and Cumulative Trends"
  ) +
  scale_x_discrete(labels = sample_names) +   # Show every sample name!
  scale_y_continuous(expand = c(0, 0), limits = c(0, y_max)) +
  theme_classic(base_size = 22) +
  theme(
    legend.title = element_blank(),
    legend.position = c(0.15, 0.15),
    axis.text.x = element_text(angle = 70, vjust = 1, hjust = 1, size = 14), # adjust text size for many labels
    axis.text.y = element_text(size = 22, color = "black", face = "bold"),
    axis.title.x = element_text(size = 24, color = "black", face = "bold", margin = margin(t = 18)),
    axis.title.y = element_text(size = 24, color = "black", face = "bold", margin = margin(r = 18)),
    plot.title   = element_text(size = 26, color = "black", face = "bold", hjust = 0.5, margin = margin(b = 22)),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(p)

# 9. Export PNG, white background, 300 dpi, auto width, height fixed (e.g., 7 inches)
ggsave("protein_IDs_per_sample_trends_fullnames.png", p, width = plot_width, height = 7, dpi = 300, bg = "white")

################################################################################
################################################################################
################################################################################

library(dplyr)
library(ggplot2)

# 1. Define sample columns in T1_1 to T10_8 order (proteins)
sample_cols <- c(
  paste0("T1_", 1:8), paste0("T2_", 1:8), paste0("T3_", 1:8), paste0("T4_", 1:8),
  paste0("T5_", 1:8), paste0("T6_", 1:8), paste0("T7_", 1:8), paste0("T8_", 1:8),
  paste0("T9_", 1:8), paste0("T10_", 1:8)
)

# 2. Convert to numeric (safeguard)
filtered_70[sample_cols] <- lapply(filtered_70[sample_cols], function(x) suppressWarnings(as.numeric(x)))

# 3. Count proteins per sample (LFQ > 0 = present)
protein_counts_mat <- sapply(split(sample_cols, rep(1:10, each=8)), function(cols) {
  sapply(cols, function(col) sum(filtered_70[[col]] > 0, na.rm = TRUE))
})

df_protein <- data.frame(
  Condition  = rep(paste0("T", 1:10), each = 8),
  Replicate  = rep(1:8, times = 10),
  ProteinCount = as.vector(protein_counts_mat)
)
df_protein$Condition <- factor(df_protein$Condition, levels = paste0("T", 1:10))

# 4. Calculate mean + SD per condition
summary_protein <- df_protein %>%
  group_by(Condition) %>%
  summarise(
    mean = mean(ProteinCount, na.rm = TRUE),
    sd   = sd(ProteinCount,   na.rm = TRUE),
    .groups = "drop"
  )

# 5. Set y-axis upper limit for space above error bar/dots
y_max_protein <- max(
  summary_protein$mean + summary_protein$sd,
  df_protein$ProteinCount
) * 1.10

# 6. Plot (open circles for samples)
plot_protein <- ggplot(summary_protein, aes(x = Condition, y = mean)) +
  geom_bar(stat = "identity", fill = "#0065bd", width = 0.85) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.3, linewidth = 1, color = "black") +
  geom_point(
    data = df_protein, aes(x = Condition, y = ProteinCount),
    position = position_jitter(width = 0.10, height = 0),
    size = 3.5, color = "black", fill = "white", shape = 21, stroke = 1.2, alpha = 1
  ) +
  labs(
    x = "Time points",
    y = "Number of identified proteins",
    title = "Total Proteins per Condition"
  ) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, y_max_protein)) +
  theme_classic(base_size = 20) +
  theme(
    axis.line        = element_line(color = "black", linewidth = 1.2),
    axis.ticks       = element_line(color = "black", linewidth = 1),
    axis.text.x      = element_text(size = 22, color = "black", face = "bold"),
    axis.text.y      = element_text(size = 22, color = "black", face = "bold"),
    axis.title.x     = element_text(size = 24, color = "black", face = "bold", margin = margin(t = 18)),
    axis.title.y     = element_text(size = 24, color = "black", face = "bold", margin = margin(r = 18)),
    plot.title       = element_text(size = 26, color = "black", face = "bold", hjust = 0.5, margin = margin(b = 22)),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

print(plot_protein)

# 7. Export PNG, white background, 300 dpi, size 10x7 inches
ggsave("protein_totals_T1_T10_opencircle.png", plot_protein, width = 10, height = 7, dpi = 300, bg = "white")

################################################################################
################################################################################
################################################################################

################################################################################
# Perseus-style imputation (Total matrix) with fixed seed for reproducibility
# Normal( μ - 1.8·σ , (0.3·σ)^2 )
################################################################################

# 1) Select LFQ columns (T#_#)
expr_cols <- grep("^T\\d+_\\d+$", colnames(filtered_70), value = TRUE)
if (!length(expr_cols)) stop("No LFQ columns (T#_#) found for imputation.")

# 2) Build numeric matrix; treat non-finite as missing
X <- as.matrix(filtered_70[expr_cols])
storage.mode(X) <- "double"
X[!is.finite(X)] <- NA_real_

# 3) Compute μ and σ from ALL valid values across selected columns (Total matrix)
vals <- X[!is.na(X)]
if (length(vals) < 2) stop("Not enough valid values to estimate μ and σ for imputation.")
mu  <- mean(vals)
sig <- sd(vals)
if (!is.finite(sig) || sig == 0) stop("σ is not finite or zero; cannot perform Perseus-style imputation.")

# 4) Define imputation Gaussian (Width = 0.3, Down shift = 1.8)
mu_imp <- mu - 1.8 * sig
sd_imp <- 0.3 * sig

# 5) Impute all missing values with fixed seed
set.seed(1)  # change to any integer for different reproducible draws
miss <- which(is.na(X))
if (length(miss)) {
  X[miss] <- rnorm(length(miss), mean = mu_imp, sd = sd_imp)
}

# 6) Replace and save
imputed_data <- filtered_70
imputed_data[expr_cols] <- X

# Save TXT
write.table(imputed_data,
            "P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70_imputed_fixedseed.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# Save Excel
library(openxlsx)
write.xlsx(imputed_data,
           "P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70_imputed_fixedseed.xlsx",
           overwrite = TRUE)

################################################################################
################################################################################
################################################################################

################################################################################
# Quantile normalization after imputation (base R) + save as Excel
################################################################################

# 1) Select LFQ columns
expr_cols <- grep("^T\\d+_\\d+$", colnames(imputed_data), value = TRUE)

# 2) Quantile normalization function (base R)
qn <- function(m) {
  r <- apply(m, 2, rank, ties.method = "min")
  s <- apply(m, 2, sort)
  m_avg <- rowMeans(s)
  m[] <- m_avg[r]
  m
}

# 3) Apply quantile normalization
quantile_norm_data <- imputed_data
quantile_norm_data[expr_cols] <- qn(as.matrix(imputed_data[expr_cols]))

# 4) Save as TXT
write.table(quantile_norm_data,
            "P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70_imputed_fixedseed_quantilenorm.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

# 5) Save as Excel
library(openxlsx)
write.xlsx(quantile_norm_data,
           "P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70_imputed_fixedseed_quantilenorm.xlsx",
           overwrite = TRUE)

################################################################################
################################################################################
################################################################################

library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)

sample_cols <- c(
  paste0("T1_", 1:8), paste0("T2_", 1:8), paste0("T3_", 1:8), paste0("T4_", 1:8),
  paste0("T5_", 1:8), paste0("T6_", 1:8), paste0("T7_", 1:8), paste0("T8_", 1:8),
  paste0("T9_", 1:8), paste0("T10_", 1:8)
)

# Color palette for T1–T10 (same as your sPLS-DA)
group_colors <- setNames(
  c("#66c2a5", "#ffd92f", "#8da0cb", "#fc8d62", "#a6cee3",
    "#1f78b4", "#ffb347", "#b3de69", "#bdbdbd", "#bc80bd"),
  paste0("T", 1:10)
)

# Prepare long-form data for both
lfq_long_filtered <- filtered_70 %>%
  as.data.frame() %>%
  dplyr:: select(all_of(sample_cols)) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "LFQ") %>%
  mutate(Group = sub("_.*", "", Sample))
lfq_long_filtered$Group <- factor(lfq_long_filtered$Group, levels = paste0("T", 1:10))

lfq_long_norm <- quantile_norm_data %>%
  as.data.frame() %>%
  dplyr::select(all_of(sample_cols)) %>%
  pivot_longer(cols = everything(), names_to = "Sample", values_to = "LFQ") %>%
  mutate(Group = sub("_.*", "", Sample))
lfq_long_norm$Group <- factor(lfq_long_norm$Group, levels = paste0("T", 1:10))

# y limits across both datasets
all_min <- min(c(lfq_long_filtered$LFQ, lfq_long_norm$LFQ), na.rm = TRUE)
all_max <- max(c(lfq_long_filtered$LFQ, lfq_long_norm$LFQ), na.rm = TRUE)
my_y_limits <- c(all_min, all_max)

# Official boxplot style function with custom fill color
make_boxplot <- function(df, title_text, ylimits) {
  ggplot(df, aes(x = Group, y = LFQ, fill = Group)) +
    geom_boxplot(
      stat = "boxplot",
      position = "dodge2",
      outlier.colour = "#0065bd",
      outlier.color = "#0065bd",
      outlier.fill = "#0065bd",
      outlier.shape = 19,
      outlier.size = 1.5,
      outlier.stroke = 0.5,
      outlier.alpha = 0.7,
      notch = FALSE,
      notchwidth = 0.5,
      staplewidth = 0,
      varwidth = FALSE,
      na.rm = FALSE,
      orientation = NA,
      show.legend = NA,
      inherit.aes = TRUE,
      coef = 1.5
    ) +
    scale_fill_manual(values = group_colors) +   # <-- Use your custom colors!
    labs(
      title = title_text,
      x = "Group",
      y = "LFQ Intensity"
    ) +
    theme_classic(base_size = 20) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
      axis.text.y = element_text(size = 18, color = "black"),
      axis.title.x = element_text(size = 20, color = "black", face = "bold"),
      axis.title.y = element_text(size = 20, color = "black", face = "bold"),
      plot.title   = element_text(size = 22, color = "black", face = "bold", hjust = 0.5),
      legend.position = "none"
    ) +
    coord_cartesian(ylim = ylimits)
}

p1 <- make_boxplot(lfq_long_filtered, "Filtered Only: LFQ Intensity Distribution", my_y_limits)
p2 <- make_boxplot(lfq_long_norm, "Quantile Normalized: LFQ Intensity Distribution", my_y_limits)

p1 + p2 + plot_layout(ncol = 2)

ggsave("boxplot_LFQ_filtered_only_official.png", p1, width = 10, height = 7, dpi = 300, bg = "white")
ggsave("boxplot_LFQ_quantile_norm_official.png", p2, width = 10, height = 7, dpi = 300, bg = "white")

################################################################################
################################################################################
################################################################################

################################################################################
################################################################################
################################################################################

# ---- [Load Libraries] ----
library(mixOmics)
library(ggplot2)
library(dplyr)
library(stringr)
library(ggrepel)

# ---- [Create Output Directory] ----
output_dir <- "sPLSDA_2D_Plots_T1_T10"
dir.create(output_dir, showWarnings = FALSE)

# ---- [Prepare Data: quantile_norm_data, columns T1_1 ... T10_8] ----
sample_cols <- c(
  paste0("T1_", 1:8), paste0("T2_", 1:8), paste0("T3_", 1:8), paste0("T4_", 1:8),
  paste0("T5_", 1:8), paste0("T6_", 1:8), paste0("T7_", 1:8), paste0("T8_", 1:8),
  paste0("T9_", 1:8), paste0("T10_", 1:8)
)
X <- as.data.frame(t(quantile_norm_data[, sample_cols]))  # samples x proteins
sample_names <- rownames(X)
Y <- factor(sub("_.*", "", sample_names), levels = paste0("T", 1:10))

# ---- [Color palette for T1–T10] ----
group_colors <- setNames(
  c("#66c2a5", "#ffd92f", "#8da0cb", "#fc8d62", "#a6cee3",
    "#1f78b4", "#ffb347", "#b3de69", "#bdbdbd", "#bc80bd"),
  paste0("T", 1:10)
)

custom_theme_splsda <- function(base_size = 14) {
  theme_minimal(base_size = base_size) +
    theme(
      panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
      plot.title = element_text(size = base_size + 4, face = "bold", hjust = 0.5),
      axis.title = element_text(size = base_size + 2, face = "bold"),
      axis.text = element_text(size = base_size),
      legend.title = element_text(size = base_size, face = "bold"),
      legend.text = element_text(size = base_size),
      legend.position = "bottom"
    )
}

# ---- [Tuning sPLS-DA] ----
set.seed(1)
optimal_ncomp <- 9
list_keepX <- c(25, 50, 100)

cat("🔍 Tuning sPLS-DA for T1-T10 groups...\n")
tune_result <- tune.splsda(X, Y, ncomp = optimal_ncomp,
                           validation = "Mfold", folds = 5,
                           dist = "centroids.dist", measure = "BER",
                           test.keepX = list_keepX, nrepeat = 10, progressBar = TRUE)
optimal_keepX <- tune_result$choice.keepX[1:optimal_ncomp]
print(optimal_keepX)

# ---- [keepX Barplot] ----
keepx_df <- data.frame(
  Component = factor(paste0("Comp", seq_along(optimal_keepX))),
  keepX = optimal_keepX
)
gg_keepx <- ggplot(keepx_df, aes(x = Component, y = keepX)) +
  geom_bar(stat = "identity", fill = "#1f78b4", alpha = 0.85, width = 0.6) +
  geom_text(aes(label = keepX), vjust = -0.5, size = 5) +
  labs(title = "Optimal keepX per Component", x = "Component", y = "Variables Selected") +
  theme_minimal(base_size = 14) +
  theme(panel.border = element_rect(color = "black", fill = NA))
print(gg_keepx)
ggsave(file.path(output_dir, "keepX_per_component_T1_T10.png"),
       plot = gg_keepx, dpi = 300, width = 6, height = 5, bg = "white")

# ---- [Fit Final Model] ----
splsda_model <- splsda(X, Y, ncomp = optimal_ncomp, keepX = optimal_keepX)

# ---- [Explained Variance] ----
expl_var <- apply(splsda_model$variates$X^2, 2, sum) / sum(splsda_model$X^2)

# ---- [Manual Ellipse Function] ----
desired_ellipse_level <- 0.95
compute_ellipse <- function(mean, cov, level = desired_ellipse_level, npoints = 100) {
  angles <- seq(0, 2 * pi, length.out = npoints)
  radius <- sqrt(qchisq(level, df = 2))
  eig <- eigen(cov)
  axes <- radius * t(eig$vectors %*% diag(sqrt(eig$values)))
  ellipse <- t(axes %*% rbind(cos(angles), sin(angles))) + matrix(rep(mean, each = npoints), ncol = 2, byrow = FALSE)
  df <- as.data.frame(ellipse)
  colnames(df) <- c("comp1", "comp2")
  return(df)
}

# ---- [Plot All Component Pairs] ----
for (i in 1:(optimal_ncomp - 1)) {
  for (j in (i + 1):optimal_ncomp) {
    
    plot_data <- data.frame(
      comp1 = splsda_model$variates$X[, i],
      comp2 = splsda_model$variates$X[, j],
      Group = Y,
      Sample = rownames(X)
    )
    
    group_centroids <- plot_data %>%
      group_by(Group) %>%
      summarise(comp1 = mean(comp1), comp2 = mean(comp2), count = n(), .groups = "drop")
    
    ellipse_data <- plot_data %>%
      group_by(Group) %>%
      do({
        group_data <- .[, c("comp1", "comp2")]   # <-- Fixed: base R select
        ell <- compute_ellipse(colMeans(group_data), cov(group_data))
        ell$Group <- unique(.$Group)
        ell
      }) %>% ungroup()
    
    x_lab <- paste0("Component ", i, " (", round(expl_var[i] * 100, 1), "%)")
    y_lab <- paste0("Component ", j, " (", round(expl_var[j] * 100, 1), "%)")
    
    p <- ggplot(plot_data, aes(x = comp1, y = comp2, color = Group, fill = Group)) +
      geom_point(size = 4, alpha = 0.9) +
      geom_polygon(data = ellipse_data, aes(group = Group), alpha = 0.18, color = NA) +
      geom_path(data = ellipse_data, aes(group = Group), linewidth = 1) +
      geom_text_repel(data = group_centroids,
                      aes(label = paste0(Group, "\n(n=", count, ")")),
                      color = "black", size = 5, fontface = "bold",
                      max.overlaps = 100, box.padding = 0.6, point.padding = 0.6) +
      scale_color_manual(values = group_colors) +
      scale_fill_manual(values = group_colors) +
      labs(
        title = paste0("sPLS-DA: (Comp ", i, " vs ", j, ")"),
        x = x_lab, y = y_lab
      ) +
      custom_theme_splsda()
    
    print(p)
    ggsave(file.path(output_dir, paste0("sPLS-DA_T1-T10_Comp", i, "_vs_Comp", j, ".png")),
           plot = p, dpi = 300, width = 10, height = 10, bg = "white")
  }
}

cat("\n✅ All sPLS-DA 2D plots saved in folder:", output_dir, "\n")

# ---- [Summary Report] ----
cat("\n📊 sPLS-DA Summary Report\n")
cat("──────────────────────────────\n")
cat("✔ Number of Components Used:", optimal_ncomp, "\n")
cat("✔ keepX per Component:", paste(optimal_keepX, collapse = ", "), "\n")
cat("✔ Confidence Level for Ellipses:", desired_ellipse_level * 100, "%\n")
n_to_report <- min(2, length(expl_var))
cat("✔ Explained Variance (First", n_to_report, "components):",
    paste0(round(expl_var[1:n_to_report] * 100, 1), collapse = "%, "), "%\n")
cat("✔ Output folder:", output_dir, "\n")

################################################################################
################################################################################
################################################################################

library(ggplot2)
library(dplyr)
library(ggrepel)

# 1. Run PCA
pca_res <- prcomp(X, center = TRUE, scale. = TRUE)
var_explained <- (pca_res$sdev)^2 / sum(pca_res$sdev^2) * 100

# 2. Prepare Data Frame
plot_data <- data.frame(
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Group = Y,
  Sample = rownames(X)
)

# 3. Group colors (same as sPLS-DA)
group_colors <- c(
  "T1" = "#66c2a5", "T2" = "#ffd92f", "T3" = "#8da0cb", "T4" = "#fc8d62", "T5" = "#a6cee3",
  "T6" = "#1f78b4", "T7" = "#ffb347", "T8" = "#b3de69", "T9" = "#bdbdbd", "T10" = "#bc80bd"
)

# 4. Compute centroids for labeling
group_centroids <- plot_data %>%
  group_by(Group) %>%
  summarise(PC1 = mean(PC1), PC2 = mean(PC2), count = n(), .groups = "drop")

# 5. Plot
p <- ggplot(plot_data, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  geom_point(size = 4, alpha = 0.9) +
  stat_ellipse(aes(group = Group, fill = Group), 
               type = "norm", level = 0.95, alpha = 0.2, geom = "polygon", color = NA) +
  stat_ellipse(aes(group = Group, color = Group), 
               type = "norm", level = 0.95, size = 1, fill = NA, geom = "path") +
  geom_text_repel(data = group_centroids,
                  aes(label = paste0(Group, "\n(n=", count, ")")),
                  color = "black", size = 5, fontface = "bold",
                  max.overlaps = 100, box.padding = 0.6, point.padding = 0.6, segment.size = 0.5) +
  scale_color_manual(values = group_colors) +
  scale_fill_manual(values = group_colors) +
  labs(
    title = "PCA: (PC1 vs PC2)",
    x = paste0("PC1 (", round(var_explained[1], 1), "%)"),
    y = paste0("PC2 (", round(var_explained[2], 1), "%)")
  ) +
  theme_minimal(base_size = 22) +
  theme(
    panel.border = element_rect(color = "black", fill = NA, linewidth = 1),
    axis.title = element_text(face = "bold"),
    plot.title = element_text(size = 24, face = "bold", hjust = 0.5),
    legend.title = element_text(size = 20, face = "bold"),
    legend.text = element_text(size = 18),
    legend.position = "bottom"
  )

print(p)
ggsave("PCA_T1-T10_PC1_vs_PC2_sPLSDAstyle.png", p, dpi = 300, width = 8, height = 8, bg = "white")



################################################################################
################################################################################
################################################################################
################################################################################
################################################################################

###############  Volcano plot with enrichment analysis  ########################

################################################################################
################################################################################
################################################################################
################################################################################
################################################################################


################################################################################
# Add Gene_name via left join (Protein IDs ↔ UniProt Entry)
# New matrix name: gene_labeled_imputed_data
################################################################################

# --- Load UniProt mapping (Entry -> Gene Names) ---
uni_file <- "C:/Users/ga53hil/Desktop/Granit_proteomics/28.08.25_result/FASTA_uniprotkb_taxonomy_id_9913_Bovine.tsv"
uniprot <- read.table(uni_file, header = TRUE, sep = "\t", quote = "",
                      check.names = FALSE, stringsAsFactors = FALSE,
                      comment.char = "", fill = TRUE)

# Safety checks
if (!"Entry" %in% colnames(uniprot)) stop("Column 'Entry' not found in UniProt TSV.")
if (!"Gene Names" %in% colnames(uniprot)) stop("Column 'Gene Names' not found in UniProt TSV.")
if (!"Protein IDs" %in% colnames(imputed_data)) stop("Column 'Protein IDs' not found in imputed_data.")

# Keep only needed columns and drop empty gene-name rows
uni_map <- uniprot[, c("Entry", "Gene Names")]
uni_map$`Gene Names`[uni_map$`Gene Names` == ""] <- NA
uni_map <- uni_map[!is.na(uni_map$`Gene Names`), ]
uni_map <- uni_map[!duplicated(uni_map$Entry), ]

# Named lookup: Entry -> Gene Names
map <- setNames(uni_map$`Gene Names`, uni_map$Entry)

# Map "P12345;Q9XXXX" -> "GENE1;GENE2" (unique, semicolon-separated)
map_gene_names <- function(prot_ids) {
  if (is.na(prot_ids) || prot_ids == "") return(NA_character_)
  ids <- trimws(strsplit(prot_ids, ";", fixed = TRUE)[[1]])
  genes <- unname(map[ids])
  genes <- genes[!is.na(genes) & genes != ""]
  if (!length(genes)) return(NA_character_)
  paste(unique(genes), collapse = ";")
}

# --- Create new matrix with gene names ---
gene_labeled_imputed_data <- imputed_data
gene_labeled_imputed_data$Gene_name <- vapply(
  gene_labeled_imputed_data$`Protein IDs`, map_gene_names, character(1)
)

# --- Save (TXT + Excel) ---
write.table(gene_labeled_imputed_data,
            "P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70_imputed_fixedseed_gene_labeled.txt",
            sep = "\t", quote = FALSE, row.names = FALSE)

library(openxlsx)
write.xlsx(gene_labeled_imputed_data,
           "P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70_imputed_fixedseed_gene_labeled.xlsx",
           overwrite = TRUE)

################################################################################
################################################################################
################################################################################

# --- [Volcano Plot Analysis: Pairwise T1–T10 Comparisons] ---

# Build an "express_matrix" from gene_labeled_imputed_data (T#_# columns),
# set rownames to cleaned Gene_name (fallback to Protein IDs), and ensure uniqueness/numeric.
expr_cols_final <- grep("^T\\d+_\\d+$", colnames(gene_labeled_imputed_data), value = TRUE)
if (!length(expr_cols_final)) stop("No expression columns (T#_#) found for volcano analysis.")

express_matrix <- as.data.frame(gene_labeled_imputed_data[expr_cols_final], check.names = FALSE)

# --- Clean & normalize gene names ---
gene_ids_raw <- ifelse(
  is.na(gene_labeled_imputed_data$Gene_name) | gene_labeled_imputed_data$Gene_name == "",
  gene_labeled_imputed_data$`Protein IDs`,
  gene_labeled_imputed_data$Gene_name
)

gene_ids <- sapply(gene_ids_raw, function(g) {
  parts <- unlist(strsplit(g, "[ ;]+"))   # split on space or semicolon
  parts <- unique(trimws(parts))          # trim spaces & remove duplicates
  parts <- parts[nzchar(parts)]           # remove empty entries
  paste(parts, collapse = ";")            # join with semicolon
}, USE.NAMES = FALSE)

rownames(express_matrix) <- make.unique(gene_ids)

# coerce all to numeric
express_matrix[] <- lapply(express_matrix, function(x) suppressWarnings(as.numeric(as.character(x))))

# Where to place volcano outputs (defaults to current working directory)
output_dir <- getwd()

# --- Required packages for volcano workflow ---
suppressPackageStartupMessages({
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
})

cat("📊 Starting volcano plot analysis for all T1–T10 comparisons...\n")

# Define group names and columns
group_names <- paste0("T", 1:10)
group_columns <- lapply(group_names, function(g) grep(paste0("^", g, "_"), colnames(express_matrix), value = TRUE))
names(group_columns) <- group_names

# Prepare a proper data.frame
data_id <- data.frame(Gene = rownames(express_matrix), express_matrix, check.names = FALSE)

# Parameters
pval_threshold <- 0.05
fc_threshold <- log2(2)

# Output directory for volcano analysis
volcano_output_dir <- file.path(output_dir, "Volcano_Results")
dir.create(volcano_output_dir, showWarnings = FALSE)

# Run pairwise comparisons: always "early" vs "late", i < j
for (i in 1:(length(group_names) - 1)) {
  for (j in (i + 1):length(group_names)) {
    group_early <- group_names[i]
    group_late  <- group_names[j]
    
    early_cols <- group_columns[[group_early]]
    late_cols  <- group_columns[[group_late]]
    
    # Skip if either group has no columns present
    if (length(early_cols) == 0 || length(late_cols) == 0) next
    
    df <- data.frame(
      Gene = data_id$Gene,
      log2FC = NA_real_,
      p_value = NA_real_
    )
    
    for (k in 1:nrow(data_id)) {
      vals_early <- suppressWarnings(as.numeric(as.character(unlist(data_id[k, early_cols]))))
      vals_late  <- suppressWarnings(as.numeric(as.character(unlist(data_id[k, late_cols]))))
      
      if (all(is.na(vals_early)) || all(is.na(vals_late))) next
      
      # Welch t-test (late vs early)
      t_test <- tryCatch(stats::t.test(vals_late, vals_early, var.equal = FALSE), error = function(e) NULL)
      if (!is.null(t_test)) {
        df$log2FC[k] <- mean(vals_late, na.rm = TRUE) - mean(vals_early, na.rm = TRUE)
        df$p_value[k] <- t_test$p.value
      }
    }
    
    df <- stats::na.omit(df)
    if (!nrow(df)) next
    
    df$negLog10P <- -log10(df$p_value)
    
    # Label directions
    label_up   <- paste("Higher in", group_late)   # late = right side
    label_down <- paste("Higher in", group_early)  # early = left side
    
    df$group <- "Non-significant"
    df$group[df$log2FC > fc_threshold & df$p_value < pval_threshold] <- label_up
    df$group[df$log2FC < -fc_threshold & df$p_value < pval_threshold] <- label_down
    
    # Top genes for plotting (not needed for file export)
    top_up <- df %>% dplyr::filter(group == label_up) %>% dplyr::arrange(p_value) %>% head(10)
    top_down <- df %>% dplyr::filter(group == label_down) %>% dplyr::arrange(p_value) %>% head(10)
    
    # Output folder for this comparison
    comp <- paste0(group_early, "_vs_", group_late)
    comp_dir <- file.path(volcano_output_dir, comp)
    dir.create(comp_dir, showWarnings = FALSE)
    
    # --- Output only the three desired files ---
    # 1. All results
    write.table(df, file = file.path(comp_dir, paste0("All_proteins_", comp, ".txt")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    # 2. Up in late, down in early (T10_up_T1_down_T1_vs_T10.txt)
    write.table(df[df$group == label_up, ], 
                file = file.path(comp_dir, paste0(group_late, "_up_", group_early, "_down_", comp, ".txt")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    # 3. Up in early, down in late (T1_up_T10_down_T1_vs_T10.txt)
    write.table(df[df$group == label_down, ], 
                file = file.path(comp_dir, paste0(group_early, "_up_", group_late, "_down_", comp, ".txt")),
                sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Volcano plot (optional: can comment out if not needed)
    color_map <- setNames(c("red", "blue", "black"),
                          c(label_up, label_down, "Non-significant"))
    volcano_plot <- ggplot(df, aes(x = log2FC, y = negLog10P)) +
      geom_point(aes(color = group), alpha = 0.7, size = 3.5) +
      scale_color_manual(values = color_map) +
      geom_hline(yintercept = -log10(pval_threshold), linetype = "dashed", color = "black", size = 1) +
      geom_vline(xintercept = c(-fc_threshold, fc_threshold), linetype = "dashed", color = "black", size = 1) +
      ggrepel::geom_text_repel(
        data = top_up, aes(label = Gene), color = "red", size = 10,
        box.padding = 0.25, point.padding = 0.25, max.overlaps = Inf
      ) +
      ggrepel::geom_text_repel(
        data = top_down, aes(label = Gene), color = "blue", size = 10,
        box.padding = 0.25, point.padding = 0.25, max.overlaps = Inf
      ) +
      labs(
        title = paste("Volcano Plot:", comp),
        x = "Log2 Fold Change",
        y = "-Log10 p-value",
        color = "Expression"
      ) +
      theme_minimal(base_size = 34) +
      theme(
        plot.title = element_text(size = 45, face = "bold", hjust = 0.5),
        axis.title = element_text(size = 45, face = "bold"),
        axis.text  = element_text(size = 32),
        legend.title = element_text(size = 36),
        legend.text  = element_text(size = 36),
        legend.position = "bottom"
      )
    
    ggsave(
      file.path(comp_dir, paste0("Volcano_", comp, ".png")),
      plot = volcano_plot, dpi = 300, width = 18, height = 18, bg = "white"
    )
    
    cat("✅ Volcano plot saved for", comp, "\n")
    
  }
}

cat("🎉 Volcano plot analysis complete. Results in:", normalizePath(volcano_output_dir), "\n")


################################################################################
################################################################################
################################################################################

# ============================
# Volcano Result Annotation Script (Updated for New File Names)
# ============================

# --- Load required libraries ---
suppressPackageStartupMessages({
  library(org.Bt.eg.db)
  library(AnnotationDbi)
  library(dplyr)
  library(tidyr)
  library(ReactomePA)
  library(biomaRt)
  library(KEGGREST)
  library(GO.db)
  library(tools)
})

# --- Check that volcano_output_dir exists ---
if (!exists("volcano_output_dir")) {
  stop("❌ 'volcano_output_dir' not found in environment. Run the volcano script first.")
}

# --- Helper function: split multi-gene names and keep mapping back ---
split_genes <- function(gene_vec) {
  gene_list <- strsplit(gene_vec, ";")
  data.frame(
    Original = rep(gene_vec, times = sapply(gene_list, length)),
    Gene_Single = unlist(gene_list),
    stringsAsFactors = FALSE
  )
}

# --- Get list of all comparison folders ---
comparison_folders <- list.dirs(volcano_output_dir, recursive = FALSE, full.names = TRUE)

# --- Loop through each volcano comparison folder ---
for (comp_folder in comparison_folders) {
  comp_name <- basename(comp_folder)
  
  # Load the full volcano table
  all_file <- list.files(comp_folder, pattern = "^All_proteins_.*\\.txt$", full.names = TRUE)
  if (!length(all_file)) {
    cat("❌ No All_proteins file found for", comp_name, "\n")
    next
  }
  Whole_data <- read.table(all_file[1], header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)
  
  # Find ONLY the two up/down files with the new pattern
  updown_files <- list.files(
    comp_folder,
    pattern = "^(T\\d+_up_T\\d+_down_.*\\.txt|T\\d+_up_T\\d+_down_.*\\.txt)$",
    full.names = TRUE
  )
  if (!length(updown_files)) {
    cat("❌ No up/down files found for", comp_name, "\n")
    next
  }
  
  for (gene_file in updown_files) {
    set_label <- tools::file_path_sans_ext(basename(gene_file))
    cat("🔍 Annotating set:", set_label, "in", comp_name, "...\n")
    
    gene_df <- read.table(gene_file, header = TRUE, sep = "\t", quote = "", fill = TRUE,
                          comment.char = "", stringsAsFactors = FALSE)
    if (!"Gene" %in% colnames(gene_df)) {
      cat("❌ 'Gene' column not found in", gene_file, "\n")
      next
    }
    
    # --- Split multi-gene entries ---
    gene_map <- split_genes(gene_df$Gene)
    
    # --- Gene to Entrez mapping ---
    gene_entrez <- AnnotationDbi::select(org.Bt.eg.db,
                                         keys = unique(gene_map$Gene_Single),
                                         columns = "ENTREZID",
                                         keytype = "SYMBOL")
    combined_entrez <- unique(na.omit(gene_entrez$ENTREZID))
    
    # --- GO Annotation ---
    go_results <- tryCatch({
      go_anno <- AnnotationDbi::select(org.Bt.eg.db, keys = combined_entrez,
                                       columns = c("GO", "ONTOLOGY"), keytype = "ENTREZID")
      go_terms <- AnnotationDbi::select(GO.db, keys = unique(go_anno$GO),
                                        columns = "TERM", keytype = "GOID")
      go_merged <- merge(go_anno, go_terms, by.x = "GO", by.y = "GOID", all.x = TRUE)
      go_merged <- left_join(go_merged, gene_entrez, by = "ENTREZID", relationship = "many-to-many")
      go_merged$GO_Pathway <- paste(go_merged$GO, go_merged$ONTOLOGY, go_merged$TERM, sep = " | ")
      go_merged %>% dplyr::select(SYMBOL, GO_Pathway)
    }, error = function(e) NULL)
    
    go_agg <- if (!is.null(go_results)) {
      go_results %>%
        group_by(SYMBOL) %>%
        summarise(GO_Pathway = paste(unique(GO_Pathway), collapse = "; "), .groups = "drop")
    } else data.frame(SYMBOL = character(), GO_Pathway = character())
    
    # --- Reactome via Human Orthologs ---
    cow_mart <- useMart("ensembl", dataset = "btaurus_gene_ensembl",
                        host = "https://dec2021.archive.ensembl.org")
    human_mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl",
                          host = "https://dec2021.archive.ensembl.org")
    
    gene_chunks <- split(unique(gene_map$Gene_Single), ceiling(seq_along(unique(gene_map$Gene_Single)) / 50))
    orthologs <- do.call(rbind, lapply(gene_chunks, function(chunk) {
      tryCatch(getLDS(attributes = "external_gene_name", filters = "external_gene_name",
                      values = chunk, mart = cow_mart,
                      attributesL = c("entrezgene_id", "external_gene_name"), martL = human_mart),
               error = function(e) NULL)
    }))
    colnames(orthologs) <- c("Cow_Gene", "Human_ENTREZID", "Human_Gene")
    human_entrez <- unique(na.omit(orthologs$Human_ENTREZID))
    
    reactome_results <- tryCatch({
      enr <- enrichPathway(
        gene = human_entrez,
        organism = "human",
        readable = TRUE,
        pvalueCutoff = 1,
        qvalueCutoff = 1
      )
      df <- as.data.frame(enr)
      df %>%
        dplyr::select(Reactome_Pathway = Description, geneID) %>%
        tidyr::separate_rows(geneID, sep = "/") %>%
        left_join(orthologs, by = c("geneID" = "Human_Gene"), relationship = "many-to-many") %>%
        dplyr::rename(SYMBOL = Cow_Gene) %>%
        dplyr::select(SYMBOL, Reactome_Pathway)
    }, error = function(e) NULL)
    
    reactome_agg <- if (!is.null(reactome_results)) {
      reactome_results %>%
        group_by(SYMBOL) %>%
        summarise(Reactome_Pathway = paste(unique(Reactome_Pathway), collapse = "; "), .groups = "drop")
    } else data.frame(SYMBOL = character(), Reactome_Pathway = character())
    
    # --- KEGG Annotation ---
    formatted_entrez <- paste0("bta:", combined_entrez)
    kegg_data <- tryCatch(keggLink("pathway", "bta"), error = function(e) NULL)
    kegg_data_filtered <- kegg_data[names(kegg_data) %in% formatted_entrez]
    
    gene_info <- AnnotationDbi::select(org.Bt.eg.db, keys = combined_entrez,
                                       columns = c("ENTREZID", "SYMBOL"), keytype = "ENTREZID")
    
    kegg_descriptions <- tryCatch({
      desc <- stack(keggList("pathway", "bta"))
      data.frame(KEGG_Pathway = sub("path:", "", desc$ind),
                 KEGG_Description = desc$values)
    }, error = function(e) NULL)
    
    kegg_results <- data.frame(
      ENTREZID = sub("bta:", "", names(kegg_data_filtered)),
      SYMBOL = gene_info$SYMBOL[match(sub("bta:", "", names(kegg_data_filtered)), gene_info$ENTREZID)],
      KEGG_Pathway = sub("path:", "", kegg_data_filtered)
    )
    
    kegg_results <- left_join(kegg_results, kegg_descriptions, by = "KEGG_Pathway", relationship = "many-to-many")
    kegg_results$KEGG_Annotation <- paste(kegg_results$KEGG_Description, kegg_results$KEGG_Pathway, sep = " | ")
    
    kegg_agg <- kegg_results %>%
      group_by(SYMBOL) %>%
      summarise(kegg_results = paste(unique(KEGG_Annotation), collapse = "; "), .groups = "drop")
    
    # --- Merge annotations back to original multi-gene entries ---
    combined_annots <- gene_map %>%
      left_join(go_agg, by = c("Gene_Single" = "SYMBOL")) %>%
      left_join(reactome_agg, by = c("Gene_Single" = "SYMBOL")) %>%
      left_join(kegg_agg, by = c("Gene_Single" = "SYMBOL")) %>%
      group_by(Original) %>%
      summarise(
        GO_Pathway = paste(na.omit(unique(GO_Pathway)), collapse = "; "),
        Reactome_Pathway = paste(na.omit(unique(Reactome_Pathway)), collapse = "; "),
        kegg_results = paste(na.omit(unique(kegg_results)), collapse = "; "),
        .groups = "drop"
      )
    
    annotated_set <- gene_df %>%
      left_join(combined_annots, by = c("Gene" = "Original"))
    
    # Save annotated set
    annotated_file <- file.path(comp_folder, paste0(set_label, "_annotated.txt"))
    write.table(annotated_set, annotated_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    # Merge with full volcano dataset
    annotated_set$Annotation_Label <- "+"
    merged_matrix <- left_join(Whole_data, annotated_set, by = "Gene", relationship = "many-to-many")
    merged_matrix[merged_matrix == "" | merged_matrix == " "] <- NA
    
    merged_file <- file.path(comp_folder, paste0("Merged_matrix_", set_label, ".txt"))
    write.table(merged_matrix, merged_file, sep = "\t", row.names = FALSE, quote = FALSE)
    
    cat("✅ Annotation complete for", set_label, "in", comp_name, "\n")
  }
}

cat("\n🎉 Annotation complete for all volcano comparison sets.\n")


################################################################################
################################################################################
################################################################################

# ==========================================================
# Fisher Exact Test Enrichment for Volcano Results (Cleaned for new naming)
# ==========================================================

library(dplyr)
library(tidyr)
library(stringr)
library(openxlsx)

# --- Working directory ---
setwd("C:/Users/ga53hil/Desktop/Granit_proteomics/28.08.25_result")

# --- Background file ---
bg_file <- "P337_02B_proteinGroups_filtered_pepGT1_log2LFQ_valid70_imputed_fixedseed_gene_labeled.txt"
background_data <- read.table(bg_file, header = TRUE, sep = "\t", quote = "", check.names = FALSE, stringsAsFactors = FALSE)
background_data$Gene_name[background_data$Gene_name == "" | is.na(background_data$Gene_name)] <- "Unknown"
background_data$Gene_name <- make.unique(background_data$Gene_name)
background_genes <- unique(background_data$Gene_name)

# --- Volcano results root ---
volcano_root <- file.path(getwd(), "Volcano_Results")

# --- Fisher test function ---
fisher_test_annotation <- function(df, source_col, sig_genes, background_genes, output_file) {
  df <- df %>%
    filter(!is.na(.data[[source_col]]) & .data[[source_col]] != "") %>%
    mutate(Pathway = strsplit(as.character(.data[[source_col]]), ";")) %>%
    unnest(Pathway) %>%
    mutate(Pathway = str_trim(Pathway)) %>%
    filter(!is.na(Pathway) & Pathway != "" & Pathway != "NA")
  
  # --- Remove "NA | NA | NA" from GO ---
  df <- df[!grepl("^NA\\s*\\|\\s*NA\\s*\\|\\s*NA$", df$Pathway, ignore.case = TRUE), ]
  
  # --- Clean pathway display ---
  df <- df %>%
    mutate(Pathway_Name = case_when(
      source_col == "kegg_results" ~ sub(" - Bos.*", "", Pathway),
      grepl("\\|", Pathway) ~ sub(".*\\|\\s*", "", Pathway),
      TRUE ~ Pathway
    ))
  
  fisher_results <- data.frame()
  
  for (i in unique(df$Pathway_Name)) {
    term_genes <- unique(df$Gene[df$Pathway_Name == i & df$Gene %in% background_genes])
    sig_in_pathway <- intersect(term_genes, sig_genes)
    
    a <- length(sig_in_pathway)
    b <- length(sig_genes) - a
    c <- length(term_genes) - a
    d <- length(background_genes) - a - b - c
    
    mat <- matrix(c(a, b, c, d), nrow = 2, byrow = TRUE)
    if (any(mat < 0) || any(!is.finite(mat))) next
    
    test <- tryCatch(fisher.test(mat), error = function(e) NULL)
    if (!is.null(test)) {
      gene_ratio <- ifelse((a + b) > 0, a / (a + b), 0)
      bg_ratio   <- ifelse((a + b + c + d) > 0, (a + c) / (a + b + c + d), 0)
      enrichment_factor <- ifelse(bg_ratio > 0, gene_ratio / bg_ratio, NA)
      
      fisher_results <- rbind(fisher_results, data.frame(
        Pathway = i,
        Sig_protein_volcano = a,
        Sig_NotIn_Pathway   = b,
        Nonsig_In_Pathway   = c,
        Nonsig_NotIn_Pathway = d,
        Genes = paste(sig_in_pathway, collapse = ";"),
        p_value = test$p.value,
        enrichment_factor = enrichment_factor
      ))
    }
  }
  
  if (nrow(fisher_results) > 0) {
    fisher_results$p_adj <- p.adjust(fisher_results$p_value, method = "BH")
    fisher_results <- fisher_results[order(fisher_results$p_adj), ]
    write.xlsx(fisher_results, output_file, rowNames = FALSE)
  }
}

# --- Loop through Volcano comparison folders ---
comparisons <- list.dirs(volcano_root, recursive = FALSE, full.names = TRUE)

for (comp in comparisons) {
  comp_name <- basename(comp)
  
  # Find only new annotated up/down files (T#_up_T#_down_T#_vs_T#_annotated.txt)
  annotated_files <- list.files(comp, 
                                pattern = "^T\\d+_up_T\\d+_down_T\\d+_vs_T\\d+_annotated\\.txt$", 
                                full.names = TRUE)
  
  if (!length(annotated_files)) {
    cat("❌ No annotated volcano sets for", comp_name, "\n")
    next
  }
  
  for (annot_file in annotated_files) {
    set_label <- tools::file_path_sans_ext(basename(annot_file))
    cat("🔍 Running Fisher for:", set_label, "in", comp_name, "\n")
    
    ann_df <- read.table(annot_file, header = TRUE, sep = "\t", quote = "", stringsAsFactors = FALSE)
    if (!"Gene" %in% colnames(ann_df)) {
      cat("❌ No Gene column in", annot_file, "\n")
      next
    }
    
    sig_genes <- unique(ann_df$Gene)
    
    # Output folder
    out_dir <- file.path(comp, "Fisher")
    dir.create(out_dir, showWarnings = FALSE)
    
    # Run Fisher for GO, KEGG, Reactome
    if ("GO_Pathway" %in% colnames(ann_df)) {
      fisher_test_annotation(ann_df, "GO_Pathway", sig_genes, background_genes,
                             file.path(out_dir, paste0("Fisher_GO_", set_label, ".xlsx")))
    }
    if ("kegg_results" %in% colnames(ann_df)) {
      fisher_test_annotation(ann_df, "kegg_results", sig_genes, background_genes,
                             file.path(out_dir, paste0("Fisher_KEGG_", set_label, ".xlsx")))
    }
    if ("Reactome_Pathway" %in% colnames(ann_df)) {
      fisher_test_annotation(ann_df, "Reactome_Pathway", sig_genes, background_genes,
                             file.path(out_dir, paste0("Fisher_Reactome_", set_label, ".xlsx")))
    }
  }
}

cat("\n🎉 Fisher exact test complete for all Volcano results.\n")


################################################################################
################################################################################
################################################################################

# Volcano Enrichment Heatmaps (GO + KEGG + Reactome) [MERGED SOURCES]
# - Merges pathway names with same text, shows all sources in []
# - Filters: adj p < 0.01, protein count > 5
# - Readable A4 output with paging and keyword filter
# ==========================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(ggplot2)
  library(openxlsx)
  library(viridisLite)   # for colors
  library(scales)        # for rescale
})

# ---- User Parameters ----
volcano_root <- "Volcano_Results"
global_keywords <- c(
  # Metabolism
  "amino acid metabolism", "translation initiation", "rRNA processing",
  "mitochondrial translation", "lipid metabolism", "carbohydrate metabolism",
  "central carbon metabolism", "nucleotide metabolism", "glucose metabolism",
  "pyruvate metabolism", "Citrate cycle (TCA cycle)", "glycolysis", "gluconeogenesis", 
  "beta oxidation", "fatty acid degradation", "steroid metabolism", 
  "cholesterol metabolism", "glutathione metabolism", "methionine metabolism", 
  "polyamine metabolism", "biosynthesis of amino acids",
  # Signaling
  "WNT signaling", "MAPK signaling", "interleukin signaling", "TNF signaling",
  "JAK-STAT signaling", "TGF-beta signaling", "interferon signaling", 
  "NOTCH signaling", "EGFR signaling", "HIF-1 signaling", "TP53 signaling",
  "TCR signaling", "B cell receptor", "chemokine signaling",
  # Stress & Apoptosis
  "UPR", "autophagy", "apoptosis", "programmed cell death", 
  "cellular senescence", "hypoxia response", "DNA repair", "DNA replication",
  "mitotic cell cycle", "cell cycle checkpoints",
  # Immune / Matrix / Transport
  "neutrophil degranulation", "extracellular matrix", "platelet activation",
  "lysosome", "endosome", "vesicle-mediated transport", "deubiquitination",
  "proteolysis", "biological oxidations", "organelle biogenesis",
  "ECM-receptor interaction",
  # Development / Morphogenesis
  "gastrulation", "somitogenesis", "chromatin organization", "histone modification",
  # Reproduction-specific
  "oocyte", "spermatogenesis", "meiosis", "gamete generation", 
  "gonad development", "reproductive system", "reproductive development",
  "sexual reproduction", "fertilization", "androgen", "estrogen", "luteinizing hormone", 
  "follicle-stimulating hormone", "hormone signaling"
)

P_ADJ_CUTOFF <- 0.01              # Adjusted p-value cutoff for pathway filtering
MIN_PROTEINS <- 6                 # Minimum number of proteins required for a pathway to be shown
MAX_PATHWAYS_PER_PAGE <- 15       # Max number of pathways per page/plot
MAX_GENES_PER_PAGE <- 60         # Max number of gene/protein columns per page/plot (wider plot!)
PATHWAY_WRAP_WIDTH <- 40          # Width to wrap pathway names (characters before breaking line)
A4_W <- 16                        # Width of output plot (inches) -- increase to fit more gene labels
A4_H <- 10                        # Height of output plot (inches)
BASE_FONTSIZE <- 14               # Base font size for theme_minimal()
AXIS_FONTSIZE_X <- 8              # Font size for gene/protein names (x-axis)
AXIS_FONTSIZE_Y <- 20             # Font size for pathway names (y-axis)
TITLE_FONTSIZE <- 20              # Font size for plot title
LEGEND_TITLE_SIZE <- 16           # Font size for legend title
LEGEND_TEXT_SIZE <- 16            # Font size for legend text
LEGEND_KEY_HEIGHT_PT <- 12        # Height of legend color keys (pt)
LEGEND_KEY_WIDTH_PT <- 10         # Width of legend color keys (pt)


# ---- Helper functions ----
first_col <- function(df, cand) {
  c <- cand[cand %in% names(df)]
  if (length(c)) c[1] else NA_character_
}
count_genes <- function(x) {
  if (is.null(x) || is.na(x)) return(0L)
  toks <- unlist(stringr::str_split(as.character(x), "[;,/\\s]+"))
  toks <- stringr::str_trim(toks)
  toks <- toks[toks != ""]
  length(unique(toks))
}

# ---- Find all Fisher/ folders ----
all_dirs    <- list.dirs(volcano_root, recursive = TRUE, full.names = TRUE)
fisher_dirs <- all_dirs[grepl("[/\\\\]Fisher$", all_dirs)]  # ends with /Fisher or \Fisher

if (!length(fisher_dirs)) {
  stop("❌ No 'Fisher' folders found under: ", normalizePath(volcano_root))
}

# ---- Main loop: each Fisher directory (comparison) ----
for (fdir in fisher_dirs) {
  cat("\n📂 Processing:", fdir, "\n")
  files <- list.files(
    fdir,
    pattern = "^Fisher_(GO|KEGG|Reactome)_T\\d+_up_T\\d+_down_T\\d+_vs_T\\d+_annotated\\.xlsx$",
    full.names = TRUE
  )
  if (!length(files)) {
    cat("⚠️  No annotated Fisher enrichment files found in", fdir, "\n")
    next
  }
  
  # Parse file info for direction/comparison grouping
  info <- stringr::str_match(
    basename(files),
    "^Fisher_(GO|KEGG|Reactome)_(T\\d+_up_T\\d+_down_T\\d+_vs_T\\d+)_annotated\\.xlsx$"
  )
  colnames(info) <- c("full","Source","Direction")
  if (any(is.na(info))) next
  directions <- unique(info[, "Direction"])
  
  for (direction in directions) {
    cat("   🔎 Direction:", direction, "\n")
    f_direction <- files[info[, "Direction"] == direction]
    src         <- info[info[, "Direction"] == direction, "Source"]
    
    # Parse comp for titling
    direction_parts <- stringr::str_match(direction, "^(T\\d+)_up_(T\\d+)_down_(T\\d+_vs_T\\d+)$")
    side <- direction_parts[2]
    comp <- direction_parts[4]
    
    pieces <- list()
    for (k in seq_along(f_direction)) {
      fpath       <- f_direction[k]
      source_type <- src[k]
      df0 <- tryCatch(openxlsx::read.xlsx(fpath), error = function(e) NULL)
      if (is.null(df0) || !"Pathway" %in% names(df0)) {
        cat("    ⚠️  Skipping (no Pathway):", basename(fpath), "\n")
        next
      }
      pcol <- first_col(df0, c("p_adj", "p.adj", "padj", "adj_p_value", "adj.P.Val", "p_adjust"))
      gcol <- first_col(df0, c("Genes", "Gene", "Gene_List", "GeneID", "Gene_Names", "Proteins", "Protein_List"))
      if (is.na(pcol) || is.na(gcol)) {
        cat("    ⚠️  Skipping (missing p-adj or gene column):", basename(fpath), "\n")
        next
      }
      # --- Filtering and explosion ---
      df <- df0 %>%
        dplyr::mutate(
          .p         = suppressWarnings(as.numeric(.data[[pcol]])),
          .genes_raw = as.character(.data[[gcol]]),
          .gene_n    = vapply(.genes_raw, count_genes, integer(1))
        ) %>%
        dplyr::filter(is.finite(.p), .p < P_ADJ_CUTOFF, .gene_n >= MIN_PROTEINS) %>%
        dplyr::filter(
          stringr::str_detect(
            tolower(Pathway),
            paste(tolower(global_keywords), collapse = "|")
          )
        ) %>%
        dplyr::mutate(
          Pathway = stringr::str_trim(Pathway),
          Source  = source_type,
          log10_p = -log10(.p)
        ) %>%
        tidyr::separate_rows(dplyr::all_of(gcol), sep = "[;,/\\s]+") %>%
        dplyr::mutate(Gene = stringr::str_trim(.data[[gcol]])) %>%
        dplyr::filter(Gene != "")
      pieces[[length(pieces) + 1]] <- df
    }
    
    side_df <- dplyr::bind_rows(pieces)
    if (!nrow(side_df)) {
      cat("   ⚠️  No GO/KEGG/Reactome rows for", direction, "after filtering\n")
      next
    }
    
    # ---- Combine same-named pathways and collect sources ----
    # Merge on pathway name (no [Source]), collect all sources, keep all genes
    source_levels <- c("KEGG", "GO", "Reactome")
    side_df <- side_df %>%
      dplyr::mutate(
        Pathway_Only = stringr::str_remove(Pathway, "\\s*\\[.*\\]$"),
        Source = factor(Source, levels = source_levels)
      ) %>%
      dplyr::group_by(Pathway_Only, Gene) %>%
      dplyr::summarise(
        log10_p = max(log10_p, na.rm = TRUE),
        SourceList = paste(intersect(source_levels, unique(as.character(Source))), collapse = ", "),
        .groups = "drop"
      ) %>%
      dplyr::mutate(
        Pathway = paste0(stringr::str_wrap(Pathway_Only, width = PATHWAY_WRAP_WIDTH),
                         " [", SourceList, "]")
      )
    
    # ---- Limit genes for readability: most frequent across pathways ----
    gene_order_by_freq <- side_df %>% dplyr::count(Gene, sort = TRUE) %>% dplyr::pull(Gene)
    keep_genes <- head(gene_order_by_freq, MAX_GENES_PER_PAGE)
    side_df <- side_df %>% dplyr::filter(Gene %in% keep_genes)
    
    mat   <- side_df %>%
      dplyr::select(Pathway, Gene, log10_p) %>%
      tidyr::complete(Pathway, Gene = keep_genes)
    # ---- Pathway ordering: by number of proteins (highest first) ----
    ord <- mat %>%
      dplyr::filter(!is.na(log10_p)) %>%
      dplyr::group_by(Pathway) %>%
      dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
      dplyr::arrange(dplyr::desc(n)) %>%
      dplyr::pull(Pathway)
    mat <- mat %>%
      dplyr::mutate(
        Pathway = factor(Pathway, levels = ord),
        Gene    = factor(Gene, levels = keep_genes)
      )
    
    out_dir <- file.path(fdir, "Combined_Enrichment_Plots", direction)
    dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
    
    # ---- Page split for A4 output ----
    pths_all <- levels(mat$Pathway)
    chunks <- split(pths_all, ceiling(seq_along(pths_all) / MAX_PATHWAYS_PER_PAGE))
    for (i in seq_along(chunks)) {
      sel <- chunks[[i]]
      d   <- mat %>%
        dplyr::filter(Pathway %in% sel) %>%
        dplyr::mutate(Pathway = factor(Pathway, levels = rev(sel)))
      title_txt <- paste0(direction, " (", comp, ") — Page ", i)
      caption_txt <- paste0("Filters: adj p < ", P_ADJ_CUTOFF, "; proteins > ", MIN_PROTEINS - 1, "; genes shown capped at ", MAX_GENES_PER_PAGE, " by frequency.")
      p <- ggplot(d, aes(x = Gene, y = Pathway, fill = log10_p)) +
        geom_tile(color = "white", linewidth = 0.25, na.rm = FALSE) +
        scale_fill_gradientn(
          colours = c("green","red"),
          na.value = "black",
          name = "-log10(p.adj)"
        ) +
        scale_x_discrete(position = "top", guide = guide_axis(check.overlap = TRUE)) +
        labs(title = title_txt, x = "Significant Gene", y = "Enriched Pathway", caption = caption_txt) +
        theme_minimal(base_size = BASE_FONTSIZE) +
        theme(
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x  = element_text(angle = 60, hjust = 0, vjust = 0.9, size = AXIS_FONTSIZE_X, margin = margin(b = 2)),
          axis.text.y  = element_text(size = AXIS_FONTSIZE_Y),
          axis.title.x = element_text(size = AXIS_FONTSIZE_Y, margin = margin(t = 6)),
          axis.title.y = element_text(size = AXIS_FONTSIZE_Y, margin = margin(r = 8)),
          plot.title   = element_text(face = "bold", hjust = 0.5, size = TITLE_FONTSIZE, margin = margin(b = 10)),
          legend.position = "right",
          legend.title = element_text(size = LEGEND_TITLE_SIZE),
          legend.text  = element_text(size = LEGEND_TEXT_SIZE),
          legend.key.height = unit(LEGEND_KEY_HEIGHT_PT, "pt"),
          legend.key.width  = unit(LEGEND_KEY_WIDTH_PT,  "pt"),
          plot.caption = element_text(size = 10, hjust = 1),
          plot.margin = margin(t = 16, r = 16, b = 16, l = 16)
        ) +
        coord_cartesian(clip = "off", expand = FALSE)
      out_base <- file.path(out_dir, paste0("Enrichment_", direction, "_", comp, "_Page", i))
      ggsave(paste0(out_base, ".pdf"), plot = p, width = A4_W, height = A4_H, units = "in", device = cairo_pdf, bg = "white")
      ggsave(paste0(out_base, ".png"), plot = p, width = A4_W, height = A4_H, units = "in", dpi = 600, bg = "white", limitsize = FALSE)
      cat("   ✅ Saved A4 PDF+PNG:", out_base, "(.pdf/.png)\n")
    }
  }
}
cat("\n🎉 All enrichment plots generated with merged pathway sources and A4-optimized readability.\n")

################################################################################
################################################################################
################################################################################
################################################################################

# ----------------------------------------------------------
# Volcano Sankey: Top 5 Pathways Per Timepoint (Hybrid: Each Pathway at Peak)
# Each pathway assigned to timepoint where it is strongest (peak)
# No repeats, Up to 5 per TP, ordered by timepoint (T1–T10)
# Timepoints appear INSIDE the node rectangles (networkD3 + CSS)
# ----------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(stringr)
  library(openxlsx)
  library(RColorBrewer)
  library(networkD3)
  library(htmlwidgets)
  library(webshot2)
})

cat("🔎 Starting Volcano Sankey: Hybrid (Best of Unique + Distinctive) Per T1-T10\n")

volcano_root <- file.path(getwd(), "Volcano_Results")
if (!dir.exists(volcano_root)) stop("❌ 'Volcano_Results' folder not found at: ", volcano_root)

timepoints <- paste0("T", 1:10)
tp_colors  <- RColorBrewer::brewer.pal(10, "Paired")
names(tp_colors) <- timepoints

FDR_CUTOFF <- 0.01
MIN_PROTEINS <- 20

global_keywords <- c(
  # Metabolism
  "amino acid metabolism", "translation initiation", "rRNA processing",
  "mitochondrial translation", "lipid metabolism", "carbohydrate metabolism",
  "central carbon metabolism", "nucleotide metabolism", "glucose metabolism",
  "pyruvate metabolism", "Citrate cycle (TCA cycle)", "glycolysis", "gluconeogenesis", 
  "beta oxidation", "fatty acid degradation", "steroid metabolism", 
  "cholesterol metabolism", "glutathione metabolism", "methionine metabolism", 
  "polyamine metabolism", "biosynthesis of amino acids",
  
  # Signaling
  "WNT signaling", "MAPK signaling", "interleukin signaling", "TNF signaling",
  "JAK-STAT signaling", "TGF-beta signaling", "interferon signaling", 
  "NOTCH signaling", "EGFR signaling", "HIF-1 signaling", "TP53 signaling",
  "TCR signaling", "B cell receptor", "chemokine signaling",
  
  # Stress & Apoptosis
  "UPR", "autophagy", "apoptosis", "programmed cell death", 
  "cellular senescence", "hypoxia response", "DNA repair", "DNA replication",
  "mitotic cell cycle", "cell cycle checkpoints",
  
  # Immune / Matrix / Transport
  "neutrophil degranulation", "extracellular matrix", "platelet activation",
  "lysosome", "endosome", "vesicle-mediated transport", "deubiquitination",
  "proteolysis", "biological oxidations", "organelle biogenesis",
  "ECM-receptor interaction",
  
  # Development / Morphogenesis
  "gastrulation", "somitogenesis", "chromatin organization", "histone modification",
  
  # Reproduction-specific
  "oocyte", "spermatogenesis", "meiosis", "gamete generation", 
  "gonad development", "reproductive system", "reproductive development",
  "sexual reproduction", "fertilization", "androgen", "estrogen", "luteinizing hormone", 
  "follicle-stimulating hormone", "hormone signaling"
)

file_pattern <- "^Fisher_(GO|KEGG|Reactome)_(T\\d+)_up_(T\\d+)_down_(T\\d+_vs_T\\d+)_annotated\\.xlsx$"

# --- Helper: Read Fisher and tidy ---
read_and_tidy <- function(filepath) {
  bn <- basename(filepath)
  m  <- stringr::str_match(bn, file_pattern)
  if (any(is.na(m))) {
    cat("⚠️  Skipping unmatched filename: ", bn, "\n")
    return(NULL)
  }
  src <- m[2]
  first_tp <- m[3]
  second_tp <- m[4]
  comp <- m[5]
  side <- if (startsWith(comp, first_tp)) "First" else "Second"
  up_in <- first_tp
  df <- tryCatch(openxlsx::read.xlsx(filepath), error = function(e) NULL)
  if (is.null(df)) {
    cat("❌ Could not read file: ", bn, "\n")
    return(NULL)
  }
  gene_col <- if ("Genes" %in% names(df)) "Genes" else if ("Gene" %in% names(df)) "Gene" else NA_character_
  if (is.na(gene_col) || !"Pathway" %in% names(df)) {
    cat("⚠️  File lacks Pathway or Gene(s): ", bn, "\n")
    return(NULL)
  }
  if ("p_adj" %in% names(df)) df$p_adj <- suppressWarnings(as.numeric(df$p_adj)) else df$p_adj <- NA_real_
  df2 <- df %>%
    dplyr::mutate(Source = src, Up_in = up_in, Side = side, Comparison = comp) %>%
    dplyr::filter(!is.na(.data[[gene_col]]), .data[[gene_col]] != "") %>%
    dplyr::mutate(Pathway = stringr::str_trim(as.character(Pathway))) %>%
    dplyr::filter(!is.na(Pathway), Pathway != "", !(Pathway %in% c("NA", "NA | NA | NA")))
  if (nrow(df2) == 0) {
    cat("⚠️  File has no usable rows after filter: ", bn, "\n")
    return(NULL)
  }
  if (!is.null(global_keywords) && length(global_keywords) > 0) {
    pat <- paste0("(", paste0(global_keywords, collapse = "|"), ")")
    df2  <- df2 %>% dplyr::filter(stringr::str_detect(stringr::str_to_lower(Pathway), stringr::str_to_lower(pat)))
    if (nrow(df2) == 0) {
      cat("⚠️  No keyword-matched pathways in: ", bn, "\n")
      return(NULL)
    }
  }
  df_long <- df2 %>%
    tidyr::separate_rows(.data[[gene_col]], sep = "[;,/\\s]+") %>%
    dplyr::mutate(Gene = stringr::str_trim(.data[[gene_col]])) %>%
    dplyr::filter(Gene != "") %>%
    dplyr::distinct(Pathway, Gene, Up_in, Side, Source, Comparison, p_adj)
  if (nrow(df_long) == 0) {
    cat("⚠️  No usable long-form rows: ", bn, "\n")
    return(NULL)
  }
  df_long <- df_long %>% dplyr::mutate(Timepoint = factor(Up_in, levels = timepoints))
  return(df_long)
}

# --- Load Fisher files ---
cat("🔎 Scanning all Fisher enrichment files...\n")
all_fisher_dirs <- list.dirs(volcano_root, recursive = TRUE, full.names = TRUE)
fdirs <- all_fisher_dirs[grepl(paste0(.Platform$file.sep, "Fisher$"), all_fisher_dirs)]
if (length(fdirs) == 0) stop("❌ No Fisher folders found under: ", volcano_root)
for (d in fdirs) cat("   Scanning:", normalizePath(d), "\n")

all_rows <- list()
for (d in fdirs) {
  files <- list.files(d, pattern = file_pattern, full.names = TRUE)
  if (length(files) == 0) next
  for (f in files) {
    r <- read_and_tidy(f)
    if (!is.null(r)) all_rows[[length(all_rows)+1]] <- r
  }
}
if (length(all_rows) == 0) stop("❌ No enrichment rows collected. Check file names, contents or keyword filter.")
dat <- dplyr::bind_rows(all_rows)

# --- Deduplicate and filter by pathway stats ---
dat <- dat %>%
  dplyr::group_by(Pathway, Gene, Timepoint, Side, Source, Comparison) %>%
  dplyr::summarise(p_adj = ifelse(all(is.na(p_adj)), NA_real_, min(p_adj, na.rm = TRUE)), .groups = "drop")

# --- Filtering (protein count > 5, FDR < 0.01) ---
path_stats_tp <- dat %>%
  dplyr::group_by(Pathway, Timepoint) %>%
  dplyr::summarise(
    min_p = ifelse(all(is.na(p_adj)), NA_real_, min(p_adj, na.rm = TRUE)),
    gene_count = dplyr::n_distinct(Gene),
    .groups = "drop"
  ) %>%
  dplyr::filter(gene_count >= MIN_PROTEINS, !is.na(min_p), min_p < FDR_CUTOFF) %>%
  dplyr::mutate(tp_num = as.numeric(stringr::str_replace(Timepoint, "T", "")))

# --- Assign each pathway to its peak TP (highest gene_count, best p-value, earliest TP) ---
peak_tp <- path_stats_tp %>%
  dplyr::group_by(Pathway) %>%
  dplyr::filter(gene_count == max(gene_count, na.rm = TRUE)) %>%   # max gene_count per pathway
  dplyr::arrange(min_p, tp_num) %>%     # break ties by p-value, then earliest TP
  dplyr::slice(1) %>%                  # pick only one timepoint per pathway
  dplyr::ungroup()

# --- For each TP, select up to top 5 pathways "peaking" there ---
final_top_paths <- peak_tp %>%
  dplyr::group_by(Timepoint) %>%
  dplyr::arrange(desc(gene_count), min_p) %>%
  dplyr::slice_head(n = 5) %>%
  dplyr::ungroup() %>%
  dplyr::pull(Pathway) %>%
  unique()

cat(sprintf("✅ Selected %d pathways (top 5 per timepoint, each at peak only).\n", length(final_top_paths)))

dat_sel <- dat %>% dplyr::filter(Pathway %in% final_top_paths)

# --- Pathway node order: by peak timepoint (T1-T10, then by strength) ---
pathways_ordered <- peak_tp %>%
  dplyr::filter(Pathway %in% final_top_paths) %>%
  dplyr::mutate(tp_num = as.numeric(stringr::str_replace(Timepoint, "T", ""))) %>%
  dplyr::arrange(tp_num, desc(gene_count), min_p) %>%
  dplyr::pull(Pathway)

# --- Sankey weighting as before ---
path_counts <- dat_sel %>%
  dplyr::group_by(Pathway) %>%
  dplyr::summarise(n_link = dplyr::n_distinct(Timepoint), .groups = "drop")
dat_w <- dat_sel %>%
  dplyr::left_join(path_counts, by = "Pathway") %>%
  dplyr::group_by(Pathway, Timepoint, Side, Comparison) %>%
  dplyr::summarise(weight = sum(1 / ifelse(is.na(n_link) | n_link == 0, 1, n_link)), .groups = "drop") %>%
  dplyr::mutate(Timepoint = factor(as.character(Timepoint), levels = timepoints))

# --- Sankey plot function: pathway order by peak TP ---
plot_hybrid_sankey <- function(df, label, outdir, pathways_ordered) {
  if (is.null(df) || nrow(df) == 0) { cat("❌ No data for ", label, "\n"); return(NULL) }
  pathways <- pathways_ordered
  tps_order <- timepoints
  nodes <- rbind(
    data.frame(name = pathways, type = "Pathway", stringsAsFactors = FALSE),
    data.frame(name = tps_order,  type = "Timepoint", stringsAsFactors = FALSE)
  )
  df <- df %>%
    dplyr::mutate(
      Pathway = factor(Pathway, levels = pathways),
      Timepoint = factor(Timepoint, levels = tps_order)
    )
  links <- df %>%
    dplyr::mutate(
      source = match(as.character(Pathway), nodes$name) - 1,
      target = match(as.character(Timepoint), nodes$name) - 1,
      value = weight,
      group = as.character(Timepoint)
    ) %>%
    dplyr::select(source, target, value, group)
  missing_tps <- setdiff(tps_order, unique(as.character(df$Timepoint)))
  if (length(missing_tps) > 0) {
    dummy_links <- data.frame(
      source = 0,
      target = match(missing_tps, nodes$name) - 1,
      value = 0.0001,
      group = missing_tps,
      stringsAsFactors = FALSE
    )
    links <- dplyr::bind_rows(links, dummy_links)
  }
  js_colors <- sprintf('d3.scaleOrdinal().domain(%s).range(%s)',
                       jsonlite::toJSON(timepoints),
                       jsonlite::toJSON(unname(tp_colors)))
  sankey <- sankeyNetwork(
    Links = links,
    Nodes = nodes,
    Source = "source",
    Target = "target",
    Value = "value",
    NodeID = "name",
    fontSize = 36,
    nodeWidth = 36,
    sinksRight = TRUE,
    LinkGroup = "group",
    colourScale = js_colors,
    width = 1200, height = 800,
    nodePadding = 20,
    iterations = 0
  )
  outfile <- file.path(outdir, paste0('Sankey_', label, '.html'))
  htmlwidgets::saveWidget(sankey, file = outfile, selfcontained = TRUE, background = "#fff")
  # CSS to remove node rectangles
  lines <- readLines(outfile, warn = FALSE)
  css_patch <- '
  <style>
    .node rect {
      fill: none !important;
      stroke: none !important;
    }
    .node text {
      font-weight: bold;
      fill: #222;
      text-anchor: middle !important;
      x: 0 !important;
      transform: none !important;
      font-size: 32px !important;
    }
    body, html { background: #fff !important; }
  </style>
  '
  idx <- grep("</head>", lines, fixed = TRUE)
  if (length(idx) > 0) {
    lines <- append(lines, css_patch, after = idx[1] - 1)
    writeLines(lines, outfile)
  }
  cat("✅ Saved Sankey HTML: ", outfile, "\n")
  outfile_png <- file.path(outdir, paste0('Sankey_', label, '.png'))
  webshot2::webshot(
    url = outfile,
    file = outfile_png,
    vwidth = 1200,
    vheight = 750,
    delay = 3,
    zoom = 4  # <--- Try 2, 3, or 4 for higher DPI/print resolution!
  )
  cat("✅ Also saved as PNG: ", outfile_png, "\n")
  invisible(sankey)
}

# --- Output directory and plot ---
out_root <- file.path(volcano_root, "Sankey_Top5HybridPerTimepoint_T1toT10_Manuscript")
if (!dir.exists(out_root)) dir.create(out_root, recursive = TRUE, showWarnings = FALSE)

cat("🔎 Plotting Top 5 Pathways per Timepoint (Hybrid)...\n")
plot_hybrid_sankey(dat_w, "Top5HybridPerTimepoint_T1toT10", out_root, pathways_ordered)

cat("\n🎉 Sankey plot (Top 5 *hybrid* Pathways/TP, T1–T10, FDR<0.01 & protein>5) as HTML and PNG in:\n", normalizePath(out_root), "\n")
cat("🎉 Done.\n")




################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
################################################################################
