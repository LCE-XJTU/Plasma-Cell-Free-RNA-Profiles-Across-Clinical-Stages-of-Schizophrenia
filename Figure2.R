# ============================================================================
# Figure 2A 
# ============================================================================

work_dir <- "Figure2"
setwd(work_dir)
OUT <- work_dir 

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(scales)
  library(grid)
  library(patchwork)
})

deg_dir <- file.path("DEG")
matrix_file <- file.path("fig2A_upset_matrix_new.csv")

m <- read.csv(matrix_file, check.names = FALSE, stringsAsFactors = FALSE)
rownames(m) <- m$gene

comp_labels <- c("SCZ vs HC","FES vs HC", "MCS vs HC", "LTS vs HC")
comp_pairs <- c("HC  SCZ","HC  FES", "HC  MCS", "HC  LTS")

mat <- as.matrix(m[, comp_labels, drop = FALSE])
stopifnot(colnames(mat) == comp_labels)
mode(mat) <- "numeric"
mat <- mat[rowSums(mat) > 0, ]
cat("Genes in matrix:", nrow(mat), "\n")

dir_df_list <- list()
for (i in seq_along(comp_pairs)) {
  src <- file.path(deg_dir, comp_pairs[i], "diffgene.csv")
  if (!file.exists(src)) {
    cat("Warning: Cannot find", src, "\n")
    next
  }
  
  df <- read.csv(src, header = TRUE, row.names = 1, check.names = FALSE)
  df$gene <- rownames(df)
  
  df_sig <- df %>% 
    filter(!is.na(change) & change %in% c("up", "down")) %>%
    dplyr::select(gene, change) %>%
    mutate(Group = comp_labels[i])
  
  dir_df_list[[i]] <- df_sig
}
dir_df <- bind_rows(dir_df_list)
dir_df <- dir_df %>% distinct(gene, Group, .keep_all = TRUE)

gene_comb <- apply(mat, 1, function(x) paste(x, collapse = ""))
comb_stats <- data.frame()

for (comb_str in unique(gene_comb)) {
  genes_in_comb <- names(gene_comb[gene_comb == comb_str])
  s <- as.numeric(unlist(strsplit(comb_str, "")))
  sets_in <- comp_labels[which(s == 1)]
  
  sub_dir <- dir_df %>% filter(gene %in% genes_in_comb & Group %in% sets_in)
  
  if (nrow(sub_dir) > 0) {
    gene_dirs <- sub_dir %>% 
      group_by(gene) %>% 
      summarise(
        n_dirs = n_distinct(change),
        dir = first(change),
        .groups = "drop"
      )
    n_up <- sum(gene_dirs$n_dirs == 1 & gene_dirs$dir == "up")
    n_down <- sum(gene_dirs$n_dirs == 1 & gene_dirs$dir == "down")
    n_discordant <- sum(gene_dirs$n_dirs > 1)
  } else {
    n_up <- 0; n_down <- 0; n_discordant <- 0
  }
  
  comb_stats <- rbind(comb_stats, data.frame(
    Comb = comb_str,
    Sets = paste(sets_in, collapse = " + "),
    N_genes = length(genes_in_comb),
    N_up = n_up,
    N_down = n_down,
    N_discordant = n_discordant,
    stringsAsFactors = FALSE
  ))
}

comb_stats <- comb_stats %>% arrange(desc(N_genes))
full_comb <- "11111"
if (full_comb %in% comb_stats$Comb) {
  other <- comb_stats[comb_stats$Comb != full_comb, ]
  full_row <- comb_stats[comb_stats$Comb == full_comb, ]
  comb_stats <- rbind(other, full_row)
}
comb_stats$Comb <- factor(comb_stats$Comb, levels = comb_stats$Comb)

max_pos <- max(comb_stats$N_up + comb_stats$N_discordant)
max_neg <- max(comb_stats$N_down)

ylim_max <- max_pos * 1.25
ylim_min <- -max_neg * 1.25

comb_stats <- comb_stats %>%
  mutate(
    Y_up_label = case_when(
      N_up == 0 ~ NA_real_,
      N_discordant > 0 ~ N_up / 2,
      TRUE ~ N_up + max_pos * 0.03
    ),
    Y_disc_label = N_up + (N_discordant / 2),
    Y_down_label = -N_down - max_neg * 0.03
  )
p_bar <- ggplot(comb_stats, aes(x = Comb)) +

  geom_col(aes(y = N_up + N_discordant), fill = "#B0B0B0", width = 0.7, alpha = 0.85) +
  geom_col(aes(y = N_up), fill = "#c82423", width = 0.7, alpha = 0.85) +
  geom_col(aes(y = -N_down), fill = "#2878b5", width = 0.7, alpha = 0.85) +
  geom_hline(yintercept = 0, linewidth = 0.3) +
  geom_text(data = comb_stats %>% filter(N_up > 0 & N_discordant == 0),
            aes(y = Y_up_label, label = N_up), size = 2.8, color = "#c82423", vjust = 0) +
  geom_text(data = comb_stats %>% filter(N_up > 0 & N_discordant > 0),
            aes(y = Y_up_label, label = N_up), size = 2.8, color = "#c82423", fontface = "bold") +
  geom_text(data = comb_stats %>% filter(N_discordant > 0),
            aes(y = Y_disc_label, label = N_discordant), 
            size = 2.3, color = "black", fontface = "bold",
            nudge_x = 0.12) +
  geom_text(data = comb_stats %>% filter(N_down > 0),
            aes(y = Y_down_label, label = N_down), size = 2.8, color = "#2878b5", vjust = 1) +
  
  scale_y_continuous(limits = c(ylim_min, ylim_max), expand = c(0, 0)) +
  labs(x = NULL, y = "DEG count") +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_text(size = 7),
    axis.title.y = element_text(size = 8),
    plot.margin = margin(5, 5, 0, 5),
    panel.grid = element_blank()
  )
dot_levels <- c("SCZ vs HC", "FES vs HC", "MCS vs HC", "LTS vs HC")
n_sets <- length(comp_labels)
dot_data <- data.frame()

for (i in 1:nrow(comb_stats)) {
  comb_str <- as.character(comb_stats$Comb[i])
  s <- as.numeric(unlist(strsplit(comb_str, "")))
  idx_map <- match(dot_levels, comp_labels) 
  s_reordered <- s[idx_map]
  
  for (j in 1:n_sets) {
    dot_data <- rbind(dot_data, data.frame(
      Comb = comb_stats$Comb[i],
      Set = dot_levels[j],
      Present = s_reordered[j],
      stringsAsFactors = FALSE
    ))
  }
}
dot_data$Comb <- factor(dot_data$Comb, levels = levels(comb_stats$Comb))
dot_data$Set <- factor(dot_data$Set, levels = rev(dot_levels))

line_data <- dot_data %>%
  filter(Present == 1) %>%
  group_by(Comb) %>%
  summarise(
    ymin = as.numeric(min(as.numeric(Set))),
    ymax = as.numeric(max(as.numeric(Set))),
    .groups = "drop"
  )

p_dot <- ggplot(dot_data, aes(x = Comb, y = Set)) +
  geom_linerange(data = line_data,
                 aes(x = Comb, ymin = ymin, ymax = ymax),
                 linewidth = 0.8, color = "#2B2B2B", inherit.aes = FALSE) +
  geom_point(aes(size = factor(Present), color = factor(Present)), alpha = 0.85) +
  scale_size_manual(values = c("0" = 0.5, "1" = 3), guide = "none") +
  scale_color_manual(values = c("0" = "#E0E0E0", "1" = "#2B2B2B"), guide = "none") +
  labs(x = NULL, y = NULL) +
  theme_classic(base_size = 8) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 9),
    axis.text.y = element_text(size = 9, face = "bold"),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_line(linewidth = 0.3),
    plot.margin = margin(0, 5, 5, 5),
    panel.grid = element_blank()
  )

p_combined <- p_bar / p_dot +
  plot_layout(heights = c(1, 1)) +
  plot_annotation(title = "DEG overlap and direction across clinical stages",
                  theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0.5)))

ggsave(file.path(OUT, "fig2A_upset_direction2.png"), p_combined, width = 5, height = 5, dpi = 300)
ggsave(file.path(OUT, "fig2A_upset_direction2.pdf"), p_combined, width = 5, height = 5)

# ============================================================================
# Figure 2B
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(readr)
})

get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) > 0) {
    return(dirname(normalizePath(sub("^--file=", "", file_arg[1]), winslash = "/")))
  }

  if (requireNamespace("rstudioapi", quietly = TRUE) &&
      rstudioapi::isAvailable()) {
    path <- rstudioapi::getSourceEditorContext()$path
    if (nzchar(path)) {
      return(dirname(normalizePath(path, winslash = "/")))
    }
  }
  
  normalizePath(getwd(), winslash = "/")
}
script_dir <- get_script_dir()
project_root <- normalizePath(file.path(script_dir, "..", ".."), winslash = "/")

input_file <- file.path("06_Stage_vs_Matched_HC_adjusted_model.csv")
output_pdf <- file.path(script_dir, "Figure2B_stage_module_grouped_forest_R_matchedHC.pdf")
output_png <- file.path(script_dir, "Figure2B_stage_module_grouped_forest_R_matchedHC.png")

if (!file.exists(input_file)) {
  stop(
    "Input file was not found:\n", input_file,
    "\nIf the script was moved, set project_root manually near the top of the script."
  )
}

stage_levels <- c("FES vs HC", "MCS vs HC", "LTS vs HC")
stage_colors <- c(
  "FES vs HC" = "#E31A1C",
  "MCS vs HC" = "#FF7F00",
  "LTS vs HC" = "#33A02C"
)

module_info <- tibble::tribble(
  ~Module,       ~module_order, ~module_label, ~function_label,
  "MEblue",              4,    "Blue",        "Platelet / adhesion / cytoskeleton",
  "MEbrown",             3,    "Brown",       "Neurodevelopment / cell fate",
  "MEturquoise",         2,    "Turquoise",   "Membrane potential / synaptic",
  "MEyellow",            1,    "Yellow",      "Chromatin organization"
)

stage_info <- tibble::tribble(
  ~Comparison, ~stage,       ~y_offset,
  "StageFES",  "FES vs HC",  0.10,
  "StageMCS",  "MCS vs HC",  0.00,
  "StageLTS",  "LTS vs HC", -0.10
)

plot_data <- read_csv(input_file, show_col_types = FALSE) %>%
  filter(Comparison %in% stage_info$Comparison) %>%
  inner_join(stage_info, by = "Comparison") %>%
  inner_join(module_info, by = "Module") %>%
  mutate(
    stage = factor(stage, levels = stage_levels),
    significant = factor(
      if_else(adj.P.Val < 0.05, "FDR < .05", "FDR >= .05"),
      levels = c("FDR < .05", "FDR >= .05")
    ),
    y = module_order + y_offset
  )

row_bands <- tibble::tibble(
  ymin = c(3.55, 1.55),
  ymax = c(4.45, 2.45)
)

module_names <- module_info %>%
  transmute(x = -0.108, y = module_order + 0.09, label = module_label)

module_functions <- module_info %>%
  transmute(x = -0.108, y = module_order - 0.10, label = function_label)

cap_height <- 0.035

p <- ggplot(plot_data) +
  geom_rect(
    data = row_bands,
    aes(xmin = -Inf, xmax = Inf, ymin = ymin, ymax = ymax),
    inherit.aes = FALSE,
    fill = "#FAFAFA"
  ) +
  geom_vline(
    xintercept = seq(-0.06, 0.04, by = 0.02),
    color = "#E6E6E6",
    linewidth = 0.35
  ) +
  geom_vline(
    xintercept = 0,
    color = "#777777",
    linewidth = 0.55,
    linetype = "dashed"
  ) +
  geom_segment(
    aes(x = CI.L, xend = CI.R, y = y, yend = y, color = stage),
    linewidth = 0.9
  ) +
  geom_segment(
    aes(
      x = CI.L, xend = CI.L,
      y = y - cap_height, yend = y + cap_height,
      color = stage
    ),
    linewidth = 0.7
  ) +
  geom_segment(
    aes(
      x = CI.R, xend = CI.R,
      y = y - cap_height, yend = y + cap_height,
      color = stage
    ),
    linewidth = 0.7
  ) +
  geom_point(
    aes(x = logFC, y = y, color = stage, shape = significant),
    size = 3.5,
    stroke = 1.15
  ) +
  geom_text(
    data = module_names,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    fontface = "bold",
    size = 3.7,
    color = "#202020"
  ) +
  geom_text(
    data = module_functions,
    aes(x = x, y = y, label = label),
    inherit.aes = FALSE,
    hjust = 0,
    size = 2.75,
    color = "#666666"
  ) +
  scale_color_manual(
    name = NULL,
    values = stage_colors,
    breaks = stage_levels
  ) +
  scale_shape_manual(
    name = NULL,
    values = c("FDR < .05" = 16, "FDR >= .05" = 1)
  ) +
  scale_x_continuous(
    breaks = seq(-0.06, 0.04, by = 0.02),
    labels = function(x) sprintf("%.2f", x),
    limits = c(-0.11, 0.055),
    expand = expansion(mult = c(0, 0))
  ) +
  scale_y_continuous(limits = c(0.55, 4.55), breaks = NULL) +
  labs(
    x = "Adjusted eigengene difference (SCZ subgroup - HC)",
    y = NULL,
    caption = "Points represent adjusted differences; horizontal lines indicate 95% CIs."
  ) +
  guides(
    color = guide_legend(order = 1, nrow = 1, override.aes = list(shape = 16, size = 3.2)),
    shape = guide_legend(order = 2, nrow = 1)
  ) +
  theme_classic(base_family = "Arial", base_size = 10) +
  theme(
    plot.title = element_text(face = "bold", size = 16, margin = margin(b = 3)),
    plot.subtitle = element_text(color = "#555555", size = 9.5, margin = margin(b = 6)),
    plot.caption = element_text(color = "#666666", size = 8, hjust = 0.5, margin = margin(t = 6)),
    axis.title.x = element_text(size = 10.5, margin = margin(t = 8)),
    axis.text.x = element_text(size = 9, color = "#333333"),
    axis.line.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.position = "top",
    legend.justification = "right",
    legend.box = "horizontal",
    legend.margin = margin(t = -4, b = 2),
    legend.spacing.x = grid::unit(5, "pt"),
    plot.margin = margin(16, 18, 10, 18)
  )

pdf_device <- if (capabilities("cairo")) grDevices::cairo_pdf else grDevices::pdf

ggsave(
  filename = output_pdf,
  plot = p,
  width = 12,
  height = 10,
  units = "in",
  device = pdf_device
)

ggsave(
  filename = output_png,
  plot = p,
  width = 12,
  height = 10,
  units = "in",
  dpi = 300,
  bg = "white"
)

# ============================================================================
# Figure 2C
# ============================================================================

nes_wide <-  read_excel("fig2C_heatmap.xlsx", sheet = "NES")
padj_wide <-  read_excel("fig2C_heatmap.xlsx",  sheet = "Padj")

compare_cols <- c("Overall SCZ vs HC","FES vs HC","MCS vs HC","LTS vs HC")
comp_labels <- c("Overall_SCZ","FES","MCS","LTS")

nes_long <- nes_wide %>%
  dplyr::select(pathway=Description, Theme, all_of(compare_cols)) %>%
  pivot_longer(
    cols = all_of(compare_cols),
    names_to = "comparison",
    values_to = "NES"
  ) %>%
  mutate(NES = ifelse(is.na(NES), "None", as.character(NES)))

padj_long <- padj_wide %>%
  dplyr::select(pathway=Description, Theme, all_of(compare_cols)) %>%
  pivot_longer(
    cols = all_of(compare_cols),
    names_to = "comparison",
    values_to = "padj"
  ) %>%
  mutate(padj = ifelse(is.na(padj), "None", as.character(padj)))

df_long <- inner_join(
  nes_long,
  padj_long,
  by = c("pathway","Theme","comparison")
) %>%
  mutate(
    comp_label = case_when(
      comparison == "Overall SCZ vs HC" ~ "Overall_SCZ",
      comparison == "FES vs HC" ~ "FES",
      comparison == "MCS vs HC" ~ "MCS",
      comparison == "LTS vs HC" ~ "LTS"
    ),
    theme = Theme
  ) %>%
  mutate(
    sig = case_when(
      padj == "None" ~ FALSE,
      as.numeric(padj) < 0.05 ~ TRUE,
      TRUE ~ FALSE
    )
  ) %>%
  dplyr::select(pathway, comparison, comp_label, theme, NES, padj, sig)

write.csv(df_long,file = "fig1d_long.csv",row.names = F,quote = F)


suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(scales)
  library(patchwork)
  library(cowplot)
  library(svglite)
  library(ragg)
})

d <- read.csv("fig1d_long.csv",
              fileEncoding = "UTF-8", stringsAsFactors = FALSE,
              na.strings = c("None", "NA", ""))
d <- subset(d,comparison%in%c( "FES","MCS","LTS"))

d <- d %>%
  group_by(pathway) %>%
  filter( !all(is.na(NES)) ) %>%
  ungroup()


d$NES <- as.numeric(d$NES)
d$padj <- as.numeric(d$padj)
d$sig <- !is.na(d$padj) & d$padj < 0.05

d$sig_fill <- ifelse(d$sig & !is.na(d$NES), d$NES, NA_real_)
d$label    <- ifelse(!is.na(d$NES), sprintf("%.2f", d$NES), "")
d$theme <- sub("Gene-expression machinery", "Gene-expression\nmachinery", d$theme)

d$theme    <- factor(d$theme,
                     levels = c("Immune", "Synaptic", "Mitochondria",
                                "Platelet", "Gene-expression\nmachinery"))


path_order <- unique(d$pathway)
d$pathway  <- factor(d$pathway, levels = rev(path_order))

comp_disp <- c(FES = "FES vs HC",
               MCS = "MCS vs HC",
               LTS  = "LTS vs HC")
d$comp_label <- comp_disp[d$comparison]
d$comparison <- factor(d$comp_label, levels = comp_disp)

nes_lim <- ceiling(max(abs(d$NES), na.rm = TRUE))

comp_levels <- c(
  "FES vs HC",
  "MCS vs HC",
  "LTS vs HC")

comp_color <- c(
  "FES vs HC"         = "#E31A1C",
  "MCS vs HC"         = "#FF7F00",
  "LTS vs HC"          = "#33A02C")

type_color <- c(
  Immune                      = "#FBB4AE", 
  Synaptic                    = "#F1B56D",
  Mitochondria                = "#82C4B6", 
  Platelet                    = "#88D8D9",
  "Gene-expression\nmachinery" = "#71B8EC"
)
nes_scale <- scale_fill_gradient2(
  name = "NES",
  low = "#2878b5", mid = "#F7F7F7", high = "#c82423",
  midpoint = 0, limits = c(-nes_lim, nes_lim),
  na.value = "#E8E8E8"
)

strip_data <- d %>% distinct(pathway, theme)

p_left_strip <- ggplot(strip_data, aes(x = 1, y = pathway, fill = theme)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = type_color, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = levels(d$pathway), drop = FALSE) +
  theme_void() +
  theme(
    plot.margin = margin(2, 0, 2, 2)
  )

p_main <- ggplot(d, aes(x = comparison, y = pathway)) +
  geom_tile(aes(fill = sig_fill), colour = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 4, colour = "black") +
  geom_vline(xintercept = 1.5, linewidth = 0.4, colour = "black") +
  geom_vline(xintercept = 5.5, linewidth = 0.9, colour = "black") +
  nes_scale +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_y_discrete(limits = levels(d$pathway), position = "right", drop = FALSE) +
  guides(fill = guide_colourbar(order = 1, barheight = unit(35, "mm"),
                                barwidth  = unit(3,  "mm"))) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    axis.text.x.top = element_blank(),
    axis.ticks.x.top = element_blank(),
    axis.text.y.right = element_text(size = 14, colour = "black", hjust = 0),
    axis.title  = element_blank(),
    panel.grid  = element_blank(),
    legend.position = "right",
    legend.title    = element_text(size = 12, face = "bold"),
    legend.text     = element_text(size = 12),
    plot.margin     = margin(2, 2, 2, 2)
  )

top_tile_data <- data.frame(
  comparison = factor(comp_levels, levels = comp_levels),
  dummy = "Top"
)

p_comp_tile <- ggplot(top_tile_data, aes(x = comparison, y = dummy, fill = comparison)) +
  geom_tile(color = "white", linewidth = 0.4) +
  geom_vline(xintercept = 1.5, linewidth = 0.4, colour = "black") +
  geom_vline(xintercept = 5.5, linewidth = 0.9, colour = "black") +
  scale_fill_manual(values = comp_color, guide = "none") +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_y_discrete(drop = FALSE) +
  theme_void() +
  theme(
    plot.margin = margin(0, 2, 0, 2),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )

legend_nes <- get_legend(p_main)

p_comp_legend_dummy <- ggplot(top_tile_data, aes(x = comparison, y = dummy, fill = comparison)) +
  geom_tile() +
  scale_fill_manual(values = comp_color, name = "Comparisons") +
  theme_minimal() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10)
  )
legend_comp_box <- get_legend(p_comp_legend_dummy)

p_type_legend <- ggplot(strip_data, aes(x = theme, y = 1, fill = theme)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = type_color, name = "Pathway type") +
  theme_void() +
  theme(
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(0, 0, 0, 0)
  )
legend_type <- get_legend(p_type_legend)

p_right_leg <- plot_grid(
  legend_comp_box, 
  legend_nes, 
  legend_type, 
  ncol = 1, 
  align = "v", 
  axis = "l",
  rel_heights = c(1.2, 1.2, 1)
)

p_main_noleg <- p_main + theme(legend.position = "none")
left_heatmap_block <- plot_grid(p_left_strip, p_main_noleg, ncol = 2, rel_widths = c(0.4, 8.5), align = "h")
left_col <- plot_grid(p_comp_tile, left_heatmap_block, ncol = 1, rel_heights = c(0.35, 9.5), align = "v")

p_final <- plot_grid(left_col, p_right_leg, ncol = 2, rel_widths = c(8.5, 1.8))

out <- "fig2C_gsea_heatmap"
dir.create(dirname(out), showWarnings = FALSE)

w_in <- 210 / 25.4
h_in <- 120 / 25.4

svglite::svglite(paste0(out, ".svg"), width = w_in, height = h_in)
print(p_final); dev.off()

grDevices::cairo_pdf(paste0(out, ".pdf"), width = w_in, height = h_in)
print(p_final); dev.off()

ragg::agg_tiff(paste0(out, ".tiff"), width = w_in, height = h_in,
               units = "in", res = 600)
print(p_final); dev.off()

ragg::agg_png(paste0(out, ".png"), width = w_in, height = h_in,
              units = "in", res = 300)
print(p_final); dev.off()


# ============================================================================
# Figure 2D
# ============================================================================

base_dir <- "Results_OptimizedNestedCV"
subgroups <- c("HC_FES","HC_MCS","HC_LTS")

all_roc_objs <- list()
auc_summary <- data.frame()

for (sg in subgroups) {
  
  rds_path <- file.path(base_dir, sg, "RDS_Objects", "LOOCV_Fold_Objects.rds")
  oof_path <- file.path(base_dir, sg, "RDS_Objects", "NestedCV_OOF_Predictions.rds")
  
  if (!file.exists(rds_path) || !file.exists(oof_path)) {
    cat("NO FILE\n")
    next
  }
  
  fold_objs <- readRDS(rds_path)
  oof_pred <- readRDS(oof_path)

  roc_obj <- roc(oof_pred$Observed, oof_pred$Probability, 
                 levels = c("Control", "Case"), direction = "<", quiet = TRUE)
  
  ci_val <- ci.auc(roc_obj)
  all_roc_objs[[sg]] <- roc_obj

  cat(sprintf("    AUC: %.3f (95%% CI: %.3f - %.3f)\n", ci_val[2], ci_val[1], ci_val[3]))
  
  auc_summary <- rbind(auc_summary, data.frame(
    Subgroup = sg,
    AUC = ci_val[2],
    CI_Lower = ci_val[1],
    CI_Upper = ci_val[3]
  ))
}

if (length(all_roc_objs) > 0) {

  roc_plot_data <- data.frame()
  for (sg in names(all_roc_objs)) {
    r <- all_roc_objs[[sg]]
    sg_name <- gsub("HC__", "", sg) 
    temp_df <- data.frame(
      Specificity = r$specificities,
      Sensitivity = r$sensitivities,
      Subgroup = sg_name
    )
    ci_info <- subset(auc_summary, Subgroup == sg)
    label <- sprintf("%s (AUC = %.3f, 95%% CI: %.3f-%.3f)", 
                     sg_name, ci_info$AUC, ci_info$CI_Lower, ci_info$CI_Upper)
    temp_df$LegendLabel <- label
    roc_plot_data <- rbind(roc_plot_data, temp_df)
  }
  
  roc_plot_data$OneMinusSpec <- 1 - roc_plot_data$Specificity
  
  
  line_colors <- c("HC_FES" = "#E41A1C", "HC_MCS" = "#FF7F00", "HC_LTS" = "#4DAF4A")
  legend_map <- unique(roc_plot_data[, c("Subgroup", "LegendLabel")])
  legend_labels <- setNames(legend_map$LegendLabel, legend_map$Subgroup)
  
  p_oof_roc <- ggplot(roc_plot_data, aes(x = OneMinusSpec, y = Sensitivity, color = Subgroup)) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "grey50", size = 0.8) +
    geom_path(size = 1.2) +
  scale_color_manual(values = line_colors, labels = legend_labels) +
    scale_x_continuous(expand = c(0.01, 0.01), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    scale_y_continuous(expand = c(0.01, 0.01), limits = c(0, 1), breaks = seq(0, 1, 0.25)) +
    labs(
      title = "Nested Out-of-Fold (OOF) ROC Curves",
      x = "1 - Specificity",
      y = "Sensitivity",
      color = NULL
    ) +
    theme_bw(base_size = 14) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 16),
      legend.position = c(0.65, 0.25),
      legend.background = element_rect(fill = "white", color = "black"),
      legend.title = element_blank(),
      panel.grid.minor = element_blank()
    )
  
  ggsave(file.path(base_dir, "F2D_Nested_OOF_ROC.png"), p_oof_roc, width = 6, height = 6, dpi = 300)
  ggsave(file.path(base_dir, "F2D_Nested_OOF_ROC.pdf"), p_oof_roc, width = 6, height = 6, dpi = 300)
}

