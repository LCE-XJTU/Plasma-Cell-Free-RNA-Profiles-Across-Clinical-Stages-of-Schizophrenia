# ============================================================================
# Figure 4A
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2); library(dplyr); library(tidyr); library(readxl)
  library(patchwork); library(ggrepel); library(scales)
})

fdr_star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
comps <- c("SCZ vs HC","FES vs HC","MCS vs HC","LTS vs HC","LTS-MDI vs HC","LTS-nMDI vs HC", "LTS-MDI vs LTS-nMDI")
comp_groups <- list(
  "SCZ vs HC"=c("SCZ","HC"),
  "FES vs HC"=c("FES","HC"),
  "MCS vs HC"=c("MCS","HC"),
  "LTS vs HC"=c("LTS","HC"),
  "LTS-MDI vs HC"=c("LTS-MDI","HC"),
  "LTS-nMDI vs HC"=c("LTS-nMDI","HC"),
  "LTS-MDI vs LTS-nMDI"=c("LTS-LTS-MDI","LTS-nMDI")
)
comp_col <- c(
  "SCZ vs HC"="#c82423",
  "FES vs HC"="#E31A1C",
  "MCS vs HC"="#FF7F00",
  "LTS vs HC"="#33A02C",
  "LTS-MDI vs HC"="#1F78B4",
  "LTS-nMDI vs HC"="#984EA3",
  "LTS-MDI vs LTS-nMDI"="#FFED6F"
)
brain_family <- c(
  "L2/3-6 intratelencephalic projecting glutamatergic neuron" = "Excitatory neurons",
  "L5 extratelencephalic projecting glutamatergic cortical neuron" = "Excitatory neurons",
  "L6b glutamatergic cortical neuron" = "Excitatory neurons",
  "corticothalamic-projecting glutamatergic cortical neuron" = "Excitatory neurons",
  "near-projecting glutamatergic cortical neuron" = "Excitatory neurons",
  "chandelier pvalb GABAergic cortical interneuron" = "Inhibitory neurons",
  "lamp5 GABAergic cortical interneuron" = "Inhibitory neurons",
  "pvalb GABAergic cortical interneuron" = "Inhibitory neurons",
  "sncg GABAergic cortical interneuron" = "Inhibitory neurons",
  "sst GABAergic cortical interneuron" = "Inhibitory neurons",
  "VIP GABAergic cortical interneuron" = "Inhibitory neurons",
  "caudal ganglionic eminence derived cortical interneuron" = "Inhibitory neurons",
  "astrocyte of the cerebral cortex" = "Glial cells",
  "microglial cell" = "Glial cells",
  "oligodendrocyte" = "Glial cells",
  "oligodendrocyte precursor cell" = "Glial cells",
  "cerebral cortex endothelial cell" = "vascular cells",
  "vascular leptomeningeal cell" = "vascular cells"
)

fam_col <- c("vascular cells"="#88D8D9","Glial cells"="#82C4B6","Inhibitory neurons"="#F1B56D","Excitatory neurons"="#FBB4AE")
fam_order <- c("vascular cells","Glial cells","Inhibitory neurons","Excitatory neurons")

a_list <- list()
for (c in comps) {
  grp <- strsplit(c, " vs ")[[1]]
  folder_name <- paste(rev(grp), collapse = " vs ")
  path <- file.path("diff_analysis_result.xlsx")
  if (!file.exists(path)) next
  df <- read_excel(path)
  df$Comparison <- c
  a_list[[c]] <- df
}
a_df <- bind_rows(a_list)

a_df$Family <- brain_family[a_df$Cell_Type]
a_df$Star <- fdr_star(a_df$FDR)
a_df$Sig <- a_df$FDR < 0.05 & abs(a_df$Log2FC) > 1
a_df$Direction <- ifelse(a_df$Log2FC > 0, "Up", "Down")
a_df$Direction[!a_df$Sig] <- "NS"

cs_l2fc <- a_df %>%
  filter(Comparison == "LTS vs HC") %>%
  dplyr::select(Cell_Type, Log2FC) %>%
  rename(CS_l2fc = Log2FC)
a_df <- a_df %>% left_join(cs_l2fc, by = "Cell_Type")
a_df <- a_df %>% mutate(Family = factor(Family, levels = fam_order))
a_df <- a_df %>% arrange(Family, desc(CS_l2fc))
cell_order <- unique(a_df$Cell_Type)
a_df$Cell_Type <- factor(a_df$Cell_Type, levels = cell_order)
a_df$Comparison <- factor(a_df$Comparison, levels = comps)

p_main <- ggplot(a_df, aes(x = Comparison, y = Cell_Type)) +
  geom_point(aes(size = abs(Log2FC), color = Direction, alpha = ifelse(Sig, 1, 0.35))) +
  scale_color_manual(
    values = c("Up" = "#E24B4A", "Down" = "#378ADD", "NS" = "#B4B2A9"),
    name = "Direction",
    guide = guide_legend(override.aes = list(size = 4))
  ) +
  scale_size_continuous(range = c(1, 7), name = "|Log2FC|",limits = c(0, 8) ) +
  scale_alpha_identity() +
  labs(x = "", y = "", title = "A-DFC Brain Reference") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(linewidth = 0.2),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(5, 2, 5, 5)
  )
fam_anno <- a_df %>%
  distinct(Cell_Type, Family) %>%
  mutate(
    Cell_Type = factor(Cell_Type, levels = cell_order),
    Family = factor(Family, levels = fam_order)
  ) %>%
  arrange(Cell_Type)

p_fam <- ggplot(fam_anno, aes(x = 1, y = Cell_Type, fill = Family)) +
  geom_tile(width = 1, height = 0.9, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = fam_col, name = "Family") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0.5)) +
  labs(x = "", y = "") +
  theme_void(base_size = 8) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 0)
  )

pA <- p_main + p_fam + plot_layout(widths = c(10, 0.5))
ggsave("PA.png", pA, width = 9, height = 7, dpi = 300)
ggsave("PA.pdf", pA, width = 9, height = 7, dpi = 300)


# ============================================================================
# Figure 4B
# ============================================================================

fdr_star <- function(p) ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
comps <- c("SCZ vs HC","FES vs HC","MCS vs HC","LTS vs HC","LTS-MDI vs HC","LTS-nMDI vs HC", "LTS-MDI vs LTS-nMDI")
comp_groups <- list(
  "SCZ vs HC"=c("SCZ","HC"),
  "FES vs HC"=c("FES","HC"),
  "MCS vs HC"=c("MCS","HC"),
  "LTS vs HC"=c("LTS","HC"),
  "LTS-MDI vs HC"=c("LTS-MDI","HC"),
  "LTS-nMDI vs HC"=c("LTS-nMDI","HC"),
  "LTS-MDI vs LTS-nMDI"=c("LTS-LTS-MDI","LTS-nMDI")
)
comp_col <- c(
  "SCZ vs HC"="#c82423",
  "FES vs HC"="#E31A1C",
  "MCS vs HC"="#FF7F00",
  "LTS vs HC"="#33A02C",
  "LTS-MDI vs HC"="#1F78B4",
  "LTS-nMDI vs HC"="#984EA3",
  "LTS-MDI vs LTS-nMDI"="#FFED6F"
)
brain_family <- c(
  "intermediate monocyte" = "Myeloid/DC",
  "CD141-positive myeloid dendritic cell" = "Myeloid/DC",
  "transitional stage B cell" = "B/plasma",
  "naive B cell" = "B/plasma",
  "Be cell" = "B/plasma",
  "plasma cell" = "B/plasma",
  "CD16-positive, CD56-dim natural killer cell, human" = "NK/innate-like",
  "CD16-negative, CD56-bright natural killer cell, human" = "NK/innate-like",
  "mucosal invariant T cell" = "NK/innate-like",
  "mature gamma-delta T cell" = "NK/innate-like",
  "naive thymus-derived CD4-positive, alpha-beta T cell" = "CD4/Treg",
  "CD4-positive, alpha-beta memory T cell" = "CD4/Treg",
  "effector memory CD4-positive, alpha-beta T cell" = "CD4/Treg",
  "central memory CD4-positive, alpha-beta T cell" = "CD4/Treg",
  "T cell" = "CD4/Treg",
  "naive regulatory T cell" = "CD4/Treg",
  "memory regulatory T cell" = "CD4/Treg",
  "naive thymus-derived CD8-positive, alpha-beta T cell" = "CD8",
  "CD8-positive, alpha-beta memory T cell" = "CD8",
  "effector memory CD8-positive, alpha-beta T cell" = "CD8",
  "central memory CD8-positive, alpha-beta T cell" = "CD8",
  "CD8-positive, alpha-beta regulatory T cell" = "CD8",
  "CD8aa(I) thymocyte" = "CD8"
  
)

fam_col <- c("CD8"="#88D8D9","CD4/Treg"="#82C4B6","NK/innate-like"="#FAF39B","B/plasma"="#F1B56D","Myeloid/DC"="#FBB4AE")
fam_order <- c("CD8","CD4/Treg","NK/innate-like","B/plasma","Myeloid/DC")

a_list <- list()
for (c in comps) {
  grp <- strsplit(c, " vs ")[[1]]
  folder_name <- paste(rev(grp), collapse = " vs ")
  path <- file.path("diff_analysis_result.xlsx")
  if (!file.exists(path)) next
  df <- read_excel(path)
  df$Comparison <- c
  a_list[[c]] <- df
}
a_df <- bind_rows(a_list)

valid_cell <- names(brain_family)

a_df <- a_df %>%
  filter(Cell_Type %in% valid_cell)

a_df$Family <- brain_family[a_df$Cell_Type]
a_df$Star <- fdr_star(a_df$FDR)
a_df$Sig <- a_df$FDR < 0.05 & abs(a_df$Log2FC) > 1
a_df$Direction <- ifelse(a_df$Log2FC > 0, "Up", "Down")
a_df$Direction[!a_df$Sig] <- "NS"
a_df <- subset(a_df,Cell_Type!="erythrocyte")

cs_l2fc <- a_df %>%
  filter(Comparison == "LTS vs HC") %>%
  select(Cell_Type, Log2FC) %>%
  rename(CS_l2fc = Log2FC)
a_df <- a_df %>% left_join(cs_l2fc, by = "Cell_Type")
a_df <- a_df %>% mutate(Family = factor(Family, levels = fam_order))
a_df <- a_df %>% arrange(Family, desc(CS_l2fc))
cell_order <- unique(a_df$Cell_Type)
a_df$Cell_Type <- factor(a_df$Cell_Type, levels = cell_order)
a_df$Comparison <- factor(a_df$Comparison, levels = comps)

p_main <- ggplot(a_df, aes(x = Comparison, y = Cell_Type)) +
  geom_point(aes(size = abs(Log2FC), color = Direction, alpha = ifelse(Sig, 1, 0.35))) +
  scale_color_manual(
    values = c("Up" = "#E24B4A", "Down" = "#378ADD", "NS" = "#B4B2A9"),
    name = "Direction",
    guide = guide_legend(override.aes = list(size = 4))
  ) +
  scale_size_continuous(range = c(1, 7),name = "|Log2FC|", limits = c(0, 8)) +
  scale_alpha_identity() +
  labs(x = "", y = "", title = "B-Peripheral Blood Reference") +
  theme_bw(base_size = 10) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 10),
    axis.text.y = element_text(size = 10),
    panel.grid.major = element_line(linewidth = 0.2),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 12),
    plot.margin = margin(5, 2, 5, 5)
  )

fam_anno <- a_df %>%
  distinct(Cell_Type, Family) %>%
  mutate(
    Cell_Type = factor(Cell_Type, levels = cell_order),
    Family = factor(Family, levels = fam_order)
  ) %>%
  arrange(Cell_Type)

p_fam <- ggplot(fam_anno, aes(x = 1, y = Cell_Type, fill = Family)) +
  geom_tile(width = 1, height = 0.9, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = fam_col, name = "Family") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0.5)) +
  labs(x = "", y = "") +
  theme_void(base_size = 8) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right",
    plot.margin = margin(5, 5, 5, 0)
  )
pB <- p_main + p_fam + plot_layout(widths = c(10, 0.5))
ggsave("PB.png", pB, width = 9, height = 7, dpi = 300)
ggsave("PB.pdf", pB, width = 9, height = 7, dpi = 300)


# ============================================================================
# Figure  4C 
# ============================================================================
suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

data_dir <- "Figure4"
brain_file <- file.path(data_dir, "brain_cell_type.txt")
blood_file <- file.path(data_dir, "blood_cell_type.txt")
meta_file <- file.path(data_dir, "metadata.txt")
OUT <- file.path(data_dir)
meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   comment.char = "", stringsAsFactors = FALSE)
gm <- setNames(meta$Group, meta$ID)
mat_b <- read.csv(brain_file, header = TRUE, check.names = FALSE, row.names = 1)

cbrain <- c(
  "oligodendrocyte",
  "VIP GABAergic cortical interneuron",
  "astrocyte of the cerebral cortex",
  "L5 extratelencephalic projecting glutamatergic cortical neuron"
)
cbrain <- cbrain[cbrain %in% colnames(mat_b)]

brain_long <- as.data.frame(mat_b[, cbrain, drop = FALSE])
brain_long$Sample <- rownames(brain_long)
brain_long <- pivot_longer(brain_long, -Sample, names_to = "CellType", values_to = "Proportion")
brain_long$Source <- "Brain"

all_long <- brain_long

all_long$Group <- gm[all_long$Sample]
all_long <- all_long %>% filter(Group %in% c("HC", "LTS-MDI"))
all_long$Group <- factor(all_long$Group, levels = c("HC", "LTS-MDI"))

all_long <- all_long %>%
  mutate(CellType_short = case_when(
    CellType == "oligodendrocyte" ~ "Oligodendrocyte",
    CellType == "VIP GABAergic cortical interneuron" ~ "VIP GABA-IN",
    CellType == "astrocyte of the cerebral cortex" ~ "Cortical Astrocyte",
    CellType == "L5 extratelencephalic projecting glutamatergic cortical neuron" ~ "L5 ET Glu-Neuron",
    TRUE ~ CellType
  ))

ct_order <- c( "Oligodendrocyte", "VIP GABA-IN", "Cortical Astrocyte","L5 ET Glu-Neuron")
all_long$CellType_short <- factor(all_long$CellType_short, levels = ct_order)

comp_col <- c("HC"= "#1F78B4", "LTS-MDI" = "#E31A1C")

pC <- ggplot(all_long, aes(x = Group, y = Proportion, fill = Group)) +
  geom_violin(trim = TRUE, alpha = 0.4, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.08, alpha = 0.15, size = 0.4) +
  facet_wrap(~ CellType_short, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = comp_col) +
  labs(x = "", y = "inferred fraction") +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
    axis.text.x = element_text(angle = 0, size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )

out_png_local <- file.path(data_dir, "PC.png")
out_pdf_local <- file.path(data_dir, "PC.pdf")
ggsave(out_png_local, pC, width = 10, height = 5, dpi = 300)
ggsave(out_pdf_local, pC, width = 10, height = 5, dpi = 300)


# ============================================================================
# Figure 4D
# ============================================================================
data_dir <- "Figure4"
brain_file <- file.path(data_dir, "brain_cell_type.txt")
blood_file <- file.path(data_dir, "blood_cell_type.txt")
meta_file <- file.path(data_dir, "metadata.txt")

OUT <- file.path(data_dir)
meta <- read.table(meta_file, header = TRUE, sep = "\t",
                   comment.char = "", stringsAsFactors = FALSE)
gm <- setNames(meta$Group, meta$ID)
mat_bl <- read.csv(blood_file, header = TRUE, check.names = FALSE, row.names = 1)

cblood <- c(
  "transitional stage B cell",
  "CD4-positive, alpha-beta memory T cell",
  "central memory CD8-positive, alpha-beta T cell",
  "naive thymus-derived CD8-positive, alpha-beta T cell"
)
cblood <- cblood[cblood %in% colnames(mat_bl)]

blood_long <- as.data.frame(mat_bl[, cblood, drop = FALSE])
blood_long$Sample <- rownames(blood_long)
blood_long <- pivot_longer(blood_long, -Sample, names_to = "CellType", values_to = "Proportion")
blood_long$Source <- "Blood"

all_long <- blood_long

all_long$Group <- gm[all_long$Sample]
all_long <- all_long %>% filter(Group %in% c("HC", "LTS-MDI"))
all_long$Group <- factor(all_long$Group, levels = c("HC", "LTS-MDI"))
all_long <- all_long %>%
  mutate(CellType_short = case_when(
    CellType == "transitional stage B cell" ~ "Transitional B",
    CellType == "CD4-positive, alpha-beta memory T cell" ~ "CD4 Tmem",
    CellType == "central memory CD8-positive, alpha-beta T cell" ~ "CD8 Tcm",
    CellType == "naive thymus-derived CD8-positive, alpha-beta T cell" ~ "CD8 naive T",
    TRUE ~ CellType
  ))

ct_order <- c( "Transitional B", "CD4 Tmem","CD8 Tcm", "CD8 naive T")

all_long$CellType_short <- factor(all_long$CellType_short, levels = ct_order)

comp_col <- c("HC"= "#1F78B4", "LTS-MDI" = "#E31A1C")


pD <- ggplot(all_long, aes(x = Group, y = Proportion, fill = Group)) +
  geom_violin(trim = TRUE, alpha = 0.4, scale = "width") +
  geom_boxplot(width = 0.15, outlier.shape = NA, alpha = 0.8) +
  geom_jitter(width = 0.08, alpha = 0.15, size = 0.4) +
  facet_wrap(~ CellType_short, scales = "free_y", nrow = 1) +
  scale_fill_manual(values = comp_col) +
  labs(x = "", y = "inferred fraction") +
  theme_bw(base_size = 14) +
  theme(
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "white", colour = "black", linewidth = 0.6),
    axis.text.x = element_text(angle = 0, size = 14),
    axis.text.y = element_text(size = 14),
    legend.position = "right",
    plot.title = element_text(face = "bold", size = 14),
    panel.grid.minor = element_blank()
  )
out_png_local <- file.path(data_dir, "PD.png")
out_pdf_local <- file.path(data_dir, "PD.pdf")

ggsave(out_png_local, pD, width = 10, height = 5, dpi = 300)
ggsave(out_pdf_local, pD, width = 10, height = 5, dpi = 300)
