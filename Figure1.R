
############################################################
## Figure 1A
############################################################
library(tidyverse)
library(ggpubr)
library(FactoMineR)
library(factoextra)
library(edgeR)
library(rstatix)

pl_list=c("HC  SCZ")
for(pl in pl_list){
  pl2 <- as.character(unlist(strsplit(pl, split = "  ")))
  metadata <- read.table(paste0("metadata.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE)
  metadata <- subset(metadata, Group %in% pl2)
  metadata$Group <- factor(metadata$Group,levels = pl2) 

  expr_count <- read.table("expr_count.tsv", header = TRUE, row.names = 1, sep = "\t", check.names = FALSE)
  expr_count <- expr_count[, metadata$ID]

  count_matrix <- round(as.matrix(expr_count))

  dds <- DESeqDataSetFromMatrix(
    countData = count_matrix,
    colData = metadata,
    design = ~ Group
  )
  vsd <- vst(dds, blind = TRUE)
  vst_matrix <- assay(vsd)
  
  set.seed(123) 
  permanova_res <- adonis2(t(vst_matrix) ~ Group, data = metadata, method = "euclidean", permutations = 999)
  p_val <- permanova_res$`Pr(>F)`[1]
  r_squared <- permanova_res$R2[1]
  if (p_val < 0.001) {
    p_label <- "P < 0.001"
  } else {
    p_label <- sprintf("P = %.3f", p_val)
  }
  anno_text <- sprintf("PERMANOVA\n%s\nR² = %.3f", p_label, r_squared)
  pca_expr <- prcomp(t(vst_matrix), scale. = FALSE)
  
  p_pca_expr <- fviz_pca_ind(
    pca_expr,
    geom = "point",
    habillage = metadata$Group,
    palette = c("#2878b5", "#c82423"),
    addEllipses = TRUE,        
    ellipse.level = 0.95,
    title = "PCA of cfRNA Gene Expression"
  ) +
    theme_bw() +
    annotate(
      "text", 
      x = Inf, y = Inf,         
      label = anno_text, 
      hjust = 1.1, vjust = 1.5,  
      size = 4,                 
      fontface = "bold",       
      color = "black"
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold") # 标题居中
    )
  print(p_pca_expr)
  ggsave(paste0(pl,"_Gene_Expression_PCA.pdf"), p_pca_expr, width = 5, height = 4)
  ggsave(paste0(pl,"_Gene_Expression_PCA.png"), p_pca_expr, width = 5, height = 4)
  write.csv(pca_expr$x, paste0(pl,"_Gene_Expression_PCA_Coordinates.csv"))
}  
  
############################################################
## Figure 1B-left
############################################################
library(ggplot2)
library(ggrepel)

Res <- read.table("Differential_analysis_results.csv",header=T, sep=",", comment.char="",check.names=F,quot="",stringsAsFactors = F)
colnames(Res)[1] <- "Gene"
feature <- read.table("feature.tsv", header=T, sep="\t", comment.char="", check.names=F, quote="", stringsAsFactors=F) %>%
  .[, c(1,2)] %>%
  `colnames<-`(c("Gene", colnames(.)[2])) %>%
  {`[<-`(., is.na(.[,2]), 2, value = .[is.na(.[,2]), 1])}

Res <- left_join(Res,feature,by="Gene")

log2FoldChange_cutoff <- log2(2)
p_cutoff <- 0.05

out <-Res %>%
  mutate(change = case_when(
    (padj < p_cutoff & log2FoldChange > log2FoldChange_cutoff) ~"up",  
    (padj < p_cutoff & log2FoldChange < -log2FoldChange_cutoff) ~"down",
    .default ="not sig"
  ))

outSorted = out[!is.na(out$padj), ]
outSorted = outSorted[order(outSorted$padj), ]

out_up <- subset(outSorted, change =="up")
out_down <- subset(outSorted, change =="down")
out_diff <- rbind(out_up,out_down)
DEG <- outSorted

DEG$change <- factor(DEG$change, levels =  c("up", "not sig", "down"))
max_lfc <- max(DEG$log2FoldChange, na.rm=T)
min_lfc <- min(DEG$log2FoldChange, na.rm=T)

mark_gene <- c("HBA2","","S100A9","","AQP9","TMEM88B","VTRNA1-2", "FCGR3B")
deg_vol_plot <- ggplot(data = DEG,
                       aes(x = log2FoldChange,
                           y = -log10(padj), 
                           colour = change)) +
  scale_color_manual(values = c("#c82423", "#7b7c7d","#2878b5")) +
  geom_point(size = 1.5, alpha = 0.5, na.rm=T) +

geom_text_repel(
  data = subset(DEG, gene_name %in% mark_gene), 
  aes(label = gene_name),
  size = 3.5,             
  force = 2,              
  min.segment.length = 0, 
  segment.size = 0.4,     
  segment.color = "black",
  box.padding = 0.4,      
  point.padding = 0.2,    
  max.overlaps = Inf,     
  fontface = "bold",
  show.legend = FALSE
) +
scale_x_continuous(limits = c(min_lfc, max_lfc)) +
  geom_hline(yintercept = -log10(0.05),lty = 4,col = "darkgray",lwd = 0.6) +
  geom_vline(xintercept = c(-log2FoldChange_cutoff, log2FoldChange_cutoff),lty = 4,col = "darkgray",lwd = 0.6) +
  theme_bw(base_size = 12) +
  theme(legend.justification = c(0,1),
        legend.position = c(0,1),
        legend.background = element_rect(fill = "white", color = "black", size = 0.2),
        panel.grid = element_blank(),
        legend.title = element_blank(),
        legend.text = element_text(face="bold",color="black",size=12),
        plot.title = element_text(hjust = 0.5, face = "bold",color = "black",size = 14),
        axis.text.x = element_text(face = "bold",color = "black",size = 12),
        axis.text.y = element_text(face = "bold",color = "black",size = 12),
        axis.title.x = element_text(face = "bold",color = "black",size = 12),
        axis.title.y = element_text(face = "bold",color = "black",size = 12),
        plot.subtitle = element_text(hjust = 0.5,size = 12, face = "italic", colour = "black")) +
  labs(title = paste0("CD vs  CS"),
       x=expression(paste("log"[2], "(FoldChange)")), 
       y=expression(paste("-log"[10], "(padj)")),
       subtitle = paste(sprintf('padj:  %.2f;', 0.05),
                        sprintf("log2FC: %.2f;", log2FoldChange_cutoff),
                        sprintf('Up: %1.0f; Down: %1.0f;', nrow(out_up), nrow(out_down)),
                        sprintf('Total: %1.0f', nrow(out_diff))))
print(deg_vol_plot)
ggsave(filename = paste0("volcano_label.pdf"), deg_vol_plot,width=150, height=150, units="mm")
ggsave(filename = paste0("volcano_label.png"), deg_vol_plot, width=150, height=150, units="mm")

############################################################
## Figure 1B-right
############################################################
library(ggplot2)
library(dplyr)
TMPDIR <- "Figure1"
OUTDIR <- TMPDIR

COMPARISONS <- c("HC__SCZ")
COMP_LABELS <- c("HC  SCZ")
comparison_config <- list("HC__SCZ" = list(col = "Group", ctrl = "HC", case = "SCZ"))

module_display <- c(
  "blue" = "Blue", "brown" = "Brown", "turquoise" = "Turquoise", "yellow" = "Yellow"
)
meta_all <- read.table(file.path(TMPDIR, "metadata.txt"), header = TRUE, sep = "\t", stringsAsFactors = FALSE, comment.char = "", quote = "")

all_results <- list()
for (i in seq_along(COMPARISONS)) {
  comp <- COMPARISONS[i]
  comp_label <- COMP_LABELS[i]
  cfg <- comparison_config[[comp]]
  cat("Processing:", comp_label, "\n")
  eigengene_path <- file.path(paste0("02_Module_Eigengenes.txt"))
  if (!file.exists(eigengene_path)) {
    cat("  File not found:", eigengene_path, "\n")
    next
  }
  
  MEs <- read.table(eigengene_path, header = TRUE, sep = "\t", 
                    stringsAsFactors = FALSE, check.names = FALSE)
  MEs$ID <- rownames(MEs)
  comp_meta <- meta_all[meta_all[[cfg$col]] %in% c(cfg$ctrl, cfg$case), ]
  common_samples <- intersect(comp_meta$ID, MEs$ID)
  if (length(common_samples) < 10) {
    cat("  Too few common samples:", length(common_samples), "\n")
    next
  }
  MEs <- MEs[MEs$ID %in% common_samples, ]
  comp_meta <- comp_meta[comp_meta$ID %in% common_samples, ]
  
  MEs <- MEs[order(MEs$ID), ]
  comp_meta <- comp_meta[order(comp_meta$ID), ]
  reg_data <- merge(MEs, comp_meta, by = "ID")
  reg_data$Group <- factor(reg_data$Group, levels = c(cfg$ctrl, cfg$case))
  reg_data$age <- as.numeric(reg_data$age)
  reg_data$BMI <- as.numeric(reg_data$BMI)
  reg_data$gender <- factor(reg_data$gender)
  reg_data$center <- factor(reg_data$center)
  n1 <- sum(reg_data$Group == cfg$ctrl)
  n2 <- sum(reg_data$Group == cfg$case)
  module_names <- setdiff(colnames(MEs), "ID")
  group_term <- paste0("Group", cfg$case) 
  for (mod in module_names) {
    mod_color <- gsub("^ME", "", mod)
    if (n1 < 3 || n2 < 3) next
    fit <- tryCatch({
      lm(as.formula(paste(mod, "~ Group + age + gender + BMI + center")), data = reg_data)
    }, error = function(e) NULL)
    
    if (is.null(fit)) next
    if (!group_term %in% rownames(summary(fit)$coefficients)) next
    
    coef_tab <- summary(fit)$coefficients
    ci_tab <- confint(fit)
    ME_Diff <- coef_tab[group_term, "Estimate"]
    ME_Diff_lower <- ci_tab[group_term, 1]
    ME_Diff_upper <- ci_tab[group_term, 2]
    P_value <- coef_tab[group_term, "Pr(>|t|)"]

    m1 <- mean(reg_data[reg_data$Group == cfg$ctrl, mod], na.rm = TRUE)
    m2 <- mean(reg_data[reg_data$Group == cfg$case, mod], na.rm = TRUE)
    s1 <- sd(reg_data[reg_data$Group == cfg$ctrl, mod], na.rm = TRUE)
    s2 <- sd(reg_data[reg_data$Group == cfg$case, mod], na.rm = TRUE)
    sp <- sqrt(((n1 - 1) * s1^2 + (n2 - 1) * s2^2) / (n1 + n2 - 2))
    Cohen_d <- (m2 - m1) / sp
    
    all_results[[length(all_results) + 1]] <- data.frame(
      Comparison = comp_label,
      Module = mod,
      Module_Color = mod_color,
      Module_Display = ifelse(mod_color %in% names(module_display), 
                              module_display[mod_color], mod_color),
      n_Control = n1,
      n_Disease = n2,
      Mean_Control = m1,
      Mean_Disease = m2,
      ME_Diff = ME_Diff,                    
      ME_Diff_lower = ME_Diff_lower,  
      ME_Diff_upper = ME_Diff_upper,  
      Cohen_d = Cohen_d,
      P_value = P_value,
      stringsAsFactors = FALSE
    )
  }
}
results_df <- do.call(rbind, all_results)
results_df$FDR <- p.adjust(results_df$P_value, method = "fdr")
results_df <- results_df %>%
  mutate(
    Diff_CI_label = sprintf("%.3f (%.3f to %.3f)", ME_Diff, ME_Diff_lower, ME_Diff_upper),
    FDR_label = ifelse(FDR < 0.001, sprintf("FDR = %.2e", FDR), sprintf("FDR = %.3f", FDR))
  )

write.csv(results_df, file.path(OUTDIR, "Forest_plot_data_FDR.csv"), row.names = FALSE)

results_df <- results_df[results_df$Module_Display %in% c("Blue", "Brown", "Turquoise", "Yellow"), ]
results_df$Module_Display <- factor(results_df$Module_Display, 
                                    levels = c("Yellow", "Turquoise", "Brown", "Blue"))

module_real_colors <- c("Blue" = "#2878B5", "Brown" = "#A0522D", 
                        "Turquoise" = "turquoise", "Yellow" = "#FFD700")

x_range <- max(results_df$ME_Diff_upper) - min(results_df$ME_Diff_lower)
text_x_pos <- max(results_df$ME_Diff_upper) + x_range * 0.15 

p1 <- ggplot(results_df, aes(x = ME_Diff, y = Module_Display)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey60", linewidth = 0.6) +
  geom_errorbarh(aes(xmin = ME_Diff_lower, xmax = ME_Diff_upper, color = Module_Display), 
                 height = 0, linewidth = 1.2) +
  geom_point(data = subset(results_df, FDR >= 0.05), 
             aes(color = Module_Display), fill = "white", shape = 21, size = 5, stroke = 1.2) +
  geom_point(data = subset(results_df, FDR < 0.05), 
             aes(color = Module_Display), shape = 16, size = 5) +
  geom_text(aes(x = text_x_pos, label = Diff_CI_label), 
            hjust = 0, vjust = -0.5, color = "grey20", size = 3.5) +
  geom_text(aes(x = text_x_pos, label = FDR_label), 
            hjust = 0, vjust = 1.5, color = "grey50", size = 3.2) +
  annotate("text", x = text_x_pos, y = length(levels(results_df$Module_Display)) + 0.6, 
           label = "Adjusted difference (95% CI)", hjust = 0, fontface = "bold", size = 3.5, color = "grey20") +
  scale_color_manual(values = module_real_colors) +
  scale_x_continuous(expand = expansion(mult = c(0.1, 0.7))) +
  scale_y_discrete(expand = expansion(mult = c(0.1, 0.2))) +
  coord_cartesian(clip = "off") +
  labs(title = "Coexpression-module changes",
       subtitle = "Adjusted for age, sex, BMI, and center",
       x = "Adjusted eigengene difference (SCZ minus HC)",
       y = NULL) +
  theme_minimal(base_size = 12) +
  theme(
    panel.grid = element_blank(),
    axis.line.x = element_line(color = "grey60"),
    axis.line.y = element_line(color = "grey80", linewidth = 0.6),
    axis.ticks.x = element_line(color = "grey60"),
    axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 12, color = "grey20", margin = margin(r = 15)),
    axis.text.x = element_text(size = 11, color = "grey20"),
    axis.title.x = element_text(size = 11, color = "grey20", margin = margin(t = 15)),
    plot.title = element_text(face = "bold", size = 15, hjust = 0),
    plot.subtitle = element_text(size = 11, color = "grey50", margin = margin(b = 20)),
    legend.position = "none",
    plot.margin = margin(t = 20, r = 20, b = 20, l = 10)
  )

ggsave(file.path(OUTDIR, "Figure_WGCNA_Forest_Plot_FDR.png"),  p1, width = 4, height = 5, dpi = 300)
ggsave(file.path(OUTDIR, "Figure_WGCNA_Forest_Plot_FDR.pdf"),  p1, width = 4, height = 5)

############################################################
## Figure 1C
############################################################
library(readxl)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(cowplot)

TMPDIR <- "Figure1"
OUTDIR <- TMPDIR 
raw <- read_excel("SCZ_HC_GSEA_result.xlsx", sheet = "Sheet2")

colnames(raw)[1] <- "Source"
colnames(raw)[4] <- "Pathway_type"

cat("Total rows loaded:", nrow(raw), "\n")
print(raw[, c("Source", "Description", "Pathway_type", "NES", "p.adjust")], row.names = FALSE)

raw$Source <- gsub(" module$", "", raw$Source)  # "Blue module" -> "Blue"

type_mapping <- c(
  "Intracellular trafficking"              = "Intracellular trafficking",
  "Immune/Inflammatory"               = "Immune/Inflammatory",
  "Neuronal/Synaptic"                = "Neuronal/Synaptic",
  "Protein synthesis/Proteostasis"     = "Protein synthesis/Proteostasis",
  "Mitochondrial/Energy"               = "Mitochondrial/Energy",
  "Cell signaling"          = "Cell signaling",
  "Cytoskeletal"                       = "Cytoskeletal",
  "Platelet"                         = "Platelet"
)

raw$Pathway_type <- type_mapping[raw$Pathway_type]
raw$Pathway_type[is.na(raw$Pathway_type)] <- "Other"
raw$significant <- as.numeric(raw$p.adjust) < 0.05
raw$fill_val <- ifelse(raw$significant, as.numeric(raw$NES), NA_real_)
raw$label <- ifelse(raw$significant, sprintf("%.2f", as.numeric(raw$NES)), "")
sources <- c("Normal", "Blue")
raw$Source <- factor(raw$Source, levels = sources)
pathway_order <- c(
  "neutrophil chemotaxis",
  "chemokine activity",
  "MHC class II protein complex",
  "C-C chemokine receptor activity",
  "regulation of dendritic spine development",
  "proton-transporting ATP synthase complex",
  "U2-type catalytic step 2 spliceosome",
  "translation initiation factor activity",
  "late endosome to lysosome transport",
  "integrin-mediated signaling pathway",
  "contractile actin filament bundle",
  "platelet alpha granule lumen",
  "platelet aggregation"
)
theme_order <- c("Immune/Inflammatory", "Neuronal/Synaptic",
                 "Protein synthesis/Proteostasis",
                 "Mitochondrial/Energy",
                 "Intracellular trafficking",
                 "Cell signaling","Cytoskeletal",
                 "Platelet")
theme_map <- raw %>%
  distinct(Description, Pathway_type) %>%
  { setNames(.$Pathway_type, .$Description) }

raw <- raw %>%
  mutate(theme = Pathway_type,
         pathway = Description)
raw$pathway <- factor(raw$pathway, levels = rev(pathway_order))
raw$theme <- factor(raw$theme, levels = theme_order)

write.csv(raw, file.path(TMPDIR, "GSEA_heatmap_data.csv"), row.names = FALSE)

source_colors <- c(
  "Normal" = "#DECBE4",
  "Blue"   = "#A6CEE3"
)
theme_colors <- c(
  "Immune/Inflammatory"                = "#c82423",
  "Neuronal/Synaptic"                  = "#FBB4AE",
  "Protein synthesis/Proteostasis" = "#71B8EC",
  "Mitochondrial/Energy"    = "#82C4B6",
  "Intracellular trafficking"            = "#F1B56D",
  "Cell signaling"               = "#CAB2D6",
  "Cytoskeletal"    = "#B8AEEA",
  "Platelet" = "#F2A8D9"
)

nes_lim <- max(abs(as.numeric(raw$NES)), na.rm = TRUE)
nes_lim <- ceiling(nes_lim * 10) / 10

strip_data <- raw %>% distinct(pathway, theme)

p_left_strip <- ggplot(strip_data, aes(x = 1, y = pathway, fill = theme)) +
  geom_tile(color = "white", linewidth = 0.5) +
  scale_fill_manual(values = theme_colors, guide = "none") +
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_discrete(limits = levels(raw$pathway), drop = FALSE) +
  theme_void() +
  theme(plot.margin = margin(2, 0, 2, 2))

p_main <- ggplot(raw, aes(x = Source, y = pathway)) +
  geom_tile(aes(fill = fill_val), color = "white", linewidth = 0.4) +
  geom_text(aes(label = label), size = 3.8, color = "black", fontface = "bold") +
  scale_fill_gradient2(
    name = "NES",
    low = "#2878b5", mid = "#F7F7F7", high = "#c82423",
    midpoint = 0, limits = c(-nes_lim, nes_lim),
    na.value = "#F0F0F0"
  ) +
  scale_x_discrete(position = "top", drop = FALSE) +
  scale_y_discrete(limits = levels(raw$pathway), position = "right", drop = FALSE) +
  guides(fill = guide_colourbar(order = 1, barheight = unit(35, "mm"),
                                barwidth = unit(3, "mm"))) +
  theme_minimal(base_size = 12, base_family = "sans") +
  theme(
    axis.text.x.top = element_text(size = 12, color = "black", face = "bold", 
                                   hjust = 0.5, vjust = 0.5),
    axis.ticks.x.top = element_blank(),
    axis.text.y.right = element_text(size = 12, colour = "black", hjust = 0),
    axis.title = element_blank(),
    panel.grid = element_blank(),
    legend.position = "right",
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 10),
    plot.margin = margin(2, 2, 2, 2)
  )

top_tile_data <- data.frame(
  Source = factor(sources, levels = sources),
  dummy = "Top"
)

p_top_tile <- ggplot(top_tile_data, aes(x = Source, y = dummy, fill = Source)) +
  geom_tile(color = "white", linewidth = 0.4) +
  scale_fill_manual(values = source_colors, guide = "none") +
  scale_x_discrete(position = "top", drop = FALSE) +
  theme_void() +
  theme(plot.margin = margin(0, 2, 0, 2))

p_main_noleg <- p_main + theme(legend.position = "none")

legend_nes <- get_legend(p_main)

p_src_legend <- ggplot(top_tile_data, aes(x = Source, y = dummy, fill = Source)) +
  geom_tile() +
  scale_fill_manual(values = source_colors, name = "Source") +
  theme_minimal() +
  theme(legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))
legend_src <- get_legend(p_src_legend)

p_theme_legend <- ggplot(strip_data, aes(x = theme, y = 1, fill = theme)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_manual(values = theme_colors, name = "Pathway type") +
  theme_void() +
  theme(legend.position = "right",
        legend.title = element_text(size = 12, face = "bold"),
        legend.text = element_text(size = 10))
legend_theme <- get_legend(p_theme_legend)

p_right_leg <- plot_grid(legend_src, legend_nes, legend_theme,
                         ncol = 1, align = "v", axis = "l",
                         rel_heights = c(1, 1.2, 1.5))
left_block <- plot_grid(p_left_strip, p_main_noleg, ncol = 2, 
                        rel_widths = c(0.4, 8.5), align = "h")

left_col <- plot_grid(p_top_tile, left_block, ncol = 1,
                      rel_heights = c(0.35, 9.5), align = "v")

p_final <- plot_grid(left_col, p_right_leg, ncol = 2, rel_widths = c(8.5, 2))

w_in <- 200 / 25.4
h_in <- 140 / 25.4

ggsave(file.path(OUTDIR, "pathway_heatmap.png"),p_final, width = w_in, height = h_in, dpi = 300)
ggsave(file.path(OUTDIR, "pathway_heatmap.pdf"),p_final, width = w_in, height = h_in)


############################################################
## Figure 1D right
############################################################

total_train_stats <- paste0(output_path, "/", SAMP_GROUP, "_all_modles.train.tsv")
total_test_stats <- paste0(output_path, "/", SAMP_GROUP, "_all_models.test.tsv")
total_val_stats <- paste0(output_path, "/", SAMP_GROUP, "_all_models.val.tsv")

if(file.exists(total_train_stats) && file.exists(total_test_stats) && file.exists(total_val_stats)){
  
  total_train_outputs <- read.table(total_train_stats, sep="\t", header=T, row.names=1, check.names = FALSE)
  total_test_outputs <- read.table(total_test_stats, sep="\t", header=T, row.names=1, check.names = FALSE)
  total_val_outputs <- read.table(total_val_stats, sep="\t", header=T, row.names=1, check.names = FALSE)
  train_auc_1 <- as.data.frame(t(total_train_outputs[1,]))
  test_auc_1 <- as.data.frame(t(total_test_outputs[1,]))
  val_auc_1 <- as.data.frame(t(total_val_outputs[1,]))
  
  df <- cbind(train_auc_1, test_auc_1, val_auc_1)
  colnames(df) <- c("train set", "test set", "vaild set")
  df$`train set` <- as.numeric(df$`train set`)
  df$`test set` <- as.numeric(df$`test set`)
  df$`vaild set` <- as.numeric(df$`vaild set`)
  plot_df <- df %>% arrange(desc(`vaild set`))
  clusterCols <- c("#D51F26", "#272E6A", "#208A42")
  names(clusterCols) <- c("train set","test set","vaild set")
  
  col_ha = columnAnnotation("Cohort" = colnames(plot_df), col = list("Cohort" = clusterCols), show_annotation_name = FALSE)

  plot_heatmap <- Heatmap(as.matrix(plot_df), name = "AUC",
                          top_annotation = col_ha,
                          col = c("#4195C1", "#FFFFFF", "#CB5746"), 
                          rect_gp = gpar(col = "black", lwd = 1), 
                          cluster_columns = FALSE, cluster_rows = FALSE, 
                          show_column_names = FALSE, 
                          show_row_names = TRUE, 
                          row_names_gp = gpar(fontsize = 11), 
                          row_names_side = "left",
                          column_split = factor(colnames(plot_df), levels = colnames(plot_df)), 
                          column_title = NULL,
                          
                          heatmap_legend_param = list(
                            title_gp = gpar(fontsize = 12, fontface = "bold"), 
                            labels_gp = gpar(fontsize = 10)
                          ),
                          
                          cell_fun = function(j, i, x, y, w, h, col) { 
                            grid.text(label = format(plot_df[i, j], digits = 3, nsmall = 3),
                                      x, y, gp = gpar(fontsize = 12, col="black")) 
                          }
  )
  
  total_heatmap_png <- paste0(output_path, "/", SAMP_GROUP, "_all_auc_heatmap_resized.png")
  png(total_heatmap_png, width = 6, height = 4.5, units = "in", res = 600)
  draw(plot_heatmap)
  dev.off()
  
  total_heatmap_pdf <- paste0(output_path, "/", SAMP_GROUP, "_all_auc_heatmap_resized.pdf")
  pdf(total_heatmap_pdf, width = 6, height = 4.5)
  draw(plot_heatmap)
  dev.off()
} 

plot_heatmap <- Heatmap(as.matrix(plot_df), name = "AUC",
                        top_annotation = col_ha,
                        col = c("#4195C1", "#FFFFFF", "#CB5746"), 
                        rect_gp = gpar(col = "black", lwd = 1), 
                        cluster_columns = FALSE, cluster_rows = FALSE, 
                        
                        show_column_names = FALSE, 
                        
                        show_row_names = TRUE, 
                        row_names_gp = gpar(fontsize = 11),    
                        row_names_side = "right",       # ✅ 纵轴模型名称放到右边
                        
                        column_split = factor(colnames(plot_df), levels = colnames(plot_df)), 
                        column_title = NULL,
                        
                        # AUC色条图例参数：横向摆放
                        heatmap_legend_param = list(
                          title_gp = gpar(fontsize = 12, fontface = "bold"), 
                          labels_gp = gpar(fontsize = 10),
                          legend_direction = "horizontal", # 图例横向
                          legend_width = unit(4, "cm")     # 横向图例长度
                        ),
                        
                        cell_fun = function(j, i, x, y, w, h, col) { 
                          grid.text(label = format(plot_df[i, j], digits = 3, nsmall = 3),
                                    x, y, gp = gpar(fontsize = 12, col="black")) 
                        }
)

plot_heatmap <- draw(plot_heatmap,
                     heatmap_legend_side = "left",
                     annotation_legend_side = "left",
                     merge_legend = TRUE,
                     padding = unit(c(0.8, 5, 0.8, 1.2), "cm"))


total_heatmap_png <- paste0(output_path, "/", SAMP_GROUP, "_all_auc_heatmap_resized_V2.png")
png(total_heatmap_png, width = 8, height = 4.5, units = "in", res = 600)
draw(plot_heatmap)
dev.off()

total_heatmap_pdf <- paste0(output_path, "/", SAMP_GROUP, "_all_auc_heatmap_resized_V2.pdf")
pdf(total_heatmap_pdf, width = 8, height = 4.5)
draw(plot_heatmap)
dev.off()

############################################################
## Figure 1D left
############################################################
MODELS <- c("C5","EXTRATREES","GLM","GLMNETLasso","GLMNETRidge","KNN","LDA","NB","NNET","PAM","RF","RPART","SVMLIN","SVMRAD")
count <- 1

for ( MODEL in MODELS ) {
  if(MODEL == "C5"){ MODEL_FIT = C5_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "EXTRATREES"){ MODEL_FIT = EXTRATREES_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "GLM"){ MODEL_FIT = GLM_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "GLMNETLasso"){ MODEL_FIT = GLMNETLasso_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "GLMNETRidge"){ MODEL_FIT = GLMNETRidge_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "KNN"){ MODEL_FIT = KNN_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "LDA"){ MODEL_FIT = LDA_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "NB"){ MODEL_FIT = NB_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "NNET"){ MODEL_FIT = NNET_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "PAM"){ MODEL_FIT = PAM_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "RF"){ MODEL_FIT = RF_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "RPART"){ MODEL_FIT = RPART_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "SVMLIN"){ MODEL_FIT = SVMLIN_TRAIN(count_matrix_train_sig,meta_data_train) }
  if(MODEL == "SVMRAD"){ MODEL_FIT = SVMRAD_TRAIN(count_matrix_train_sig,meta_data_train) }
  PERF_LIST = meaure_performance(MODEL_FIT,
                                 count_matrix_train_sig,count_matrix_test_sig,count_matrix_val_sig,
                                 meta_data_train,meta_data_test,meta_data_val,
                                 meta_data_all, count_matrix_all_sig,
                                 group1_name,group2_name,
                                 "prob",MODEL)
  
  MODEL_LIST = PERF_LIST[["MODEL_LIST"]]                      
  train_outputs = PERF_LIST[["train_outputs"]]
  test_outputs = PERF_LIST[["test_outputs"]]
  val_outputs = PERF_LIST[["val_outputs"]]
  AUC_figure = PERF_LIST[["AUC_figure"]]
  CS_figure = PERF_LIST[["CS_figure"]]
  
  rds_output <- paste(output_path,paste(SAMP_GROUP,"_",MODEL,".rds",sep = ""),sep = "/")
  train_stats <- paste(output_path,paste(SAMP_GROUP,"_",MODEL,".train.tsv",sep = ""),sep = "/")
  test_stats <- paste(output_path,paste(SAMP_GROUP,"_",MODEL,".test.tsv",sep = ""),sep = "/")
  val_stats <- paste(output_path,paste(SAMP_GROUP,"_",MODEL,".val.tsv",sep = ""),sep = "/")
  auc_figure <- paste(output_path,paste(SAMP_GROUP,"_",MODEL,".auc.png",sep = ""),sep = "/")
  classifier_figure <- paste(output_path,paste(SAMP_GROUP,"_",MODEL,".cs.png",sep = ""),sep = "/")
  
  colnames(train_outputs) <- MODEL
  colnames(test_outputs) <- MODEL
  colnames(val_outputs) <- MODEL
  
  MODEL_LIST %>% saveRDS(rds_output)
  train_outputs %>% write.table(train_stats,sep="\t",row.names=T, quote=F)
  test_outputs %>% write.table(test_stats,sep="\t",row.names=T, quote=F)
  val_outputs %>% write.table(val_stats,sep="\t",row.names=T, quote=F)
  ggsave(auc_figure, height = 5, width = 6, units = "in", dpi = 600, AUC_figure)
  ggsave(classifier_figure, width = 14, height = 14, units = "cm", dpi = 600, CS_figure)
  
  
  if (count == 1){total_train_outputs <- train_outputs; total_test_outputs <- test_outputs; total_val_outputs <- val_outputs}
  else{ total_train_outputs <- cbind(total_train_outputs,train_outputs)
  total_test_outputs <- cbind(total_test_outputs,test_outputs)
  total_val_outputs <- cbind(total_val_outputs,val_outputs) }
  count = count + 1
}
