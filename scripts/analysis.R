# - Analyzing count data ------------------------------------------------------
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(patchwork)
if (!dir.exists("results")) dir.create("results")
if (!dir.exists("results/plots")) dir.create("results/plots")

# - Read in data and info -----------------------------------------------------
# read in counts
count_files <- list.files(
  pattern = "ReadsPerGene.out.tab",
  full.names = TRUE,
  recursive = TRUE
)
samples <- str_extract(count_files, "Sample_[0-9]+")
counts <- lapply(count_files, function(x) {
  read.csv(x, skip = 4, sep = "\t", header = FALSE, row.names = 1) %>% select(V4)
})
counts <- bind_cols(counts) %>% tibble::rownames_to_column("Gene")
colnames(counts) <- c("Gene", samples)

# read in metadata
metadata <- bind_rows(
  read.csv("data/Run_2425/Run_2425_seeley.csv", skip = 18),
  read.csv("data/Run_2643/Run_2643_seeley.csv", skip = 18)
) %>%
  select(Sample_ID, Description) %>%
  filter(!duplicated(Sample_ID))

# format group info
metadata <- metadata %>%
  mutate("Pair" = str_extract(Description, "^[\\S]{7}")) %>%
  mutate("Cells" = rep(c("GFP", "Sup"), times = nrow(metadata)/2)) %>%
  mutate(Cells = factor(Cells, levels = c("Sup", "GFP")))

# add food status to metadata object
path <- "data/Run_2643/IlluminaSampleCore-FACS samples Bozadjieva "
food <- bind_rows(
  readxl::read_xlsx(paste0(path, "7-9-18.xlsx"), skip = 10),
  readxl::read_xlsx(paste0(path, "12-18-18.xlsx"), skip = 10)
) %>%
  select(Mouse, `Body weight`, `Sample Name*`) %>%
  mutate("Pair" = str_extract(`Sample Name*`, "^[\\S]{7}")) %>%
  filter(!is.na(Mouse)) %>%
  mutate("Food" = ifelse(str_detect(Mouse, "Chow"), "Chow", "HFD"))
metadata <- left_join(
  metadata,
  select(food, `Body weight`, Pair, Food),
  by = "Pair"
) %>%
  mutate(Food = factor(Food, levels = c("Chow", "HFD"))) %>%
  rename("BW" = `Body weight`)

# read in gene info
genes <- read.csv("data/genes.csv") %>%
  select(gene_id, gene_name) %>%
  filter(!duplicated(gene_id))

# - Sample QC -----------------------------------------------------------------
# sample read depth
p <- colSums(counts[, 2:ncol(counts)]) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("Sample_ID") %>%
  left_join(select(metadata, Sample_ID, Description, Cells),
            by = "Sample_ID") %>%
  mutate("Counts" = `.` / 10^6) %>%
  ggplot(aes(x = Sample_ID, y = Counts, fill = Cells)) +
  geom_col(color = "black") +
  theme_minimal() +
  scale_fill_manual(values = c("gray", "#73A657")) +
  scale_y_continuous(expand = c(0, 0)) +
  labs(y = "Mapped reads (M)", x = element_blank()) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
  geom_hline(aes(yintercept = 30), linetype = "dashed")
ggsave("results/plots/depth.png", p, width = 6, height = 4, units = "in",
       dpi = 400, bg = "white")
# -- Sample_123059 is pretty low --

# enrichment of EGFP:L10a by pair
cpm <- apply(counts[, 2:ncol(counts)], 2, function(x) {
  x / sum(x) * 10^6
})
cpm <- as.data.frame(cpm) %>% mutate("Gene" = counts$Gene)
p <- cpm %>% filter(Gene == "Gfp_L10a") %>% select(-Gene) %>%
  tidyr::pivot_longer(everything(), names_to = "Sample_ID", values_to = "CPM") %>%
  left_join(select(metadata, Sample_ID, Pair, Cells), by = "Sample_ID") %>%
  ggplot(aes(x = Cells, y = CPM, group = Pair, color = Pair)) +
  geom_point() +
  geom_line(show.legend = FALSE) +
  scale_color_brewer(type = "qual", palette = "Paired",
                     name = expression(underline("Sample"))) +
  scale_y_continuous(trans = "log2") +
  theme_minimal() +
  labs(y = "Counts per million", x = element_blank()) +
  theme(panel.grid = element_blank(),
        axis.text.x = element_text(color = "black", size = 10))
ggsave("results/plots/GFP.png", p, width = 3, height = 3, units = "in",
       dpi = 400, bg = "white")
# -- good enrichment of most samples, except 15LT462 --

# heatmap plot of expressed genes
cpm <- tibble::column_to_rownames(cpm, "Gene")
expressed <- apply(cpm, 1, function(x) sum(x >= 1) >= 4)
cpm <- cpm[expressed, ]
# make log2
log2cpm <- apply(cpm, c(1, 2), function(x) log2(x + 1))
# cluster
tree1 <- hclust(dist(log2cpm))
tree2 <- hclust(dist(t(log2cpm)))
df <- log2cpm %>% as.data.frame() %>%
  tibble::rownames_to_column("Gene") %>%
  tidyr::pivot_longer(-Gene, names_to = "Sample_ID", values_to = "CPM") %>%
  mutate("Gene" = factor(Gene, levels = tree1$labels[tree1$order])) %>%
  mutate("Sample_ID" = factor(Sample_ID, levels = tree2$labels[tree2$order]))
p1 <- ggdendro::ggdendrogram(as.dendrogram(tree2)) +
  theme_void() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = tree2$labels[tree2$order])
group_colors <- c(
  "Chow Sup" = "gray70",
  "HFD Sup" = "gray30",
  "Chow GFP" = "lightgreen",
  "HFD GFP" = "darkgreen"
)
p2 <- metadata %>% mutate("Group" = paste(Food, Cells)) %>%
  mutate(Group = factor(Group, levels = names(group_colors))) %>%
  mutate(Sample_ID = factor(Sample_ID, levels = tree2$labels[tree2$order])) %>%
  ggplot(aes(x = Sample_ID, y = "", fill = Group)) +
  geom_tile() +
  scale_fill_manual(values = group_colors) +
  theme_void() +
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.width = unit(0.05, "in"),
        legend.key.height = unit(0.05, "in"))
p3 <- ggplot(df, aes(x = Sample_ID, y = Gene, fill = CPM)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  theme_void() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 8, color = "black"),
        legend.key.width = unit(0.02, "in"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))
p <- p1 + p2 + p3 + plot_layout(heights = c(0.2, 0.05, 1), ncol = 1)
ggsave("results/plots/heatmap.png", p, width = 4, height = 6, units = "in",
       dpi = 400, bg = "white")

# - MDS plot -
mds <- log2cpm %>% t() %>% dist() %>% cmdscale()
p <- mds %>% as.data.frame() %>%
  tibble::rownames_to_column("Sample_ID") %>%
  left_join(metadata, by = "Sample_ID") %>%
  ggplot(aes(x = V1, y = V2, color = paste(Food, Cells))) +
  geom_point() +
  geom_text_repel(aes(label = Sample_ID), size = 2) +
  scale_color_manual(values = group_colors, name = "Group") +
  theme_classic() +
  xlab("Dim 1") +
  ylab("Dim 2") +
  theme(axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))
ggsave("results/plots/mds.png", p, width = 5, height = 4, units = "in",
       dpi = 400, bg = "white")

# - sample correlation -
corr <- cor(log2cpm)
tree <- hclust(dist(corr))
df <- corr %>% as.data.frame() %>%
  tibble::rownames_to_column("Var1") %>%
  tidyr::pivot_longer(-Var1, names_to = "Var2", values_to = "R") %>%
  mutate_at(vars(Var1, Var2), ~ factor(.x, levels = tree$labels[tree$order]))
p1 <- ggdendro::ggdendrogram(as.dendrogram(tree)) +
  theme_void() +
  scale_y_continuous(expand = c(0, 0)) +
  scale_x_discrete(limits = tree$labels[tree$order])
p2 <- metadata %>% mutate("Group" = paste(Food, Cells)) %>%
  mutate(Group = factor(Group, levels = names(group_colors))) %>%
  mutate(Sample_ID = factor(Sample_ID, levels = tree$labels[tree$order])) %>%
  ggplot(aes(x = Sample_ID, y = "", fill = Group)) +
  geom_tile() +
  scale_fill_manual(values = group_colors) +
  theme_void() +
  theme(legend.text = element_text(size = 6),
        legend.title = element_text(size = 7),
        legend.key.width = unit(0.05, "in"),
        legend.key.height = unit(0.05, "in"))
p3 <- ggplot(df, aes(x = Var1, y = Var2, fill = R)) +
  geom_tile() +
  scale_fill_viridis_c(option = "magma") +
  theme_void() +
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5,
                                   size = 8),
        axis.text.y = element_text(size = 8, hjust = 1),
        legend.key.width = unit(0.02, "in"),
        legend.text = element_text(size = 6),
        legend.title = element_text(size = 7))
p <- p1 + p2 + p3 + plot_layout(heights = c(0.2, 0.05, 1), ncol = 1)
ggsave("results/plots/correlation.png", p, width = 4, height = 4.4, units = "in",
       dpi = 400, bg = "white")

# - Regulation ----------------------------------------------------------------
outliers <- c("Sample_123059", "Sample_113435")
regulation_info <- filter(metadata, !Sample_ID %in% outliers & Cells == "GFP")

# run DESeq2 to identify genes regulated by HFD
mtx <- select(counts, Gene, regulation_info$Sample_ID) %>%
  tibble::column_to_rownames("Gene")
dds <- DESeq2::DESeqDataSetFromMatrix(
  countData = mtx,
  colData = regulation_info,
  design = ~ Food
) %>%
  DESeq2::DESeq()

regulation <- DESeq2::results(dds, contrast = c("Food", "HFD", "Chow")) %>%
  as.data.frame() %>%
  tibble::rownames_to_column("gene_id") %>%
  arrange(desc(padj < 0.05), desc(abs(log2FoldChange)))

cpm <- DESeq2::fpm(dds)
cpm_means <- sapply(unique(regulation_info$Food), function(x) {
  rowMeans(cpm[, regulation_info$Food == x])
})
colnames(cpm_means) <- unique(regulation_info$Food)

# - MA plot and volcano plot -
ma <- data.frame(
  "M" = log2(cpm_means[, "HFD"]) - log2(cpm_means[, "Chow"]),
  "A" = rowMeans(log2(cpm_means)),
  "P" = regulation$padj[match(rownames(cpm), regulation$gene_id)],
  "Gene" = rownames(cpm_means)
) %>%
  filter(!is.na(P)) %>%
  left_join(genes, by = c("Gene" = "gene_id")) %>%
  arrange(P < 0.05)
p <- ggplot(ma, aes(x = A, y = M, color = P < 0.05)) +
  geom_point(stroke = 0) +
  geom_hline(aes(yintercept = 0)) +
  geom_text_repel(data = tail(ma, 20), aes(label = gene_name), size = 2) +
  scale_color_manual(values = c("gray80", "gray20")) +
  scale_y_continuous(expand = expansion(mult = c(0, 0.05))) +
  xlab("Expression") +
  ylab("Regulation") +
  theme_classic() +
  theme(axis.text = element_text(size = 7, color = "black"),
        axis.title = element_text(size = 8),
        legend.position = "none")
ggsave("results/plots/ma_plot.png", p, width = 4, height = 3, units = "in",
       dpi = 400, bg= "white")

# - Pathway analysis ----------------------------------------------------------
convert <- function(genes, from = "SYMBOL", to = "ENTREZID") {
  clusterProfiler::bitr(
    genes,
    fromType = from,
    toType = to,
    OrgDb = org.Mm.eg.db::org.Mm.eg.db
  )
}
get_enrichment <- function(direction) {
  sign <- switch(direction, "up" = 1, "down" = -1)
  df <- left_join(regulation, genes, by = "gene_id")
  de <- filter(df, padj < 0.05 & sign(log2FoldChange) == sign)$gene_name
  universe <- filter(df, !is.na(padj))$gene_name
  clusterProfiler::enrichKEGG(
    gene = convert(de) %>% pull(ENTREZID),
    organism = "mmu",
    keyType = "ncbi-geneid",
    universe = convert(universe) %>% pull(ENTREZID)
  ) %>%
    as.data.frame() %>%
    mutate("Direction" = direction)
}
kegg <- lapply(c("up", "down"), get_enrichment)
kegg <- bind_rows(kegg)
write.csv(kegg, "results/kegg.csv", row.names = FALSE, na = "")

plot_kegg <- function(direction = NULL) {
  if (!is.null(direction)) kegg <- filter(kegg, Direction == direction)
  kegg$Description = factor(kegg$Description, levels = kegg$Description[nrow(kegg):1])
  ggplot(kegg, aes(x = -log10(qvalue), y = Description, fill = Count)) +
    geom_col(color = "black") +
    scale_fill_gradient(low = "gray90", high = "firebrick4",
                        limits = c(0, max(kegg$Count))) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.05)), limits = c(0, NA)) +
    theme_classic() +
    xlab(expression(italic("P")*" value (-log"[10]*")")) +
    ylab(NULL) +
    theme(axis.text.x = element_text(size = 7, color = "black"),
          axis.text.y = element_text(size = 8, color = "black", hjust = 1),
          axis.title = element_text(size = 8),
          legend.text = element_text(size = 6),
          legend.title = element_text(size = 7),
          legend.key.width = unit(0.02, "in"))
}
p <- plot_kegg()
ggsave("results/plots/kegg.png", p, width = 5, height = 3, units = "in",
       dpi = 400, bg = "white")

# - Write data ----------------------------------------------------------------
regulation %>% left_join(genes, by = "gene_id") %>%
  full_join(tibble::rownames_to_column(as.data.frame(cpm_means), "gene_id")) %>%
  select(gene_name, gene_id, Chow, HFD, log2FoldChange, padj) %>%
  mutate_at(vars(Chow, HFD, log2FoldChange), ~round(.x, 2)) %>%
  mutate(padj = format(padj, scientific = TRUE, digits = 2)) %>%
  rename("FC" = log2FoldChange, "P" = padj) %>%
  write.csv("results/results.csv", row.names = FALSE, na = "")
