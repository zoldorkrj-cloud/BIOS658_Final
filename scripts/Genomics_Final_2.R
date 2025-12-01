# Set up and data import

# Core packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(edgeR)
library(pheatmap)
library(tibble)

# If needed for enrichment
# install.packages(c("fgsea"))
# BiocManager::install(c("org.Hs.eg.db","clusterProfiler","AnnotationDbi"))
library(org.Hs.eg.db)
library(clusterProfiler)
library(fgsea)

# 1) Load GEO meta
library(GEOquery)
GSE <- "GSE72536"
gds <- getGEO(GSE, GSEMatrix = FALSE, AnnotGPL = TRUE)

# Parse meta (HPV status from characteristics)
titles <- geo_accessions <- hpv_statuses <- c()
for (gsm in names(gds@gsms)) {
  titles <- c(titles, gds@gsms[[gsm]]@header$title)
  geo_accessions <- c(geo_accessions, gds@gsms[[gsm]]@header$geo_accession)
  
  chars <- gds@gsms[[gsm]]@header$characteristics_ch1
  txt <- if (is.null(chars)) "" else paste(chars, collapse = " ; ")
  hpv <- NA_character_
  if (grepl("hpv", txt, ignore.case = TRUE)) {
    if (grepl("positive|pos|detected", txt, ignore.case = TRUE)) hpv <- "HPVpos"
    if (grepl("negative|neg|not detected|undetected", txt, ignore.case = TRUE)) hpv <- "HPVneg"
    if (is.na(hpv)) hpv <- stringr::str_trim(txt)
  }
  hpv_statuses <- c(hpv_statuses, hpv)
}
meta <- data.frame(sample = geo_accessions, title = titles, hpv_status = hpv_statuses,
                   stringsAsFactors = FALSE)

# 2) Load counts from GEO
library(readr)
mtx <- read_tsv(
  "https://www.ncbi.nlm.nih.gov/geo/download/?type=rnaseq_counts&acc=GSE72536&format=file&file=GSE72536_raw_counts_GRCh38.p13_NCBI.tsv.gz",
  show_col_types = FALSE
)

# Sanity: sample names match meta
stopifnot(all.equal(colnames(mtx)[2:ncol(mtx)], meta$sample) == TRUE)

# 3) Clean gene IDs and aggregate duplicates
mtx_clean <- mtx %>%
  filter(!is.na(GeneID) & GeneID != "") %>%
  group_by(GeneID) %>%
  summarise(across(where(is.numeric), ~ sum(.x, na.rm = TRUE)), .groups = "drop")

# 4) Move GeneID to rownames (crucial to avoid treating it as a sample)
mtx_mat <- mtx_clean %>%
  as.data.frame() %>%
  column_to_rownames("GeneID")
stopifnot(!"GeneID" %in% colnames(mtx_mat))  # must be TRUE

out_dir <- "C:/Users/ryanz/OneDrive/Desktop/Genomics_Final/data"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# QA
# Library sizes (sum of counts per sample)
lib_sizes <- colSums(mtx_mat, na.rm = TRUE)
lib_df <- data.frame(sample = names(lib_sizes), lib_size = lib_sizes)

p_lib <- ggplot(lib_df, aes(x = reorder(sample, lib_size), y = lib_size)) +
  geom_col(fill = "steelblue") +
  coord_flip() +
  scale_y_continuous(labels = scales::comma) +
  theme_minimal() +
  theme(axis.text.y = element_text(size = 8)) +
  labs(title = "Library sizes", x = "Sample", y = "Total counts")
ggsave(file.path(out_dir, "qc_library_sizes.png"), p_lib, width = 8, height = 6, dpi = 300)

# Log2-TPM-like visualization (log2 of counts+1 is fine for QC)
log2_tpm <- as.data.frame(mtx_mat)
log2_tpm[] <- log2(log2_tpm + 1)

# Density per sample
counts_long <- log2_tpm %>%
  rownames_to_column("GeneID") %>%
  pivot_longer(cols = -GeneID, names_to = "sample", values_to = "log2TPM")
p_dens <- ggplot(counts_long, aes(x = log2TPM, color = sample)) +
  geom_density(linewidth = 0.5, alpha = 0.7) +
  theme_minimal() +
  theme(legend.position = "right") +
  labs(title = "Density of log2 counts+1 per sample", x = "log2(counts+1)", y = "Density")
ggsave(file.path(out_dir, "qc_density.png"), p_dens, width = 8, height = 6, dpi = 300)

# Boxplot ordered by median
sample_order <- counts_long %>%
  group_by(sample) %>%
  summarise(med = median(log2TPM, na.rm = TRUE)) %>%
  arrange(med) %>% pull(sample)
p_box <- ggplot(counts_long, aes(x = factor(sample, levels = sample_order), y = log2TPM)) +
  geom_boxplot(outlier.size = 0.4, fill = "lightgray") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8)) +
  labs(title = "Per-sample log2 counts+1 distributions", x = "Sample", y = "log2(counts+1)")
ggsave(file.path(out_dir, "qc_boxplot.png"), p_box, width = 10, height = 6, dpi = 300)

# Sample-sample correlation heatmap
log2_mat <- as.matrix(log2_tpm)
cols <- colorRampPalette(c("navy","white","firebrick3"))(50)
pheatmap(cor(log2_mat, method = "pearson"), color = cols,
         main = "Sample correlation (Pearson)",
         show_rownames = TRUE, show_colnames = TRUE, fontsize = 8, border_color = NA,
         filename = file.path(out_dir, "qc_correlation_heatmap.png"))

# PCA on top 20% variable genes (log2 counts+1)
vars <- apply(log2_mat, 1, var)
n_top <- ceiling(0.2 * length(vars))
top_genes <- names(sort(vars, decreasing = TRUE))[1:n_top]
pca <- prcomp(t(log2_mat[top_genes, ]), center = TRUE, scale. = TRUE)
pct <- round(100 * summary(pca)$importance[2, 1:2], 1)
scores <- data.frame(sample = rownames(pca$x), PC1 = pca$x[,1], PC2 = pca$x[,2])
p_pca <- ggplot(scores, aes(PC1, PC2, label = sample)) +
  geom_point(size = 2) +
  geom_text(hjust = 0, nudge_x = 0.02, size = 3) +
  labs(title = "PCA (top 20% variable genes, log2 counts+1)",
       x = paste0("PC1 (", pct[1], "%)"), y = paste0("PC2 (", pct[2], "%)")) +
  theme_minimal()
ggsave(file.path(out_dir, "qc_pca_before.png"), p_pca, width = 8, height = 6, dpi = 300)

# Preprocessing (filter and TMM)

# Prepare DGEList with clean matrix (genes as rows, samples as columns)
y <- DGEList(counts = mtx_mat)

# Filter lowly expressed genes (no design known yet → global filter)
keep <- filterByExpr(y, design = NULL)
y <- y[keep, , keep.lib.sizes = FALSE]

# TMM normalization
y <- calcNormFactors(y)

# Log2-CPM for visualization/EDA
log2cpm <- cpm(y, log = TRUE, prior.count = 1)  # genes x samples

# Batch effect correction
# Merge meta to ensure order alignment
stopifnot(identical(colnames(y), meta$sample))

# If batch is unknown: SVA to estimate surrogate variables
library(sva)
# Model including biological variable if known (HPV)
meta$HPV <- factor(meta$hpv_status)
mod <- model.matrix(~ HPV, data = meta)         # biological model
mod0 <- model.matrix(~ 1, data = meta)          # null model
svobj <- sva(log2cpm, mod, mod0)                # estimates sv (svobj$n.sv)

# Regress out surrogate variables: build augmented design
mod_sv <- cbind(mod, svobj$sv)                  # include SVs in design for downstream DE
# For visualization of "batch-corrected" expression, you can remove SV effects via limma:
library(limma)
fit <- lmFit(log2cpm, design = mod_sv)
beta_sv <- fit$coefficients[, (ncol(mod)+1):ncol(mod_sv), drop = FALSE]
# beta_sv: genes x SVs
# svobj$sv: samples x SVs
# transpose svobj$sv to SVs x samples
sv_contrib <- beta_sv %*% t(svobj$sv)   # genes x samples
# contribution of SVs
log2cpm_corrected <- log2cpm - sv_contrib


#PCA before and after correction
# Before
vars_b <- apply(log2cpm, 1, var)
top_b <- names(sort(vars_b, decreasing = TRUE))[1:ceiling(0.1*length(vars_b))]
pca_b <- prcomp(t(log2cpm[top_b, ]), center = TRUE, scale. = TRUE)
pct_b <- round(100 * summary(pca_b)$importance[2, 1:2], 1)
scores_b <- data.frame(sample = rownames(pca_b$x), PC1 = pca_b$x[,1], PC2 = pca_b$x[,2])
p_pca_b <- ggplot(scores_b, aes(PC1, PC2, color = meta$HPV, label = sample)) +
  geom_point(size = 2) + geom_text(hjust = 0, nudge_x = 0.02, size = 3) +
  labs(title = "PCA before batch correction",
       x = paste0("PC1 (", pct_b[1], "%)"), y = paste0("PC2 (", pct_b[2], "%)")) +
  theme_minimal()
ggsave(file.path(out_dir, "pca_before_batch.png"), p_pca_b, width = 8, height = 6, dpi = 300)

# After (SVA-removed)
vars_a <- apply(log2cpm_corrected, 1, var)
top_a <- names(sort(vars_a, decreasing = TRUE))[1:ceiling(0.1*length(vars_a))]
pca_a <- prcomp(t(log2cpm_corrected[top_a, ]), center = TRUE, scale. = TRUE)
pct_a <- round(100 * summary(pca_a)$importance[2, 1:2], 1)
scores_a <- data.frame(sample = rownames(pca_a$x), PC1 = pca_a$x[,1], PC2 = pca_a$x[,2])
p_pca_a <- ggplot(scores_a, aes(PC1, PC2, color = meta$HPV, label = sample)) +
  geom_point(size = 2) + geom_text(hjust = 0, nudge_x = 0.02, size = 3) +
  labs(title = "PCA after SVA correction",
       x = paste0("PC1 (", pct_a[1], "%)"), y = paste0("PC2 (", pct_a[2], "%)")) +
  theme_minimal()
ggsave(file.path(out_dir, "pca_after_batch.png"), p_pca_a, width = 8, height = 6, dpi = 300)

# Heircarchial clustering top 20% before and after
# Before
top10_b <- names(sort(vars_b, decreasing = TRUE))[1:ceiling(0.1*length(vars_b))]
pheatmap(log2cpm[top10_b, ], show_rownames = FALSE, show_colnames = TRUE,
         scale = "row", clustering_method = "complete",
         main = "HC before correction (top 10% variable genes)",
         filename = file.path(out_dir, "hc_before_correction.png"))

# After
top10_a <- names(sort(vars_a, decreasing = TRUE))[1:ceiling(0.1*length(vars_a))]
pheatmap(log2cpm_corrected[top10_a, ], show_rownames = FALSE, show_colnames = TRUE,
         scale = "row", clustering_method = "complete",
         main = "HC after SVA correction (top 10% variable genes)",
         filename = file.path(out_dir, "hc_after_correction.png"))

# EdgeR Differential Expression Analysis

# Build design with HPV (and SVs if using SVA)
group <- factor(meta$HPV)  # levels like HPVneg/HPVpos
design <- model.matrix(~ group + svobj$sv, data = meta)  # add SVs as covariates

# DGEList from original counts (filtered to expressed genes again)
y <- DGEList(counts = mtx_mat)
keep <- filterByExpr(y, design = design)
y <- y[keep, , keep.lib.sizes = FALSE]
y <- calcNormFactors(y)

# Estimate dispersion and fit GLM
y <- estimateDisp(y, design)
fit <- glmFit(y, design)
lrt <- glmLRT(fit, coef = "groupHPVpos")  # adjust depending on reference level

# Results table
de_tab <- topTags(lrt, n = Inf)$table
de_tab$FDR <- p.adjust(de_tab$PValue, method = "BH")

# Save results
write.csv(de_tab, file.path(out_dir, "edgeR_DE_HPVs_vs_HPVn.csv"), row.names = TRUE)

# Map Entrez IDs to SYMBOLs
gene_symbols <- mapIds(org.Hs.eg.db,
                       keys = rownames(de_tab),
                       column = "SYMBOL",
                       keytype = "ENTREZID",
                       multiVals = "first")

de_tab$SYMBOL <- gene_symbols
de_tab$ENTREZID <- rownames(de_tab)   # keep Entrez IDs

# Top 10 DE genes by FDR
top10 <- rownames(de_tab[order(de_tab$FDR), ])[1:10]

# Expression values for those genes
expr_logcpm <- cpm(y, log = TRUE, prior.count = 1)
expr_long <- as.data.frame(expr_logcpm[top10, ]) %>%
  rownames_to_column("ENTREZID") %>%
  mutate(SYMBOL = gene_symbols[ENTREZID]) %>%
  pivot_longer(cols = -c(ENTREZID, SYMBOL),
               names_to = "sample", values_to = "logCPM") %>%
  left_join(meta, by = c("sample" = "sample"))

# Boxplot with both SYMBOL and Entrez ID in facet labels
expr_long$GeneLabel <- paste0(expr_long$SYMBOL, " (", expr_long$ENTREZID, ")")

p_box_de <- ggplot(expr_long, aes(x = hpv_status, y = logCPM)) +
  geom_boxplot(outlier.size = 0.4, fill = "lightgray") +
  facet_wrap(~ GeneLabel, scales = "free_y", ncol = 5) +
  theme_minimal() +
  labs(title = "Top 10 DE genes: logCPM by HPV status",
       x = "HPV status", y = "logCPM")

ggsave(file.path(out_dir, "edgeR_top10_boxplots_with_entrez.png"),
       p_box_de, width = 12, height = 8, dpi = 300)

# Prepare volcano data
library(ggplot2)
library(ggrepel)

# Add Entrez IDs and SYMBOLs to DE table
de_tab$ENTREZID <- rownames(de_tab)
de_tab$SYMBOL <- mapIds(org.Hs.eg.db,
                        keys = rownames(de_tab),
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")

# Define up/down regulation relative to HPVpos vs HPVneg
volcano_df <- de_tab %>%
  mutate(direction = case_when(
    logFC > 0 & FDR < 0.05 ~ "Upregulated",
    logFC < 0 & FDR < 0.05 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Select top 5 up and down genes by FDR
top_up <- volcano_df %>%
  filter(direction == "Upregulated") %>%
  arrange(FDR) %>%
  head(5)

top_down <- volcano_df %>%
  filter(direction == "Downregulated") %>%
  arrange(FDR) %>%
  head(5)

top_labels <- bind_rows(top_up, top_down)

# Volcano plot
library(ggplot2)
library(ggrepel)

# Add SYMBOLs to DE table
de_tab$ENTREZID <- rownames(de_tab)
de_tab$SYMBOL <- mapIds(org.Hs.eg.db,
                        keys = rownames(de_tab),
                        column = "SYMBOL",
                        keytype = "ENTREZID",
                        multiVals = "first")

# Define up/down regulation relative to HPVpos vs HPVneg
volcano_df <- de_tab %>%
  mutate(direction = case_when(
    logFC > 0 & FDR < 0.05 ~ "Upregulated",
    logFC < 0 & FDR < 0.05 ~ "Downregulated",
    TRUE ~ "Not Significant"
  ))

# Select top 5 up and down genes by FDR
top_up <- volcano_df %>%
  filter(direction == "Upregulated") %>%
  arrange(FDR) %>%
  head(5)

top_down <- volcano_df %>%
  filter(direction == "Downregulated") %>%
  arrange(FDR) %>%
  head(5)

top_labels <- bind_rows(top_up, top_down)

# Volcano plot
p_volcano <- ggplot(volcano_df, aes(x = logFC, y = -log10(FDR))) +
  geom_point(aes(color = direction), alpha = 0.6) +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not Significant" = "gray")) +
  geom_text_repel(data = top_labels,
                  aes(label = SYMBOL),
                  size = 3, max.overlaps = Inf) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black") +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black") +
  theme_minimal() +
  labs(title = "Volcano plot (HPVpos vs HPVneg)",
       x = "log2 Fold Change (HPVpos vs HPVneg)",
       y = "-log10 FDR")

ggsave(file.path(out_dir, "edgeR_volcano_symbol.png"),
       p_volcano, width = 8, height = 6, dpi = 300)


# MA plot: logFC vs average expression
plotMD(lrt,
       column = "groupHPVpos",   # coefficient tested
       status = decideTests(lrt),# highlights significant genes
       main = "MA plot (HPVpos vs HPVneg)",
       xlab = "Average logCPM",
       ylab = "logFC")

# Save MA plot
png(file.path(out_dir, "edgeR_MAplot.png"), width = 800, height = 600)
plotMD(lrt, column = "groupHPVpos", status = decideTests(lrt),
       main = "MA plot (HPVpos vs HPVneg)",
       xlab = "Average logCPM", ylab = "logFC")
dev.off()

# Heatmap of top DE genes
library(pheatmap)

# Select top 50 DE genes by FDR
top50 <- rownames(de_tab[order(de_tab$FDR), ])[1:50]

# Extract normalized logCPM values for those genes
expr_top <- cpm(y, log = TRUE, prior.count = 1)[top50, ]

# Annotation for samples
annotation_col <- data.frame(HPV = meta$hpv_status)
rownames(annotation_col) <- meta$sample

# Save only the heatmap
png(file.path(out_dir, "edgeR_heatmap_top50.png"), width = 1000, height = 800)
pheatmap(expr_top,
         scale = "row",
         annotation_col = annotation_col,
         show_rownames = FALSE,
         clustering_method = "complete",
         main = "Top 50 DE genes heatmap")
dev.off()



# Functional enrichment analysis

# Prepare ranked gene list for GSEA (Entrez IDs for KEGG; SYMBOL for GO via clusterProfiler)
# Map SYMBOLs
gene_symbols <- mapIds(org.Hs.eg.db, keys = rownames(de_tab), keytype = "ENTREZID",
                       column = "SYMBOL", multiVals = "first")
de_tab$SYMBOL <- gene_symbols

# Ranked vector (named by SYMBOL)
ranked <- de_tab$logFC
names(ranked) <- de_tab$SYMBOL
ranked <- ranked[!is.na(names(ranked))]
ranked <- sort(ranked, decreasing = TRUE)

# fgsea with GO sets via clusterProfiler (convert to gene sets)
# Build GO gene sets
ego_bp <- enrichGO(gene = de_tab$SYMBOL[de_tab$FDR < 0.05],
                   OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                   ont = "BP", pAdjustMethod = "BH", qvalueCutoff = 0.25)
ego_mf <- enrichGO(gene = de_tab$SYMBOL[de_tab$FDR < 0.05],
                   OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                   ont = "MF", pAdjustMethod = "BH", qvalueCutoff = 0.25)
ego_cc <- enrichGO(gene = de_tab$SYMBOL[de_tab$FDR < 0.05],
                   OrgDb = org.Hs.eg.db, keyType = "SYMBOL",
                   ont = "CC", pAdjustMethod = "BH", qvalueCutoff = 0.25)

# ORA results plots
p_bp <- dotplot(ego_bp, showCategory = 20) + ggtitle("GO Biological Process ORA")
p_mf <- dotplot(ego_mf, showCategory = 20) + ggtitle("GO Molecular Function ORA")
p_cc <- dotplot(ego_cc, showCategory = 20) + ggtitle("GO Cellular Component ORA")
ggsave(file.path(out_dir, "ORA_GO_BP.png"), p_bp, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "ORA_GO_MF.png"), p_mf, width = 8, height = 6, dpi = 300)
ggsave(file.path(out_dir, "ORA_GO_CC.png"), p_cc, width = 8, height = 6, dpi = 300)

# KEGG ORA (requires Entrez IDs)
entrez_ids <- rownames(de_tab)
ekegg <- enrichKEGG(gene = entrez_ids[de_tab$FDR < 0.05],
                    organism = "hsa", pAdjustMethod = "BH", qvalueCutoff = 0.25)
p_kegg <- dotplot(ekegg, showCategory = 20) + ggtitle("KEGG ORA")
ggsave(file.path(out_dir, "ORA_KEGG.png"), p_kegg, width = 8, height = 6, dpi = 300)

#GSEA w KEGG

# Get KEGG pathways as gene sets via clusterProfiler helper
kegg_sets <- as.list(ekegg@result$geneID)

#KEGG PATHWAY "hsa05200" = Pathways in cancer

# Install if needed
# BiocManager::install("pathview")

library(pathview)

# Prepare a named vector of logFC values (Entrez IDs as names)
gene_fc <- de_tab$logFC
names(gene_fc) <- de_tab$ENTREZID   # must be Entrez IDs

# Pick a KEGG pathway ID (e.g., "hsa04612" = Antigen processing and presentation)
#pathway_id <- "hsa05200"  # Pathways in cancer
pathway_id <- "hsa04612"  # Antigen processing and presentation

# Generate KEGG pathway map with DE genes colored by logFC
pathview(gene.data  = gene_fc,
         pathway.id = pathway_id,
         species    = "hsa",       # human
         out.suffix = "HPV_DE",
         kegg.native = TRUE)
