install.packages("dplyr")
install.packages("RColorBrewer")
install.packages("ggrepel")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("TxDb.Mmusculus.UCSC.mm39.knownGene")
BiocManager::install("ChIPseeker")
BiocManager::install("clusterProfiler")
BiocManager::install("org.Mm.eg.db") # mouse organism annotation
BiocManager::install("biomaRt")
BiocManager::install("rtracklayer")


library(dplyr)
library(TxDb.Mmusculus.UCSC.mm39.knownGene)
library(ChIPseeker)
library(clusterProfiler)
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggrepel)
library(biomaRt)
library(rtracklayer)

projPath = "~/Charite Thesis/H3K9ac Cut&Tag/Results/"
setwd(projPath)


####################### Differential enrichment analysis #############################

# diff_files <-
#   list.files("differential_analysis", pattern = "all_prob_sep_2k_without_outliers_padj0.1.ucsc.csv", full.names = TRUE)
# diff_files <- as.list(diff_files)
# diff_files <- lapply(diff_files, read.table, sep = ',', header = TRUE)

diff_KD <- read.table(
  "differential_analysis/KD_masterPeak_prob_sep_2k_without_outliers_padj0.1.ucsc.csv",
  sep = ',',
  header = TRUE
)
enriched_KD <- diff_KD[diff_KD$diffexpressed == 'UP', ]
lost_KD <- diff_KD[diff_KD$diffexpressed == 'DOWN', ]
diff_OE <- read.table(
  "differential_analysis/OE_masterPeak_prob_sep_2k_without_outliers_padj0.1.ucsc.csv",
  sep = ',',
  header = TRUE
)
enriched_OE <- diff_OE[diff_OE$diffexpressed == 'UP', ]
lost_OE <- diff_OE[diff_OE$diffexpressed == 'DOWN', ]

diff_files <- list(enriched_KD, lost_KD, enriched_OE, lost_OE)

names(diff_files) <-
  c("Enriched_KD", "Lost_KD", "Enriched_OE", "Lost_OE")

diff_files <- lapply(diff_files, makeGRangesFromDataFrame, keep.extra.columns =
                       TRUE)

# Peak annotation

txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
peakAnnoList <-
  lapply(diff_files,
         annotatePeak,
         TxDb = txdb,
         tssRegion = c(-1000, 1000))

# Create a list with genes from each sample

genes <- lapply(peakAnnoList, function(i) {
  as.data.frame(i)$geneId
})

# Enrichement

compGO <- compareCluster(
  geneCluster = genes,
  fun = "enrichGO",
  #organism = "mouse",
  pvalueCutoff = 0.1,
  OrgDb = org.Mm.eg.db,
  pAdjustMethod = "BH"
)

dotplot(compGO,
        showCategory = 7,
        font.size = 8,
        title = "GO Enrichment Analysis")

cluster_result_df <- compGO@compareClusterResult

# Extract genes ID of the desired/dominant enrichment class

lost_KD_df <-
  subset(cluster_result_df, Cluster == "Lost_KD")
# extract class
genes_in_major_class = lost_KD_df[which.min(lost_KD_df$p.adjust), ]$geneID
genes_in_major_class
# genes_in_major_class <-
#   lost_KD_df$geneID
# genes_in_major_class
genes_in_major_class <- strsplit(genes_in_major_class, split = "/")


# Get annotation dataframe for each cluster

enriched_KD_anno <-
  data.frame(peakAnnoList[["Enriched_KD"]]@anno)
enriched_OE_anno <- data.frame(peakAnnoList[["Enriched_OE"]]@anno)
lost_KD_anno <-
  data.frame(peakAnnoList[["Lost_KD"]]@anno)
lost_OE_anno <-
  data.frame(peakAnnoList[["Lost_OE"]]@anno)


# Get gene names from gene ids

mart <-
  useEnsembl(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")
tmp_enriched_KD <- data.frame("entrezgene_id" = enriched_KD_anno$geneId)
tmp_enriched_OE <- data.frame("entrezgene_id" = enriched_OE_anno$geneId)
tmp_lost_KD <- data.frame("entrezgene_id" = lost_KD_anno$geneId)
tmp_lost_OE <- data.frame("entrezgene_id" = lost_OE_anno$geneId)

ids_enriched_KD <- getBM(
  attributes = c("entrezgene_id", "external_gene_name"),
  filters = "entrezgene_id",
  values = tmp_enriched_KD[, "entrezgene_id"],
  mart = mart
)
ids_enriched_OE <- getBM(
  attributes = c("entrezgene_id", "external_gene_name"),
  filters = "entrezgene_id",
  values = tmp_enriched_OE[, "entrezgene_id"],
  mart = mart
)
ids_lost_KD <- getBM(
  attributes = c("entrezgene_id", "external_gene_name"),
  filters = "entrezgene_id",
  values = tmp_lost_KD[, "entrezgene_id"],
  mart = mart
)
ids_lost_OE <- getBM(
  attributes = c("entrezgene_id", "external_gene_name"),
  filters = "entrezgene_id",
  values = tmp_lost_OE[, "entrezgene_id"],
  mart = mart
)

# Merge Gene Annotations

names(ids_enriched_KD)[names(ids_enriched_KD) == "entrezgene_id"] <- "geneId"
ids_enriched_KD$geneId <- as.character(ids_enriched_KD$geneId)
enriched_KD_anno <- left_join(
  x = enriched_KD_anno,
  y = ids_enriched_KD,
  by = "geneId",
  relationship = "many-to-many"
)
enriched_KD_anno <- enriched_KD_anno[, -c(7:8)]
#enriched_KD_anno_major_class <-
#subset(enriched_KD_anno, geneId %in% genes_in_major_class[[1]])

names(ids_enriched_OE)[names(ids_enriched_OE) == "entrezgene_id"] <- "geneId"
ids_enriched_OE$geneId <- as.character(ids_enriched_OE$geneId)
enriched_OE_anno <- left_join(
  x = enriched_OE_anno,
  y = ids_enriched_OE,
  by = "geneId",
  relationship = "many-to-many"
)
enriched_OE_anno <- enriched_OE_anno[, -c(7:8)]

names(ids_lost_KD)[names(ids_lost_KD) == "entrezgene_id"] <- "geneId"
ids_lost_KD$geneId <- as.character(ids_lost_KD$geneId)
lost_KD_anno <- left_join(
  x = lost_KD_anno,
  y = ids_lost_KD,
  by = "geneId",
  relationship = "many-to-many"
)
lost_KD_anno <- lost_KD_anno[, -c(7:8)]

names(ids_lost_OE)[names(ids_lost_OE) == "entrezgene_id"] <- "geneId"
ids_lost_OE$geneId <- as.character(ids_lost_OE$geneId)
lost_OE_anno <- left_join(
  x = lost_OE_anno,
  y = ids_lost_OE,
  by = "geneId",
  relationship = "many-to-many"
)
lost_OE_anno <- lost_OE_anno[, -c(7:8)]

# Volcano plot

diff_KD_anno <- rbind(enriched_KD_anno, lost_KD_anno)
diff_KD_anno_gr <- GRanges(
  seqnames = diff_KD_anno$seqnames,
  ranges = IRanges(start = diff_KD_anno$start, end = diff_KD_anno$end)
)
export(
  diff_KD_anno_gr,
  "differential_analysis/diff_KD_masterPeak_sep_2k_without_outliers.bed",
  "bed"
)
diff_KD_anno <- diff_KD_anno[c("PeakID", "external_gene_name")]

diff_OE_anno <- rbind(enriched_OE_anno, lost_OE_anno)
diff_OE_anno_gr <- GRanges(
  seqnames = diff_OE_anno$seqnames,
  ranges = IRanges(start = diff_OE_anno$start, end = diff_OE_anno$end)
)
export(
  diff_OE_anno_gr,
  "differential_analysis/diff_OE_masterPeak_sep_2k_without_outliers.bed",
  "bed"
)
diff_OE_anno <- diff_OE_anno[c("PeakID", "external_gene_name")]


diff_KD <- left_join(
  x = diff_KD,
  y = diff_KD_anno,
  by = "PeakID",
  relationship = "many-to-many"
)

diff_OE <- left_join(
  x = diff_OE,
  y = diff_OE_anno,
  by = "PeakID",
  relationship = "many-to-many"
)


p <- ggplot(
  data = diff_KD,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    col = diffexpressed,
    label = external_gene_name
  )
) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red") +
  geom_hline(yintercept = -log10(0.1), col = "red")
p

p1 <- ggplot(
  data = diff_OE,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    col = diffexpressed,
    label = external_gene_name
  )
) +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red") +
  geom_hline(yintercept = -log10(0.1), col = "red")
p1

#make change color per chr

######################## General enrichment analysis ###############################

samplefiles <-
  list.files("SEACR/UCSC", pattern = ".bed", full.names = TRUE)
samplefiles <- as.list(samplefiles)
names(samplefiles) <-
  c("Control_rep1",
    "Control_rep2",
    "TAU_KD_rep1",
    "TAU_KD_rep2",
    "TAU_OE_rep1")


# Peak annotation

txdb <- TxDb.Mmusculus.UCSC.mm39.knownGene
peakAnnoList <-
  lapply(samplefiles,
         annotatePeak,
         TxDb = txdb,
         tssRegion = c(-1000, 1000))

# Create a list with genes from each sample

genes <- lapply(peakAnnoList, function(i) {
  as.data.frame(i)$geneId
})

compKEGG <- compareCluster(
  geneCluster = genes,
  fun = "enrichKEGG",
  organism = "mouse",
  pvalueCutoff = 0.1,
  #OrgDb = org.Mm.eg.db,
  pAdjustMethod = "BH"
)

dotplot(
  compKEGG,
  showCategory = 7,
  font.size = 8,
  title = "KEGG Enrichment Analysis"
)


