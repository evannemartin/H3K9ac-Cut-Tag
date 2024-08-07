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
getOption("clusterProfiler.download.method")
library(org.Mm.eg.db)
library(RColorBrewer)
library(ggrepel)
library(biomaRt)
library(rtracklayer)
library(R.utils)

R.utils::setOption("clusterProfiler.download.method", "auto")

projPath = "~/Charite Thesis/H3K9ac Cut&Tag/Results/"
setwd(projPath)


####################### Differential enrichment analysis #############################

## Uncomment code to get dataframe of all annotations regardless of differential peaks ##


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

# diff_files <- list(diff_KD, diff_OE)

names(diff_files) <-
  c("Enriched_KD", "Lost_KD", "Enriched_OE", "Lost_OE")

# names(diff_files) <-
#   c("diff_KD", "diff_OE")

diff_files <- lapply(
  diff_files,
  makeGRangesFromDataFrame,
  keep.extra.columns =
    TRUE,
  na.rm = TRUE
)

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
        font.size = 7,
        title = "GO Enrichment Analysis Without Outliers and Low Peaks")

cluster_result_df <- compGO@compareClusterResult

# Extract genes ID of the desired/dominant enrichment class

lost_KD_df <-
  subset(cluster_result_df,
         Cluster == "Lost_KD" & Description == "clathrin binding")
# extract class
genes_in_major_class = lost_KD_df$geneID
genes_in_major_class
#which.min(lost_KD_df$p.adjust)
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

# diff_KD_anno <-
#   data.frame(peakAnnoList[["diff_KD"]]@anno)
# diff_OE_anno <- data.frame(peakAnnoList[["diff_OE"]]@anno)


# Get gene names from gene ids

mart <-
  useEnsembl(dataset = "mmusculus_gene_ensembl", biomart = "ensembl")
tmp_enriched_KD <- data.frame("entrezgene_id" = enriched_KD_anno$geneId)
tmp_enriched_OE <- data.frame("entrezgene_id" = enriched_OE_anno$geneId)
tmp_lost_KD <- data.frame("entrezgene_id" = lost_KD_anno$geneId)
tmp_lost_OE <- data.frame("entrezgene_id" = lost_OE_anno$geneId)

# tmp_diff_KD <- data.frame("entrezgene_id" = diff_KD_anno$geneId)
# tmp_diff_OE <- data.frame("entrezgene_id" = diff_OE_anno$geneId)

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

# ids_diff_KD <- getBM(
#   attributes = c("entrezgene_id", "external_gene_name"),
#   filters = "entrezgene_id",
#   values = tmp_diff_KD[, "entrezgene_id"],
#   mart = mart
# )
# ids_diff_OE <- getBM(
#   attributes = c("entrezgene_id", "external_gene_name"),
#   filters = "entrezgene_id",
#   values = tmp_diff_OE[, "entrezgene_id"],
#   mart = mart
# )

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
lost_KD_anno_major_class <- subset(lost_KD_anno, geneId %in% genes_in_major_class[[1]])

names(ids_lost_OE)[names(ids_lost_OE) == "entrezgene_id"] <- "geneId"
ids_lost_OE$geneId <- as.character(ids_lost_OE$geneId)
lost_OE_anno <- left_join(
  x = lost_OE_anno,
  y = ids_lost_OE,
  by = "geneId",
  relationship = "many-to-many"
)
lost_OE_anno <- lost_OE_anno[, -c(7:8)]

# names(ids_diff_KD)[names(ids_diff_KD) == "entrezgene_id"] <- "geneId"
# ids_diff_KD$geneId <- as.character(ids_diff_KD$geneId)
# diff_KD_anno <- left_join(
#   x = diff_KD_anno,
#   y = ids_diff_KD,
#   by = "geneId",
#   relationship = "many-to-many"
# )
# diff_KD_anno <- diff_KD_anno[, -c(7:8)]
#
# names(ids_diff_OE)[names(ids_diff_OE) == "entrezgene_id"] <- "geneId"
# ids_diff_OE$geneId <- as.character(ids_diff_OE$geneId)
# diff_OE_anno <- left_join(
#   x = diff_OE_anno,
#   y = ids_diff_OE,
#   by = "geneId",
#   relationship = "many-to-many"
# )
# diff_OE_anno <- diff_OE_anno[, -c(7:8)]

# Export differential peaks

diff_KD_anno <- rbind(enriched_KD_anno, lost_KD_anno)
write.table(
  diff_KD_anno,
  "differential_analysis/all_Tau_KD_anno_without_outliers_2k.csv",
  row.names = FALSE,
  quote = FALSE,
  sep = ";"
)
diff_KD_anno_gr <- GRanges(
  seqnames = diff_KD_anno$seqnames,
  ranges = IRanges(start = diff_KD_anno$start, end = diff_KD_anno$end)
)
export(
  diff_KD_anno_gr,
  "differential_analysis/diff_KD_masterPeak_sep_without_outliers.bed",
  "bed"
)
diff_KD_anno <- diff_KD_anno[c("PeakID", "external_gene_name")]

diff_OE_anno <- rbind(enriched_OE_anno, lost_OE_anno)
write.table(
  diff_OE_anno,
  "differential_analysis/all_Tau_OE_anno_without_outliers_2k.csv",
  row.names = FALSE,
  quote = FALSE,
  sep = ";"
)
diff_OE_anno_gr <- GRanges(
  seqnames = diff_OE_anno$seqnames,
  ranges = IRanges(start = diff_OE_anno$start, end = diff_OE_anno$end)
)
export(
  diff_OE_anno_gr,
  "differential_analysis/diff_OE_masterPeak_sep_without_outliers.bed",
  "bed"
)
diff_OE_anno <- diff_OE_anno[c("PeakID", "external_gene_name")]

# Merge gene name with statistics

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

# Volcano plot

p <- ggplot(
  data = diff_KD,
  aes(
    x = log2FoldChange,
    y = -log10(padj),
    col = diffexpressed,
    label = external_gene_name
  )
) +
  ggtitle("Tau Knockdown Differential Annotation") +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red") +
  geom_hline(yintercept = -log10(0.1), col = "red") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0, 7) +
  xlim(-5, 5)
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
  ggtitle("Tau Overexpression Differential Annotation") +
  geom_point() +
  theme_minimal() +
  geom_text_repel() +
  scale_color_manual(values = c("blue", "black", "red")) +
  geom_vline(xintercept = c(-0.5, 0.5), col = "red") +
  geom_hline(yintercept = -log10(0.1), col = "red") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylim(0, 7) +
  xlim(-9, 9)
p1


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
file <- read.table("SEACR/Tau_KD_1_top01_peaks.stringent.bed")

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
  fun = "enrichGO",
  #organism = "mouse",
  pvalueCutoff = 0.1,
  OrgDb = org.Mm.eg.db,
  pAdjustMethod = "BH"
)

dotplot(
  compKEGG,
  showCategory = 7,
  font.size = 8,
  title = "GO Enrichment Analysis"
)

