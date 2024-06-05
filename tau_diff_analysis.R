install.packages("tidyverse")
install.packages("RColorBrewer")
install.packages("ggrepel")

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("DESeq2")
BiocManager::install("chromVAR")
BiocManager::install("rtracklayer")
BiocManager::install('plyranges')
BiocManager::install('IHW')
BiocManager::install('edgeR')


library(DESeq2)
library(dplyr)
library(chromVAR)
library(rtracklayer)
library(plyranges)
library(tidyverse)
library(RColorBrewer)
library(ggrepel)
library(IHW)
library(edgeR)

projPath = "~/Charite Thesis/H3K9ac Cut&Tag/Results/"
setwd(projPath)

repL = c("1", "2")
#condL = c("Control", "Tau_KD", "Tau_OE")
condL_KD = c("Control", "Tau_KD")
condL_OE = c("Control", "Tau_OE")


# Create a master peak list merging all the peaks called for each sample

#mPeak = GRanges()
mPeak_KD = GRanges()
mPeak_OE = GRanges()

# for(cond in condL){
#   for(rep in repL){
#     if (cond=="Tau_OE" & rep=="2") {
#       break
#     }
#     peakRes = read.table(paste0(projPath, "SEACR/", cond, "_", rep, "_top01_peaks.stringent.bed"), header = FALSE, fill = TRUE)
#     mPeak = GRanges(seqnames = peakRes$V1, IRanges(start = peakRes$V2, end = peakRes$V3), strand = "*") %>% append(mPeak, .)
#   }
# }
# masterPeak = reduce(mPeak)


for (cond in condL_KD) {
  for (rep in repL) {
    peakRes = read.table(
      paste0(
        projPath,
        "SEACR/",
        cond,
        "_",
        rep,
        "_top01_peaks.stringent.bed"
      ),
      header = FALSE,
      fill = TRUE
    )
    mPeak_KD = GRanges(
      seqnames = peakRes$V1,
      IRanges(start = peakRes$V2, end = peakRes$V3),
      strand = "*"
    ) %>% append(mPeak_KD, .)
  }
}
masterPeak_KD = GenomicRanges::reduce(mPeak_KD)
export(masterPeak_KD,
       "differential_analysis/KD_masterPeak.bed",
       "bed")

for (cond in condL_OE) {
  for (rep in repL) {
    if (cond == "Tau_OE" & rep == "2") {
      break
    }
    peakRes = read.table(
      paste0(
        projPath,
        "SEACR/",
        cond,
        "_",
        rep,
        "_top01_peaks.stringent.bed"
      ),
      header = FALSE,
      fill = TRUE
    )
    mPeak_OE = GRanges(
      seqnames = peakRes$V1,
      IRanges(start = peakRes$V2, end = peakRes$V3),
      strand = "*"
    ) %>% append(mPeak_OE, .)
  }
}
masterPeak_OE = GenomicRanges::reduce(mPeak_OE)
export(masterPeak_OE,
       "differential_analysis/OE_masterPeak.bed",
       "bed")


################ clean bed file from heatmap clustering ####################

heatmapKD <- import("heatmap/KD_H3K9ac_sorted_regions_4kmeans.bed", format =
                      "BED")
heatmapKD
avoid_clusters = c('cluster_1')
heatmapKD_clean <- heatmapKD %>% filter(NA. %in% avoid_clusters)
heatmapKD_clean
KD_masterPeak <- import("differential_analysis/KD_masterPeak_without_outliers.bed",
                        format = "BED")
KD_masterPeak
KD_masterPeak_clean <- KD_masterPeak %>% filter_by_non_overlaps(heatmapKD_clean)
KD_masterPeak_clean

export(
  KD_masterPeak_clean,
  "differential_analysis/KD_masterPeak_4k_without_outliers_and_c1.bed",
  "bed"
)

heatmapOE <- import("heatmap/OE_H3K9ac_sorted_regions_4kmeans.bed", format =
                      "BED")
heatmapOE
avoid_clusters = c('cluster_1')
heatmapOE_clean <- heatmapOE %>% filter(NA. %in% avoid_clusters)
heatmapOE_clean
OE_masterPeak <- import("differential_analysis/OE_masterPeak_without_outliers.bed",
                        format = "BED")
OE_masterPeak
OE_masterPeak_clean <- OE_masterPeak %>% filter_by_non_overlaps(heatmapOE_clean)
OE_masterPeak_clean

export(
  OE_masterPeak_clean,
  "differential_analysis/OE_masterPeak_4k_without_outliers_and_c1.bed",
  "bed"
)


################## start here with clean #################


masterPeak_KD <- import("differential_analysis/KD_masterPeak_4k_without_outliers_and_c1.bed",
                        format = "BED")
masterPeak_OE <- import("differential_analysis/OE_masterPeak_4k_without_outliers_and_c1.bed",
                        format = "BED")

# Get the fragment counts for each peak in the master peak list

bamDir = paste0(projPath, "bam")
#countMat = matrix(NA, length(masterPeak), 5)
countMat_KD = matrix(NA, length(masterPeak_KD), 4)
countMat_OE = matrix(NA, length(masterPeak_OE), 3)

# Overlap with bam file to get count

bamFile_C1 = paste0(bamDir, "/", "Control_1.sorted.bam")
bamFile_C2 = paste0(bamDir, "/", "Control_2.sorted.bam")
bamFile_KD1 = paste0(bamDir, "/", "Tau_KD_1.sorted.bam")
bamFile_KD2 = paste0(bamDir, "/", "Tau_KD_2.sorted.bam")
bamFile_OE1 = paste0(bamDir, "/", "Tau_OE_1.sorted.bam")

# Counts from a unified KD/OE masterPeak

# bamFile_C1 = paste0(bamDir, "/", "Control_1.sorted.bam")
# fragment_counts_C1 <- getCounts(bamFile_C1, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
# countMat[, 1] = counts(fragment_counts_C1)[,1]
#
# bamFile_C2 = paste0(bamDir, "/", "Control_2.sorted.bam")
# fragment_counts_C2 <- getCounts(bamFile_C2, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
# countMat[, 2] = counts(fragment_counts_C2)[,1]
#
# bamFile_KD1 = paste0(bamDir, "/", "Tau_KD_1.sorted.bam")
# fragment_counts_KD1 <- getCounts(bamFile_KD1, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
# countMat[, 3] = counts(fragment_counts_KD1)[,1]
#
# bamFile_KD2 = paste0(bamDir, "/", "Tau_KD_2.sorted.bam")
# fragment_counts_KD2 <- getCounts(bamFile_KD2, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
# countMat[, 4] = counts(fragment_counts_KD2)[,1]
#
# bamFile_OE1 = paste0(bamDir, "/", "Tau_OE_1.sorted.bam")
# fragment_counts_OE1 <- getCounts(bamFile_OE1, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
# countMat[, 5] = counts(fragment_counts_OE1)[,1]

# code efficient but memory would crash so decided to separate
# i = 1
# for(cond in condL){
#   for(rep in repL){
#     if (cond=="Tau_OE" & rep=="2") {
#       break
#     }
#     bamFile = paste0(bamDir, "/", cond, "_", rep, ".sorted.bam")
#     fragment_counts <- getCounts(bamFile, masterPeak, paired = TRUE, by_rg = FALSE, format = "bam")
#     countMat[, i] = counts(fragment_counts)[,1]
#     i = i + 1
#   }
# }

# Counts with KD

fragment_counts_C1_KD <- getCounts(
  bamFile_C1,
  masterPeak_KD,
  paired = TRUE,
  by_rg = FALSE,
  format = "bam"
)
countMat_KD[, 1] = counts(fragment_counts_C1_KD)[, 1]

fragment_counts_C2_KD <- getCounts(
  bamFile_C2,
  masterPeak_KD,
  paired = TRUE,
  by_rg = FALSE,
  format = "bam"
)
countMat_KD[, 2] = counts(fragment_counts_C2_KD)[, 1]

fragment_counts_KD1_KD <- getCounts(
  bamFile_KD1,
  masterPeak_KD,
  paired = TRUE,
  by_rg = FALSE,
  format = "bam"
)
countMat_KD[, 3] = counts(fragment_counts_KD1_KD)[, 1]

fragment_counts_KD2_KD <- getCounts(
  bamFile_KD2,
  masterPeak_KD,
  paired = TRUE,
  by_rg = FALSE,
  format = "bam"
)
countMat_KD[, 4] = counts(fragment_counts_KD2_KD)[, 1]


# Counts with OE

fragment_counts_C1_OE <- getCounts(
  bamFile_C1,
  masterPeak_OE,
  paired = TRUE,
  by_rg = FALSE,
  format = "bam"
)
countMat_OE[, 1] = counts(fragment_counts_C1_OE)[, 1]

fragment_counts_C2_OE <- getCounts(
  bamFile_C2,
  masterPeak_OE,
  paired = TRUE,
  by_rg = FALSE,
  format = "bam"
)
countMat_OE[, 2] = counts(fragment_counts_C2_OE)[, 1]

fragment_counts_OE1_OE <- getCounts(
  bamFile_OE1,
  masterPeak_OE,
  paired = TRUE,
  by_rg = FALSE,
  format = "bam"
)
countMat_OE[, 3] = counts(fragment_counts_OE1_OE)[, 1]


# colnames(countMat) = c("Control_1", "Control_2", "Tau_KD_1", "Tau_KD_2", "Tau_OE_1")
colnames(countMat_KD) = c("Control_1", "Control_2", "Tau_KD_1", "Tau_KD_2")
colnames(countMat_OE) = c("Control_1", "Control_2", "Tau_OE_1")

#write.csv(countMat, "countMat.csv", row.names = FALSE)
write.csv(countMat_KD,
          "countMat_KD_4k_without_outliers_and_c1.csv",
          row.names = FALSE)
write.csv(countMat_OE,
          "countMat_OE_4k_without_outliers_and_c1.csv",
          row.names = FALSE)

# verify import is ok
# countMat2 = read.table("countMat.csv", sep=",", header=TRUE)

############### START HERE when countMat created ##############################

countMat_KD <- read.table("countMat_KD_4k_without_outliers_and_c1.csv",
                          sep = ",",
                          header = TRUE)
countMat_OE <- read.table("countMat_OE_4k_without_outliers_and_c1.csv",
                          sep = ",",
                          header = TRUE)

# Sequencing depth normalization and differential enriched peaks detection

# selectR = which(rowSums(countMat) > 5) ## remove low count genes
# dataS = countMat[selectR,]
selectR_KD = which(rowSums(countMat_KD) > 10) ## remove low count genes
dataS_KD = countMat_KD[selectR_KD, ]
selectR_OE = which(rowSums(countMat_OE) > 10) ## remove low count genes
dataS_OE = countMat_OE[selectR_OE, ]

# Control vs KD
#dataS1 = dataS[,0:4]
condition_KD = factor(c("Control", "Control", "Tau_KD", "Tau_KD"))
dds_KD = DESeqDataSetFromMatrix(
  countData = dataS_KD,
  colData = DataFrame(condition_KD),
  design = ~ condition_KD
)
DDS_KD = DESeq(dds_KD)
normDDS_KD = counts(DDS_KD, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS_KD) = paste0(colnames(normDDS_KD), "_norm")
res_KD = results(DDS_KD,
                 independentFiltering = FALSE,
                 altHypothesis = "greaterAbs")

countMatDiff_KD = cbind(dataS_KD, normDDS_KD, res_KD)
head(countMatDiff_KD)

# Some plots
plotMA(res_KD, ylim = c(-10, 10))
dds <- DESeqTransform(DDS_KD)
plotPCA(dds, intgroup = "condition_KD")
boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)
plotDispEsts(DDS_KD)

# The significant differentially expressed genes are the ones found in the upper-left and upper-right corners.
# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2FoldChange respectively positive or negative)

# add a column of NAs
countMatDiff_KD$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and padj < 0.1, set as "UP"
countMatDiff_KD$diffexpressed[countMatDiff_KD$log2FoldChange > 0.5 &
                                countMatDiff_KD$padj < 0.1] <- "UP"
# if log2Foldchange < -0.5 and padj < 0.1, set as "DOWN"
countMatDiff_KD$diffexpressed[countMatDiff_KD$log2FoldChange < -0.5 &
                                countMatDiff_KD$padj < 0.1] <- "DOWN"

# Volcano plot to check threshold
p <- ggplot(data = countMatDiff_KD, aes(
  x = log2FoldChange,
  y = -log10(padj),
  col = diffexpressed
)) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold
p2 <- p + geom_vline(xintercept = c(-0.5, 0.5), col = "red") +
  geom_hline(yintercept = -log10(0.1), col = "red")
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3


# Basic scatterplot

countMatDiff_KD[, c('Control_1_norm', 'Control_2_norm')]
countMatDiff_KD$mean_control <- rowMeans(countMatDiff_KD[, c('Control_1_norm', 'Control_2_norm')], na.rm =
                                           TRUE)
countMatDiff_KD$mean_KD <- rowMeans(countMatDiff_KD[, c('Tau_KD_1_norm', 'Tau_KD_2_norm')], na.rm =
                                      TRUE)

ggplot(countMatDiff_KD, aes(x = mean_control, y = mean_KD)) +
  geom_point()

# Merge differential peaks with masterPeak

KD_masterPeak_prob <- merge(
  x = data.frame(masterPeak_KD),
  y = countMatDiff_KD,
  by = "row.names",
  all = TRUE
)
names(KD_masterPeak_prob)[1] <- ('PeakID')
write.table(
  KD_masterPeak_prob,
  "differential_analysis/KD_masterPeak_prob_sep_2k_without_outliers_padj0.1.csv",
  row.names = FALSE,
  quote = FALSE
)


# Enriched/lost peaks

# enriched_KD_genes = which(countMatDiff_KD$padj < 0.1 & countMatDiff_KD$log2FoldChange > 0.5)
# countMatDiff_KD[enriched_KD_genes,]
#
# lost_KD_genes = which(countMatDiff_KD$padj < 0.1 & countMatDiff_KD$log2FoldChange < -0.5)
# countMatDiff_KD[lost_KD_genes,]

# enriched_KD_masterPeak <- masterPeak[enriched_KD_genes]
# export(enriched_KD_masterPeak, "differential_analysis/enriched_KD_masterPeak.bed", "bed")
#
# lost_KD_masterPeak <- masterPeak[lost_KD_genes]
# export(lost_KD_masterPeak, "differential_analysis/lost_KD_masterPeak.bed", "bed")
# to compare : import("differential_analysis/lost_KD_masterPeak_sep_without_outliers.bed", format="BED")


# Control vs OE

#skipping KD columns
#dataS2 = dataS[,-c(3:4)]

condition_OE = factor(c("Control", "Control", "Tau_OE"))

#same length as dataS_OE so does not change anything
keep <- filterByExpr(dataS_OE, design = condition_OE)

dds_OE = DESeqDataSetFromMatrix(
  countData = dataS_OE,
  colData = DataFrame(condition_OE),
  design = ~ condition_OE
)
DDS_OE = DESeq(dds_OE)
normDDS_OE = counts(DDS_OE, normalized = TRUE) ## normalization with respect to the sequencing depth
colnames(normDDS_OE) = paste0(colnames(normDDS_OE), "_norm")
res_OE = results(DDS_OE,
                 independentFiltering = FALSE,
                 altHypothesis = "greaterAbs")
summary(res_OE)
resultsNames(DDS_OE)

countMatDiff_OE = cbind(dataS_OE, normDDS_OE, res_OE)
head(countMatDiff_OE)

# Some plots
plotMA(res_OE, ylim = c(-10, 10))
dds <- DESeqTransform(DDS_OE)
plotPCA(dds, intgroup = "condition_OE")
boxplot(log10(assays(dds)[["cooks"]]), range = 0, las = 2)
plotDispEsts(DDS_OE)

# Doesn't increase stat power
resIHW <- results(DDS_OE, filterFun = ihw)
summary(resIHW)



# add a column of NAs
countMatDiff_OE$diffexpressed <- "NO"
# if log2Foldchange > 0.5 and padj < 0.1, set as "UP"
countMatDiff_OE$diffexpressed[countMatDiff_OE$log2FoldChange > 0.5 &
                                countMatDiff_OE$padj < 0.1] <- "UP"
# if log2Foldchange < -0.5 and padj < 0.1, set as "DOWN"
countMatDiff_OE$diffexpressed[countMatDiff_OE$log2FoldChange < -0.5 &
                                countMatDiff_OE$padj < 0.1] <- "DOWN"

# Volcano plot
p <- ggplot(data = countMatDiff_OE, aes(
  x = log2FoldChange,
  y = -log10(padj),
  col = diffexpressed
)) + geom_point() + theme_minimal()
# Add vertical lines for log2FoldChange thresholds, and one horizontal line for the p-value threshold
p2 <- p + geom_vline(xintercept = c(-0.5, 0.5), col = "red") +
  geom_hline(yintercept = -log10(0.1), col = "red")
mycolors <- c("blue", "red", "black")
names(mycolors) <- c("DOWN", "UP", "NO")
p3 <- p2 + scale_colour_manual(values = mycolors)
p3

# Merge differential peaks with masterPeak

OE_masterPeak_prob <- merge(
  x = data.frame(masterPeak_OE),
  y = countMatDiff_OE,
  by = "row.names",
  all = TRUE
)
names(OE_masterPeak_prob)[1] <- ('PeakID')
write.table(
  OE_masterPeak_prob,
  "differential_analysis/OE_masterPeak_prob_sep_2k_without_outliers_padj0.1.csv",
  row.names = FALSE,
  quote = FALSE
)


# Enriched/lost peaks

# enriched_OE_genes = which(countMatDiff_OE$padj < 0.1 & countMatDiff_OE$log2FoldChange > 0.5)
# countMatDiff_OE[enriched_OE_genes,]
#
# lost_OE_genes = which(countMatDiff_OE$padj < 0.1 & countMatDiff_OE$log2FoldChange < -0.5)
# countMatDiff_OE[lost_OE_genes,]

# enriched_OE_masterPeak <- masterPeak[enriched_OE_genes]
# export(enriched_OE_masterPeak, "differential_analysis/enriched_OE_masterPeak.bed", "bed")
#
# lost_OE_masterPeak <- masterPeak[lost_OE_genes]
# export(lost_OE_masterPeak, "differential_analysis/lost_OE_masterPeak.bed", "bed")


############### Convert bed data from Ensembl to UCSC convention needed for annotation #######################

convert_seqnames <- function(ensembl_dataframe, conversion_dataframe) {
  #colnames(ensembl_dataframe) <- c('chrom', 'chromStart', 'chromEnd')
  names(ensembl_dataframe)[2:4] <- c('chr', 'start', 'end')
  # Indices in conversion_dataframe that match ensembl_dataframe chr column values
  indices <- match(ensembl_dataframe$chr, conversion_dataframe$ensembl)
  # Replace "ensembl" column values in df1 by corresponding values from "ucsc" column in df2
  ensembl_dataframe$chr <- conversion_dataframe$ucsc[indices]
  #ucsc_gr <- GRanges(seqnames = ensembl_dataframe$chrom,
  #ranges = IRanges(start = ensembl_dataframe$chromStart, end = ensembl_dataframe$chromEnd))
  
  return(ensembl_dataframe)
}

conversion = read.table("GRCm39_ensembl2UCSC.txt", col.names = c("ensembl", "ucsc"))

# Specify folder containing bed Ensembl files to convert
ensembl_bed_folder <- "differential_analysis/"

# Specify folder for outputs
ucsc_bed_folder <- "differential_analysis/"

# List of bed Ensembl files in folder
ensembl_files <- list.files(ensembl_bed_folder,
                            pattern = "prob_sep_2k_without_outliers_padj0.1.csv$",
                            full.names = TRUE)
ensembl_files

# Loop on each file
for (ensembl_file in ensembl_files) {
  # Create output file name
  ucsc_file <-  file.path(ucsc_bed_folder,
                          paste0(tools::file_path_sans_ext(basename(ensembl_file)), ".ucsc.csv"))
  # Upload bed Ensembl file
  #ensembl_data <- read.table(ensembl_file, header = FALSE, colClasses = c("character", "numeric", "numeric", rep("NULL", 3)))
  ensembl_data <- read.table(ensembl_file, header = TRUE)
  # Convert
  ucsc_bed <- convert_seqnames(ensembl_data, conversion)
  #ucsc_bed <- subset(ucsc_bed, select = -c(width))
  # Save output
  write.table(
    ucsc_bed,
    ucsc_file,
    sep = ',',
    row.names = FALSE,
    quote = FALSE
  )
  #export(ucsc_bed, ucsc_file, "bed")
  #print(import(ucsc_file))
}
