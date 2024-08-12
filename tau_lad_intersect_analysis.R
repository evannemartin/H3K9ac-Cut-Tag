if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("regioneR")
BiocManager::install("BSgenome.Mmusculus.UCSC.mm39") # used by RegioneR

library(regioneR)
library(BSgenome.Mmusculus.UCSC.mm39)

projPath = "~/Charite Thesis/H3K9ac Cut&Tag/Results/"
setwd(projPath)

############################ Permutation test #################################

# Import
T1_LAD_peaks <- import("../LAD domains/T1_GSM5669198_D5Midbrain_LADs_mm39.bed",
                       format = "BED")
T2_LAD_peaks <- import("../LAD domains/T2_GSM5669198_D5Midbrain_LADs_mm39.bed",
                       format = "BED")
enriched_KD_peaks <- import(
  "differential_analysis/enriched_KD_masterPeak_sep_without_outliers_padj0.1.ucsc.bed",
  format = "BED"
)
lost_KD_peaks <- import(
  "differential_analysis/lost_KD_masterPeak_sep_without_outliers_padj0.1.ucsc.bed",
  format = "BED"
)
enriched_OE_peaks <- import(
  "differential_analysis/enriched_OE_masterPeak_sep_without_outliers_padj0.1.ucsc.bed",
  format = "BED"
)
lost_OE_peaks <- import(
  "differential_analysis/lost_OE_masterPeak_sep_without_outliers_padj0.1.ucsc.bed",
  format = "BED"
)

len_T1_LAD_peaks <- length(T1_LAD_peaks)
len_T2_LAD_peaks <- length(T2_LAD_peaks)
len_enriched_KD_peaks <- length(enriched_KD_peaks)
len_lost_KD_peaks <- length(lost_KD_peaks)
len_enriched_OE_peaks <- length(enriched_OE_peaks)
len_lost_OE_peaks <- length(lost_OE_peaks)

# Filter to have only canonical
T1_LAD_peaks <- filterChromosomes(T1_LAD_peaks, organism = "mm39", chr.type =
                                    "canonical")
T2_LAD_peaks <- filterChromosomes(T2_LAD_peaks, organism = "mm39", chr.type =
                                    "canonical")
enriched_KD_peaks <- filterChromosomes(enriched_KD_peaks, organism = "mm39", chr.type =
                                         "canonical")
lost_KD_peaks <- filterChromosomes(lost_KD_peaks, organism = "mm39", chr.type =
                                     "canonical")
enriched_OE_peaks <- filterChromosomes(enriched_OE_peaks, organism = "mm39", chr.type =
                                         "canonical")
lost_OE_peaks <- filterChromosomes(lost_OE_peaks, organism = "mm39", chr.type =
                                     "canonical")

### Permutation test with T1

enriched_KD_analysis <- overlapPermTest(
  A = enriched_KD_peaks,
  B = T1_LAD_peaks,
  ntimes = 500,
  genome = "mm39",
  count.once = TRUE,
  mask = NA
)
enriched_KD_analysis
plot(enriched_KD_analysis)
mean(enriched_KD_analysis$numOverlaps$permuted)

lost_KD_analysis <- overlapPermTest(
  A = lost_KD_peaks,
  B = T1_LAD_peaks,
  ntimes = 500,
  genome = "mm39",
  count.once = TRUE,
  mask = NA
)
lost_KD_analysis
plot(lost_KD_analysis)
mean(lost_KD_analysis$numOverlaps$permuted)

enriched_OE_analysis <- overlapPermTest(
  A = enriched_OE_peaks,
  B = T1_LAD_peaks,
  ntimes = 500,
  genome = "mm39",
  count.once = TRUE,
  mask = NA
)
enriched_OE_analysis
plot(enriched_OE_analysis)
mean(enriched_OE_analysis$numOverlaps$permuted)

lost_OE_analysis <- overlapPermTest(
  A = lost_OE_peaks,
  B = T1_LAD_peaks,
  ntimes = 500,
  genome = "mm39",
  count.once = TRUE,
  mask = NA
)
lost_OE_analysis
plot(lost_OE_analysis)
mean(lost_OE_analysis$numOverlaps$permuted)

### Permutation test with T2

enriched_KD_analysis_T2 <- overlapPermTest(
  A = enriched_KD_peaks,
  B = T2_LAD_peaks,
  ntimes = 500,
  genome = "mm39",
  count.once = TRUE,
  mask = NA
)
enriched_KD_analysis_T2
plot(enriched_KD_analysis_T2)
mean(enriched_KD_analysis_T2$numOverlaps$permuted)

lost_KD_analysis_T2 <- overlapPermTest(
  A = lost_KD_peaks,
  B = T2_LAD_peaks,
  ntimes = 500,
  genome = "mm39",
  count.once = TRUE,
  mask = NA
)
lost_KD_analysis_T2
plot(lost_KD_analysis_T2)
mean(lost_KD_analysis_T2$numOverlaps$permuted)

enriched_OE_analysis_T2 <- overlapPermTest(
  A = enriched_OE_peaks,
  B = T2_LAD_peaks,
  ntimes = 500,
  genome = "mm39",
  count.once = TRUE,
  mask = NA
)
enriched_OE_analysis_T2
plot(enriched_OE_analysis_T2)
mean(enriched_OE_analysis_T2$numOverlaps$permuted)

lost_OE_analysis_T2 <- overlapPermTest(
  A = lost_OE_peaks,
  B = T2_LAD_peaks,
  ntimes = 500,
  genome = "mm39",
  count.once = TRUE,
  mask = NA
)
lost_OE_analysis_T2
plot(lost_OE_analysis_T2)
mean(lost_OE_analysis_T2$numOverlaps$permuted)
