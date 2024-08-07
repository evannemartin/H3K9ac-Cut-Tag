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
862/len_enriched_KD_peaks

pvalue_enriched_KD <- enriched_KD_analysis$numOverlaps$pval

# Compute Odds Ratio

random_overlap_nb_enriched_KD <- mean(enriched_KD_analysis$numOverlaps$permuted)
observed_overlap_nb_enriched_KD <-  enriched_KD_analysis$numOverlaps$observed
random_nonoverlap_nb_enriched_KD <- len_enriched_KD_peaks-random_overlap_nb_enriched_KD
observed_nonoverlap_nb_enriched_KD <- len_enriched_KD_peaks-observed_overlap_nb_enriched_KD

OR_enriched_KD = (observed_overlap_nb_enriched_KD / random_overlap_nb_enriched_KD)/(observed_nonoverlap_nb_enriched_KD /
                                                                                      random_nonoverlap_nb_enriched_KD)
log(OR_enriched_KD)


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
1082/len_lost_KD_peaks

pvalue_lost_KD <- lost_KD_analysis$numOverlaps$pval

random_overlap_nb_lost_KD <- mean(lost_KD_analysis$numOverlaps$permuted)
observed_overlap_nb_lost_KD <-  lost_KD_analysis$numOverlaps$observed
random_nonoverlap_nb_lost_KD <- len_lost_KD_peaks-random_overlap_nb_lost_KD
observed_nonoverlap_nb_lost_KD <- len_lost_KD_peaks-observed_overlap_nb_lost_KD

OR_lost_KD = (observed_overlap_nb_lost_KD / random_overlap_nb_lost_KD)/(observed_nonoverlap_nb_lost_KD /
                                                                          random_nonoverlap_nb_lost_KD)
log(OR_lost_KD)


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
472/len_enriched_OE_peaks

pvalue_enriched_OE <- enriched_OE_analysis$numOverlaps$pval

random_overlap_nb_enriched_OE <- mean(enriched_OE_analysis$numOverlaps$permuted)
observed_overlap_nb_enriched_OE <-  enriched_OE_analysis$numOverlaps$observed
random_nonoverlap_nb_enriched_OE <- len_enriched_OE_peaks-random_overlap_nb_enriched_OE
observed_nonoverlap_nb_enriched_OE <- len_enriched_OE_peaks-observed_overlap_nb_enriched_OE

OR_enriched_OE = (observed_overlap_nb_enriched_OE / random_overlap_nb_enriched_OE)/(observed_nonoverlap_nb_enriched_OE /
                                                                                      random_nonoverlap_nb_enriched_OE)
log(OR_enriched_OE)


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
933/len_lost_OE_peaks

pvalue_lost_OE <- lost_OE_analysis$numOverlaps$pval

random_overlap_nb_lost_OE <- mean(lost_OE_analysis$numOverlaps$permuted)
observed_overlap_nb_lost_OE <-  lost_OE_analysis$numOverlaps$observed
random_nonoverlap_nb_lost_OE <- len_lost_OE_peaks-random_overlap_nb_lost_OE
observed_nonoverlap_nb_lost_OE <- len_lost_OE_peaks-observed_overlap_nb_lost_OE

OR_lost_OE = (observed_overlap_nb_lost_OE / random_overlap_nb_lost_OE)/(observed_nonoverlap_nb_lost_OE /
                                                                          random_nonoverlap_nb_lost_OE)
log(OR_lost_OE)


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
# 862/len_enriched_KD_peaks
#
# pvalue_enriched_KD <- enriched_KD_analysis$numOverlaps$pval
#
# random_overlap_nb_enriched_KD <- mean(enriched_KD_analysis$numOverlaps$permuted)
# observed_overlap_nb_enriched_KD <-  enriched_KD_analysis$numOverlaps$observed
# random_nonoverlap_nb_enriched_KD <- len_enriched_KD_peaks-random_overlap_nb_enriched_KD
# observed_nonoverlap_nb_enriched_KD <- len_enriched_KD_peaks-observed_overlap_nb_enriched_KD
#
# OR_enriched_KD = (observed_overlap_nb_enriched_KD / random_overlap_nb_enriched_KD)/(observed_nonoverlap_nb_enriched_KD /
#                                                                                       random_nonoverlap_nb_enriched_KD)
# log(OR_enriched_KD)


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
# 1082/len_lost_KD_peaks
#
# pvalue_lost_KD <- lost_KD_analysis$numOverlaps$pval
#
# random_overlap_nb_lost_KD <- mean(lost_KD_analysis$numOverlaps$permuted)
# observed_overlap_nb_lost_KD <-  lost_KD_analysis$numOverlaps$observed
# random_nonoverlap_nb_lost_KD <- len_lost_KD_peaks-random_overlap_nb_lost_KD
# observed_nonoverlap_nb_lost_KD <- len_lost_KD_peaks-observed_overlap_nb_lost_KD
#
# OR_lost_KD = (observed_overlap_nb_lost_KD / random_overlap_nb_lost_KD)/(observed_nonoverlap_nb_lost_KD /
#                                                                           random_nonoverlap_nb_lost_KD)
# log(OR_lost_KD)


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
# 472/len_enriched_OE_peaks
#
# pvalue_enriched_OE <- enriched_OE_analysis$numOverlaps$pval
#
# random_overlap_nb_enriched_OE <- mean(enriched_OE_analysis$numOverlaps$permuted)
# observed_overlap_nb_enriched_OE <-  enriched_OE_analysis$numOverlaps$observed
# random_nonoverlap_nb_enriched_OE <- len_enriched_OE_peaks-random_overlap_nb_enriched_OE
# observed_nonoverlap_nb_enriched_OE <- len_enriched_OE_peaks-observed_overlap_nb_enriched_OE
#
# OR_enriched_OE = (observed_overlap_nb_enriched_OE / random_overlap_nb_enriched_OE)/(observed_nonoverlap_nb_enriched_OE /
#                                                                                       random_nonoverlap_nb_enriched_OE)
# log(OR_enriched_OE)


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
# 933/len_lost_OE_peaks
#
# pvalue_lost_OE <- lost_OE_analysis$numOverlaps$pval
#
# random_overlap_nb_lost_OE <- mean(lost_OE_analysis$numOverlaps$permuted)
# observed_overlap_nb_lost_OE <-  lost_OE_analysis$numOverlaps$observed
# random_nonoverlap_nb_lost_OE <- len_lost_OE_peaks-random_overlap_nb_lost_OE
# observed_nonoverlap_nb_lost_OE <- len_lost_OE_peaks-observed_overlap_nb_lost_OE
#
# OR_lost_OE = (observed_overlap_nb_lost_OE / random_overlap_nb_lost_OE)/(observed_nonoverlap_nb_lost_OE /
#                                                                           random_nonoverlap_nb_lost_OE)
# log(OR_lost_OE)
