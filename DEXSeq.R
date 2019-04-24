# https://github.com/markrobinsonuzh/diff_splice_paper/blob/master/software/Rcode/dexseq_function_filter.R

library(DEXSeq)
library(BiocParallel)
BPPARAM = MulticoreParam(workers=4)

countFiles <- list.files("/home/haeun/DATA/splicing/Analysis/Quant-all/htseq_tier123_HRD_sig3_median", 
                         pattern="htseq_tier123.txt$", full.names = T)
flattenedFile <- list.files("/home/haeun/DATA/splicing/data/GENCODE", pattern="gencode.v29.protein_coding.DEXSeq.gff", full.names = T)
sampleTable <- read.delim("/home/haeun/DATA/splicing/Analysis/HR_leafviz/BRCA_free/groups/BRCA_free.HRD_sig3_median_groups.txt",
                          header = F, row.names = 1)
colnames(sampleTable) <- c("condition")

suppressPackageStartupMessages( library( "DEXSeq" ) )
dxd = DEXSeqDataSetFromHTSeq(
  countFiles,
  sampleData=sampleTable,
  design= ~ sample + exon + condition:exon,
  flattenedfile=flattenedFile )

# Normalization for defferent seqeuncing depth
dxd = estimateSizeFactors( dxd )
print(dim(featureCounts(dxd)))

filter_bin_count <- 50

## Filter by normalized bin counts
if (!is.null(filter_bin_count)) {
  norm_counts <- featureCounts(dxd, normalized = TRUE)
  keep <- which(rowSums(norm_counts) >= filter_bin_count)
  dxd <- dxd[keep, ]
}

print(dim(featureCounts(dxd)))
print(mean(colSums(featureCounts(dxd))))


# Dispersion estimation
dxd = estimateDispersions( dxd, BPPARAM=BPPARAM)
plotDispEsts( dxd )

dxd = testForDEU( dxd, BPPARAM=BPPARAM)
dxd = estimateExonFoldChanges(dxd, BPPARAM=BPPARAM)

message("Summarizing results on gene level...")
res <- DEXSeqResults(dxd)
pgq <- perGeneQValue(res, p = "pvalue")

## Save results
out_basename = "/home/haeun/DATA/splicing/Analysis/HR_leafviz/BRCA_fBRCA_free.HRD_sig3_median_groups.txt"
save(dxd, res, pgq, file = paste0(out_basename, ".Rdata"))
tmp <- cbind(gene = names(pgq), "adjP" = pgq)
colnames(tmp)[colnames(tmp) == "adjP"] <- paste0(method_name, ":adjP")
write.table(tmp, 
            file = paste0(out_basename, ".txt"), 
            col.names = TRUE, row.names = FALSE, quote = FALSE, sep = "\t")