setwd('DATA/splicing/Analysis/DEXSeq')

library(DEXSeq)
load('BRCA_free.HRD_sig3_median_groups.tier_123.dexseq.Rdata')
plotMA(res, cex=0.8)

gene_mapping = read.delim('/home/haeun/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_list.ensembl.txt', 
                          header = F, row.names = 2, colClasses = "character")
for (i in unique(res$groupID[which(res$padj < 0.1)])) {
  try({
    pdf(paste(gene_mapping[i, 'V1'], '_DEXSeq.pdf', sep=""), width = 5, height = 8)
    plotDEXSeq(res, i, displayTranscripts = T, legend = T,
               cex.axis=1.2, cex=1, lwd=1.5)
    dev.off()
  })
}

# gene level result
tmp <- cbind(gene_id = names(pgq), gene_name = gene_mapping[names(pgq), ], "adjP" = pgq)
write.table(tmp,
            file = paste0("BRCA_free.HRD_sig3_median_groups.tier_123.dexseq_gene_level", ".txt"), 
            col.names = T, row.names = F, quote = F, sep = "\t")

# exon level result
write.table(res,
            file = paste0("BRCA_free.HRD_sig3_median_groups.tier_123.dexseq_exon_level", ".txt"),
            col.names = T, row.names = T, quote = F, sep = "\t")
