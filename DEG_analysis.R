library(edgeR)
library(stringr)
library(tidyr)


setwd('/home/haeun/DATA/splicing/')
# import data

#sample_info <- read.delim('Analysis/HR_leafviz/BRCA_free/groups/BRCA_free.HRD_sig3_median_groups.txt', 
#                          header = F, check.names = F)
#clinical_info <- read.delim('Analysis/DEG/BRCA_free.HRD_sig3_median_groups.clinical_info.txt', 
#                            header = 1, check.names = F, sep='\t')
#colnames(sample_info) <- c("sample", "group")
#colnames(clinical_info)[1] <- "sample"
#merged_info <- merge(sample_info, clinical_info, by="sample")
#colnames(merged_info) <- c("sample", "HRD", "ER", "PR", "HER2", "stage", "age")
#expr <- read.delim('Analysis/DEG/BRCA_free.HRD_sig3_median_groups.protein_coding.htseq.txt', 
#                   row.names = 1, header=TRUE, check.names = FALSE)
#expr <- expr %>% drop_na()


sample_info <- read.delim('Analysis/DEG/sample_list/all_438_samples.txt', 
                          header = F, check.names = F)
clinical_info <- read.delim('Analysis/DEG/all_438.HRD_sig3_median_groups.clinical_info.txt', 
                            header = F, check.names = F, sep='\t')
colnames(sample_info) <- c("sample", "HRD", "BRCAness")
colnames(clinical_info)[1] <- "sample"
merged_info <- merge(sample_info, clinical_info, by="sample")
colnames(merged_info) <- c("sample", "HRD", "BRCAness","ER", "PR", "HER2", "stage", "age")


expr <- read.delim('data/Cell_lines/CCLE/RNAseq/CCLE-OV.read_counts.tsv', 
                   row.names = 1, header=TRUE, check.names = FALSE)
expr <- expr %>% drop_na()

genes<- read.delim('data/GENCODE/gencode.v19.protein_coding.id_name_mapping.txt', sep='\t', header=F, row.names = 1, check.names = F)
# filter out non-expressed genes (> 10 reads in >=10 samples)
filter <- apply(expr, 1, function(x) length(x[x>10])>=10)
filtered <- expr[filter,]
genes <- subset(genes,rownames(genes) %in% rownames(filtered))

# merge two dataframes
filtered$gene <- genes$V2[match(rownames(filtered),rownames(genes))]


###############################
# edgeR: DEG analysis
###############################

# BRCA-free/high vs. BRCA-event/high
# BRCA-free/low vs. BRCA-event/high

# make design matrix
Group <- factor(paste(merged_info$HRD, merged_info$BRCAness, sep="."))
HRD <- factor(merged_info$HRD)
ER <- factor(merged_info$ER)
HER2 <- factor(merged_info$HER2)
design <- model.matrix(~0+Group+ER+HER2)
colnames(design)[1:4] <- levels(Group)

my.contrasts <- makeContrasts(
  high.eventvsfree = high.BRCA_event - high.BRCA_free,
  low.eventvsfree = low.BRCA_event - low.BRCA_free,
  event.highvslow = high.BRCA_event - low.BRCA_event,
  free.highvslow = high.BRCA_free - low.BRCA_free,
  levels=design[, 1:4])

y <- DGEList(counts=subset(filtered,select=-gene), genes = filtered$gene)
y <- calcNormFactors(y)

# filter out lowly expressed genes
keep <- rowSums(cpm(y) > 1) >= 10
y <- y[keep, , keep.lib.sizes=F]
y <- calcNormFactors(y)

# MDS plot prior to test
mds <- plotMDS(y, labels=merged_info$HRD, col=as.numeric(merged_info$HER2))
mdslegend("bottomleft", legend=levels(merged_info$HER2), col=as.numeric(levels(merged_info$HER2)))

# estimate dispersion
y <- estimateDisp(y, design, robust=TRUE)
plotBCV(y)

# test DEG
fit <- glmQLFit(y, design)
for (contrast in c("high.eventvsfree", "low.eventvsfree", "event.highvslow", "free.highvslow")) {
  qlf <- glmQLFTest(fit, contrast=my.contrasts[, contrast])
  # save result
  print(summary(decideTests(qlf)))
  tab <- topTags(qlf, n=20000, adjust.method = 'fdr', p.value = 0.05)
  write.table(tab, str_interp("Analysis/DEG/${contrast}.HRD_sig3_median_groups.DEG_list.er_her2_covariate.fdr.txt"), sep="\t", quote = F)
}

write.table(cpm(y, normalized.lib.sizes = T, log = T), "Analysis/DEG/CCLE-OV.log_cpm.txt", sep="\t", quote=F)
write.table(cpm(y, normalized.lib.sizes = T, log = F), "Analysis/DEG/CCLE-OV.cpm.txt", sep="\t", quote=F)
