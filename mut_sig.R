library(BSgenome)
# download and load reference genome
ref_genome <- "BSgenome.Hsapiens.UCSC.hg38"
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(ref_genome, version = "3.8")
library(ref_genome, character.only = TRUE)
library(MutationalPatterns)


setwd('/home/haeun/DATA/splicing/data/TCGA/TCGA-BRCA/mutect2')
vcf_files <- list.files(pattern = ".vcf$", full.names = TRUE)
sample_names <- substring(vcf_files, 3, 14)
vcfs <- read_vcfs_as_granges(vcf_files, sample_names, ref_genome)
summary(vcfs)

muts = mutations_from_vcf(vcfs[[1]])
types = mut_type(vcfs[[1]])
context = mut_context(vcfs[[1]], ref_genome)
type_context = type_context(vcfs[[1]], ref_genome)
type_occurrences <- mut_type_occurrences(vcfs, ref_genome)
mut_mat <- mut_matrix(vcf_list = vcfs, ref_genome = ref_genome)
colnames(mut_mat) <- sample_names
#mut_mat <- mut_mat + 0.0001
#drop <- c('TCGA-AN-A046')
#mut_mat_wo4 <- mut_mat[ , -which(colnames(mut_mat) %in% drop)]

setwd('/home/haeun/DATA/splicing/Analysis/Signature/MutationalPatterns/')
write.table(mut_mat, "TCGA-BRCA.96subs_matrix.txt", sep="\t", quote=F)
write.table(type_occurrences, "TCGA-BRCA.6subs_matrix.txt", sep='\t', quote=F)
#####################################################
# De novo mutational signature extraction using NMF #
################################################### #

# “...a common way of deciding on the rank is to try different values, compute some quality
# measure of the results, and choose the best value according to this quality criteria. The most
# common approach is to choose the smallest rank for which cophenetic correlation coefficient
# starts decreasing. Another approach is to choose the rank for which the plot of the residual
# sum of squares (RSS) between the input matrix and its estimate shows an inflection point.”
library("NMF")
estimate <- nmf(mut_mat, rank=2:15, method="brunet", nrun=30, seed=123456)
plot(estimate)

for (num in 2:10){
  nmf_res <- extract_signatures(mut_mat, rank = num, nrun = 200)
  colnames(nmf_res$signatures) <- c(1:num)
  rownames(nmf_res$contribution) <- c(1:num)
  print(num)
  # calculate cosine similarity
  cos_sim_ori_rec <- cos_sim_matrix(mut_mat, nmf_res$reconstructed)
  cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
  print(sprintf("mean: %f median: %f", mean(cos_sim_ori_rec[, 1]), median(cos_sim_ori_rec[, 1])))
  relative_signatures <- apply(nmf_res$signatures, 2, function(x) x / sum(x) )
  for (i in 1:num){
    for (k in 1:30){

      if (cos_sim(relative_signatures[,i], cancer_signatures[,k]) > 0.8){
        print(sprintf("my %i - cosmic %i : %f", i, k, cos_sim(relative_signatures[,i], cancer_signatures[,k])))
      }
    }
  }
  
  rel_contribution <- apply(nmf_res$contribution, 2, function(x) x / sum(x) )
  abs_contribution <- nmf_res$contribution * colSums(nmf_res$signatures)
  # recurrently active mutations
  print("relative contribution")
  print(rowSums(abs_contribution) / sum(abs_contribution))
  print("abs contribution")
  print(rowSums(rel_contribution))
  print("***************************")
}


cos_sim_ori_rec <- cos_sim_matrix(mut_mat, nmf_res$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))


##################################
# Compare with COSMIC signatures #
##################################

# Download mutational signatures from the COSMIC website:
sp_url <- paste("https://cancer.sanger.ac.uk/cancergenome/assets/",
                "signatures_probabilities.txt", sep = "")
cancer_signatures = read.table(sp_url, sep = "\t", header = TRUE)
# Match the order of the mutation types to MutationalPatterns standard
new_order = match(row.names(mut_mat), cancer_signatures$Somatic.Mutation.Type)
# Reorder cancer signatures dataframe
cancer_signatures = cancer_signatures[as.vector(new_order),]
# Add trinucletiode changes names as row.names
row.names(cancer_signatures) = cancer_signatures$Somatic.Mutation.Type
# Keep only 96 contributions of the signatures in matrix
cancer_signatures = as.matrix(cancer_signatures[,4:33])


# Fit mutation matrix to the COSMIC mutational signatures:
select <- which(rowSums(fit_res2$contribution) > 10000)
consensus <- c(1, 2, 3, 5, 6, 8, 13, 17, 18, 26, 30)

fit_res <- fit_to_signatures(mut_mat, cancer_signatures[, consensus])
fit_res2 <- fit_to_signatures(mut_mat, cancer_signatures)
# Select signatures with some contribution

cos_sim_ori_rec <- cos_sim_matrix(mut_mat, fit_res$reconstructed)
cos_sim_ori_rec <- as.data.frame(diag(cos_sim_ori_rec))
cos_sim_ori_rec2 <- cos_sim_matrix(mut_mat, fit_res2$reconstructed)
cos_sim_ori_rec2 <- as.data.frame(diag(cos_sim_ori_rec2))

# total number of mutations per siganture
total_signatures = colSums(cancer_signatures)
# calculate signature contribution in absolute number of signatures
abs_contribution = fit_res2$contribution * total_signatures
# melt matrix
m_contribution = melt(abs_contribution)
colnames(m_contribution) = c("Signature", "Sample", "Contribution")

write.table(apply(fit_res2$contribution, 2, function(x) x / sum(x) ), 'TCGA-BRCA.all_cosimc.rel_contribution.txt', sep="\t", quote=F)
write.table(fit_res2$contribution, 'TCGA-BRCA.all_cosimc.abs_contribution.txt', sep="\t", quote=F)
write.table(cos_sim_ori_rec2, 'TCGA-BRCA.all_cosimc.ori_rec_cos_sim.txt', sep="\t", quote=F)
write.table(apply(fit_res$contribution, 2, function(x) x / sum(x) ), 'TCGA-BRCA.consensus_cosimc.rel_contribution.txt', sep="\t", quote=F)
write.table(fit_res$contribution, 'TCGA-BRCA.consensus_cosimc.abs_contribution.txt', sep="\t", quote=F)
write.table(cos_sim_ori_rec, 'TCGA-BRCA.consensus_cosimc.ori_rec_cos_sim.txt', sep="\t", quote=F)