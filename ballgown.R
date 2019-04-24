library(ballgown)
library(RSkittleBrewer)
library(genefilter)
library(dplyr)
library(devtools)

bg = ballgown(dataDir = "/home/haeun/DATA1/splicing/ballgown-allgene_394", 
	      samplePattern = "TCGA", meas= "all")

sample_info = read.delim("/home/haeun/DATA/splicing/Analysis/HR_leafviz/BRCA_free/groups/BRCA_free.HRD_sig3_median_groups.extended.txt", 
                         sep='\t', header = F, row.names = 1, check.names = F)
pData(bg) = sample_info


# filter out low-abundance genes
bg_filt = subset(bg,"rowVars(texpr(bg)) >1",genomesubset=TRUE)   
results_transcripts = stattest(bg_filt, feature="transcript",
			       covariate="Group", getFC=TRUE, meas="FPKM")

# add gene name
results_transcripts = data.frame(geneNames=ballgown::geneNames(bg_filt), 
                                 geneIDs=ballgown::geneIDs(bg_filt), 
                                 transcriptNames=ballgown::transcriptNames(bg_filt), 
                                 results_transcripts)    

indices = match(results_transcripts$id, texpr(bg_filt, 'all')$gene_id)
gene_names_for_result = texpr(bg_filt, 'all')$gene_name[indices]
results_transcripts2 = data.frame(geneNames=gene_names_for_result, results_transcripts)
final = subset(results_transcripts2 ,results_transcripts$qval<0.05)     # filtering
write.csv(final, "DET_results.csv", row.names=FALSE)


plotTranscripts(gene='MSTRG.67928', gown=bg, samples='TCGA-D8-A1JM-01A-11R-A13Q-07', 
                meas='cov', colorby='transcript', 
                main='transcripts from gene XLOC_000454: sample 12, FPKM')

plotMeans('MSTRG.67928', bg, groupvar='V2', meas='cov', colorby='transcript')
clusterTranscripts(gene='MSTRG.33696', gown=bg, k=5)
