library(IsoformSwitchAnalyzeR)

########## Import Data #############

stringTieQuant <- importIsoformExpression(
  parentDir = "/home/haeun/DATA/splicing/Analysis/Quant/TCGA-BRCA/ballgown-strict",
  readLength =75
)

myDesign <- data.frame(
    sampleID = colnames(salmonQuant$abundance)[-1],
    condition = gsub('_.*', '', colnames(salmonQuant$abundance)[-1])
)
myDesign
#>        sampleID   condition
#> 1 Fibroblasts_0 Fibroblasts
#> 2 Fibroblasts_1 Fibroblasts
#> 3        hESC_0        hESC
#> 4        hESC_1        hESC
#> 5         iPS_0         iPS
#> 6         iPS_1         iPS

### Create switchAnalyzeRlist
aSwitchList <- importRdata(
    isoformCountMatrix   = salmonQuant$counts,
    isoformRepExpression = salmonQuant$abundance,
    designMatrix         = myDesign,
    isoformExonAnnoation = system.file("extdata/example.gtf.gz"             , package="IsoformSwitchAnalyzeR"),
    isoformNtFasta       = system.file("extdata/example_isoform_nt.fasta.gz", package="IsoformSwitchAnalyzeR"),
    showProgress = FALSE
)
aSwitchList
#> This switchAnalyzeRlist list contains:
#>  1092 isoforms from 362 genes
#>  3 comparison from 3 conditions (in total 6 samples)
#>
#> Feature analyzed:
#> [1] "ntSequence"
