# SPLICING
Haeun Hwangbo 2019-06-17  
나중에 파일 찾는데 조금이나마 도움이 되길 바라며...  
## NOTE

BRCA-free = BRCA-active  
BRCA-event = BRCA-inactive  
HRD high = HRD  
HRD low = Non-HRD

## FILE LOCATION

### Overview

haeun@143.248.31.113  
password: hello0512  

splicing/  
├── **Analysis** : Analyzed by Haeun (StringTie Quantification, DEXSeq results, etc.)  
├── **data** : downloaded files (GENCODE, KEGG, DepMap, etc.)   
└── **src** : source codes


### RNA-seq BAM Files
- TCGA-BRCA
    - haeun@143.248.31.113:/home/haeun/DATA/splicing/data/TCGA/TCGA-BRCA/bamfiles/ (n=85)
    - omics@143.248.31.178:/home/omics/DATA2/04_hyeongu/TCGA_BAM/bamfiles/ (n=418)
    - omics@143.248.31.223:/home/omics/DATA_GPU1/06_haeun/splicing/data/TCGA-BRCA/bamfiles/ (n=107)
    - wrongly downloaded files (duplicated, normal tissue, metastatic) aer in *to-be-removed* folders or tar.gz files
- CCLE-BRCA
    - omics@143.248.31.178:/home/omics/DATA2/04_hyeongu/Cell-line/CCLE-BREAST/bamfiles/ (n=56)
- PDC
    - omics@143.248.31.178:/home/omics/DATA2/04_hyeongu/PDX_RNAseq/bamfiles/ (n=24)
- TCGA-OV
    - haeun@143.248.31.113:/home/haeun/DATA1/splicing/TCGA-OV/bamfiles/ (n=76)
    - omics@143.248.31.223:/home/omics/DATA_GPU1/06_haeun/splicing/data/TCGA-OV/bamfiles/ (n=70)
    - Note: median 기준으로 HRD / non-HRD에 속하는 sample들만 다운받았다.
- CCLE-OV
    - omics@143.248.31.178:/home/omics/DATA1/06_haeun/splicing/data/CCLE/RNAseq/CCLE-OV/bamfiles/ (n=45)
    - Note: bamfiles.old directory contains bamfiles before changing to chr~ notation (여기는 chr nubmer가 그냥 1 2 3 처럼 되어있다.)

### Target Genes (Tier)

Merged KEGG Replication and repair pathways and curated set from [here](https://www.mdanderson.org/documents/Labs/Wood-Laboratory/human-dna-repair-genes.html)  
- Main File  
    - haeun@143.248.31.113:/home/haeun/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt
    - omics@143.248.31.178:/home/omics/DATA1/06_haeun/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt (duplicate)


### StringTie (Transcript-level Quantification)

- Source Code
    - omics@143.248.31.178:/home/omics/DATA1/06_haeun/splicing/src/quant/
- Annotation File
    - GENCODE v29 (TCGA, PDC), v19 (CCLE)
    - omics@143.248.31.178:/home/omics/DATA1/06_haeun/splicing/data/GENCODE/
- TCGA-BRCA
    - omics@143.248.31.178:/home/omics/DATA1/06_haeun/splicing/Analysis/Quant/TCGA-BRCA/
- CCLE-BRCA
    - omics@143.248.31.178:/home/omics/DATA1/06_haeun/splicing/Analysis/Quant/CCLE-BRCA/
- PDC
    - omics@143.248.31.178:/home/omics/DATA1/06_haeun/splicing/Analysis/Quant/PDC/
- TCGA-OV
    - omics@143.248.31.178:/home/omics/DATA1/06_haeun/splicing/Analysis/Quant/TCGA-OV.tar.gz
- CCLE-OV
    - omics@143.248.31.178:/home/omics/DATA1/06_haeun/splicing/Analysis/Quant/CCLE-OV.tar.gz

#### Directory structure

├── **assembled** : stringtie assembly  
├── **ballgown** : old-version quantification (merge -i, -f 0.05, -a5)  
├── **ballgown-strict** : new-version quantification (stringtie default)  
├── **gffcompare** : gffcompare results  
└── **merged** : GTF and tmap_extended files used for ballgown/ballgown-strict and downstream analyses


### EdgeR (Gene-level Quantification & DEG)

- Source Code
    - haeun@143.248.31.113:/home/haeun/DATA/splicing/src/DEG_analysis.R
- Result
    - haeun@143.248.31.113:/home/haeun/DATA/splicing/Analysis/DEG/


### Ballgown (DET)

- Source Code
    - haeun@143.248.31.113:/home/haeun/DATA/splicing/src/ballgown.R
    - haeun@143.248.31.113:/home/haeun/DATA/splicing/src/annotate_DET_result.py
- Result
    - haeun@143.248.31.113:/home/haeun/DATA/splicing/Analysis/DET/allgene_DET_results.tsv
    - haeun@143.248.31.113:/home/haeun/DATA/splicing/Analysis/DET/allgene_DET_results_tier12.tsv


### DEXSeq

Differential exon usage between HRD and Non-HRD in BRCA-active  group (TCGA-BRCA)  
Only tested tier 123 genes (tier123 gene만으로만 테스트했기 때문에 결과 정확하지 않을 수 있음)

- Source Code
    - haeun@143.248.31.113:/home/haeun/DATA/splicing/src/DEXSeq.R
    - haeun@143.248.31.113:/home/haeun/DATA/splicing/src/DEXSeq_results.R
- Result
    - haeun@143.248.31.113:/home/omics/DATA1/06_haeun/splicing/Analysis/DEXSeq/
    - omics@143.248.31.178:/home/omics/DATA1/06_haeun/splicing/Analysis/DEXSeq/ (duplicate)

### Leafcutter

- Source Code
    - haeun@143.248.31.113:/home/omics/DATA1/06_haeun/splicing/Analysis/HR_leafviz/
- Result
    - haeun@143.248.31.113:/home/omics/DATA1/06_haeun/splicing/Analysis/HR_leafviz/


### Cell Lines

haeun@143.248.31.113:/home/haeun/DATA/splicing/data/Cell_lines  
├── **CCLE/**  
│   ├── CCLE_sample_info_file_2012-10-18.txt  
│   ├── Cell_lines_annotations_20181226.txt  
│   ├── CNV/  
│   ├── drug/  
│   ├── README  
│   ├── RNAseq/  
│   ├── RPPA/  
│   └── variant/ : vcfs and mutalisk results  
├── Cell_lines_CCLE_COSMIC_merged_info.txt  
├── **COSMIC/**  
│   ├── CellLinesCodingMuts.vcf  
│   ├── CellLinesNonCodingVariants.vcf  
│   ├── CosmicCLP_MutantExport.tsv  
│   ├── CosmicSample.tsv  
│   ├── genotypes  
│   ├── QC.xlsx  
│   └── README  
├── **DepMap/**  
│   ├── D2_BRCA_combined_gene_dep_scores_minus.tsv  
│   ├── D2_BRCA_combined_gene_dep_scores.tsv  
│   ├── D2_BRCA_sample_info.tsv  
│   ├── D2_combined_gene_dep_scores.csv  
│   ├── D2_combined_gene_dep_score_SDs.csv  
│   ├── D2_combined_gene_dep_scores.tsv  
│   ├── D2_OV_combined_gene_dep_scores.tsv  
│   ├── D2_OV_sample_info.tsv  
│   ├── D2_README.txt  
│   ├── D2_sample_info.csv  
│   ├── DepMap_18q3_README.txt  
│   ├── DepMap_18Q3_sample_info.csv  
│   ├── DepMap_18QC_gene_depenency.csv  
│   ├── DepMap-2018q3-celllines.csv  
│   ├── public_18Q3_sample_info.cnv  
│   └── README  
└── **GDSC/**  
    ├── Cell_Lines_Details.xlsx  
    ├── GDSC-CCLE-CTRP_conversion.xlsx  
    ├── GDSC_Fitted_Data_Description.pdf  
    ├── GDSC_Raw_Data_Description.pdf  
    ├── PARP-inhibitor  
    ├── v17.3_fitted_dose_response.tsv  
    ├── v17.3_fitted_dose_response.xlsx  
    └── v17.3_public_raw_data.tsv  
