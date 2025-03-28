#!/bin/bash

if [ $# -ne 1 ]; then
    echo "usage: $0 new_env_name"
    exit 1
fi

new_env_name=$1

echo "note: r-seurat is listed here but does not actually work"

conda create -n $new_env_name -c conda-forge -c bioconda -c jluquette -c dranew -c soil 'bcftools>=1.19' bedtools bioconductor-annotationdbi bioconductor-annotationfilter bioconductor-biovizbase bioconductor-bsgenome bioconductor-bsgenome.hsapiens.1000genomes.hs37d5 bioconductor-bsgenome.hsapiens.ucsc.hg19 bioconductor-bsgenome.hsapiens.ucsc.hg38 bioconductor-genomeinfodb bioconductor-genomeinfodbdata bioconductor-genomicalignments bioconductor-genomicfeatures bioconductor-genomicranges bioconductor-gviz bioconductor-txdb.hsapiens.ucsc.hg19.knowngene bioconductor-variantannotation bx-python mscorefonts pandas r-argparse r-curl r-devtools r-extrafont r-extrafontdb r-harmonicmeanp r-kernsmooth r-lme4 r-lmertest r-nnls r-openxlsx r-pheatmap r-rcolorbrewer r-readxl r-rttf2pt1 r-svglite s3transfer scan2 ucsc-bedgraphtobigwig ucsc-bigwigaverageoverbed r-r.matlab bioconductor-dnacopy bioconductor-ctc r-inline r-fastcluster r-heatmap3 r-rio vcf2maf vcflib r-seurat r-mclust bigtools wiggletools sigprofilerassignment r-ggpmisc r-ggplot2 r-forcats r-ggsignif r-patchwork r-broom r-ggseqlogo


echo "INSTALL MANUALLY: r-mutenrich"
echo "Also install manually: SigProfilerExtractor. See https://github.com/AlexandrovLab/SigProfilerExtractor"
