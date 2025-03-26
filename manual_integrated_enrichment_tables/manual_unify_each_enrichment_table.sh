#!/bin/bash

mkdir -p quantile
mkdir -p bed_regions

for eclass in conservation dna_repair_hotspots gencode_simplified gtex gtex_expression_mc02 nott repliseq roadmap_chromhmm_brain roadmap_histone_signal_brain scatacseq scrnaseq_expression_mc02 fragile_sites gene_length; do
echo $eclass
    for etype in quantile bed_regions; do
echo $etype
        if [ -d ../enrichment/$eclass/$etype ]; then
echo "enrichment/$eclass/$etype"
            outfile="$etype/${eclass}.csv"
            if [ ! -f $outfile ]; then
                ./unify_enrichment_tables.R $outfile ../metadata/group_metadata_table2.csv $(ls ../enrichment/$eclass/$etype/*.csv)
            fi
        fi
    done
done


echo "removing nott quantiles (not wanted, but bed_regions wanted)"
rm quantile/nott.csv

# delete old stuff
rm */*.standardized.csv

# rearrange and rename columns into a standard format
for f in */*.csv; do
    echo $f
    ./standardize_enrichment_table.R $f ${f/.csv/.standardized.csv}
done


# combine all standardized tables into one
../snakemake/scripts/combine_uniform_tables.R enrichment.csv */*.standardized.csv
