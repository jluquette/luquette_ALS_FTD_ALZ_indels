#!/bin/bash

if [ $# -ne 3 ]; then
    echo "usage: $0 binsize_string binsize_number output_prefix"
    exit 1
fi

bs_string=$1
bs_number=$2
out_prefix=$3


# Autosomes
windows=$(bedtools makewindows -w $bs_number -b hg19_autosomes.bed | bedtools subtract -a /dev/stdin -b unanalyzable_chr2_3kb_region_2\:33139000-33141948.bed | wc -l)
auto_out=${out_prefix}_chr1-22_${windows}windows_${bs_string}.3kb_unanalyzable_region_removed.txt
if [ -f $auto_out ]; then
    echo "output file $auto_out already exists, please delete it first"
    exit 1
else
    echo "making $auto_out"
    bedtools makewindows -w $bs_number -b hg19_autosomes.bed \
    | bedtools subtract -a /dev/stdin -b unanalyzable_chr2_3kb_region_2\:33139000-33141948.bed \
    > $auto_out
fi



# chrX
windows=$(bedtools makewindows -w $bs_number -b hg19_chrX_no_PARs.bed | wc -l)
x_out=${out_prefix}_chrX_noPARs_${windows}windows_${bs_string}.txt
if [ -f $x_out ]; then
    echo "output file $x_out already exists, please delete it first"
    exit 1
else
    echo "making $x_out"
    bedtools makewindows -w $bs_number -b hg19_chrX_no_PARs.bed > $x_out
fi



# chrY
windows=$(bedtools makewindows -w $bs_number -b hg19_chrY_no_PARs.bed | wc -l)
y_out=${out_prefix}_chrY_noPARs_${windows}windows_${bs_string}.txt
if [ -f $y_out ]; then
    echo "output file $y_out already exists, please delete it first"
    exit 1
else
    echo "making $y_out"
    bedtools makewindows -w $bs_number -b hg19_chrY_no_PARs.bed > $y_out
fi



# final file
windows=$(wc -l $auto_out $x_out $y_out | tail -n 1 |sed -e's/ total$//' | tr -d ' ')
final_out=${out_prefix}_chr1-22XY_${windows}windows_${bs_string}.3kb_unanalyzable_region_removed.txt
if [ -f $final_out ]; then
    echo "output file $final_out already exists, please delete it first"
    exit 1
else
    echo "making $final_out"
    cat $auto_out $x_out $y_out > $final_out
fi
