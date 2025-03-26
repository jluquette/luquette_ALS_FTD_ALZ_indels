#!/bin/bash

#if [ $# -lt 5 ]; then
    #echo "usage: $0 joint_depth_matrix.tab.gz bins.bed reference_genome.fasta.fai column_number_for_sample output.bw [chrom]"
    #exit 1
#fi

joint_depth_matrix="${snakemake_input[joint_matrix]}"
echo $joint_depth_matrix
#bins="${snakemake_input[1]}"
#echo $bins
fasta_index="${snakemake_input[fai]}"
echo $fasta_index
sample_name="${snakemake_wildcards[sample]}"
echo $sample_name
output="${snakemake_output[0]}"
echo $output
#chrom=$6
#tmpdir="${snakemake_resources[tmpdir]}"

sample_column=$(tabix -H $joint_depth_matrix | tr '\t' '\n' | awk '{{ if ($0 == "'$sample_name'") print NR; }}')

echo "sample_name=$sample_name -> sample_column=$sample_column"

#tabix $joint_depth_matrix $chrom \
gunzip -c $joint_depth_matrix \
    | tail -n +2 \
    | cut -f 1,2,$sample_column \
    | awk -v OFS='\t' '{ print $1, $2-1, $2, $3; }' \
    | bigtools bedgraphtobigwig --sorted start --nzooms 0 - $fasta_index $output

    #| bedtools map -a $bins -b /dev/stdin -o mean \
    #| sort --buffer-size=4G --temporary-directory= -k1,1 -k2,2n \
