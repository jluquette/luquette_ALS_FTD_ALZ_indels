#!/bin/bash

cat ../metadata/immutable_metadata.csv | grep -v '^donor' \
| while read line; do
    sample=$(echo $line | cut -f2 -d,)
    bam=$(echo $line | cut -f4 -d,)
    bai=$(ls ${bam/.bam/}*.bai)

    echo ln -s $bam ${sample}.bam
    echo ln -s $bai ${sample}.bai
done
