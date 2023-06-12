#!/bin/bash

outfile=snakefile.ginkgo_per_donor

if [ -f $outfile ]; then
    echo "output file $outfile already exists, please delete it first"
    exit 1
fi


cat > $outfile << EOF
# vim: syntax=python
#
# It doesn't seem possible to make the standard scan2_call_mutations rule
# create all output .rda files (one per single cell) for each donor because
# the donor is a wildcard and outputs cannot depend on wildcards.  Further,
# since one invocation of SCAN2 creates all such outputs per donor, it would
# be incorrect to allow tools further down the pipeline to ask for each
# individual output file.
#
# Using modules here as a workaround.
#
# ------------------ THIS FILE IS AUTOMATICALLY GENERATED -------------------
EOF

donors=$(tail -n +2 ../metadata/panel_metadata_162cells_66ptaols_40mdaglia_56ptaneurons_21bulks_20brains.csv | cut -f1 -d, | sort | uniq)

for donor in $donors; do
    cat >> $outfile << EOF


module ginkgo_run_$donor:
    snakefile: "snakefile.ginkgo"

use rule ginkgo_run from ginkgo_run_${donor} as ginkgo_run_${donor} with:
    input:
        mapped_beds=lambda wildcards: expand("ginkgo/$donor/{sample}.bed_mapped",
            sample=list(bams['$donor']['single_cell'].keys()) + list(bams['$donor']['bulk'].keys())),
    log:
        "ginkgo/$donor/ginkgo.log"
    benchmark:
        "ginkgo/$donor/ginkgo_benchmark.txt"
    output:
        "ginkgo/$donor/CNV1",
        "ginkgo/$donor/CNV2",
        "ginkgo/$donor/SegBreaks",
        "ginkgo/$donor/SegCopy",
        "ginkgo/$donor/SegFixed",
        "ginkgo/$donor/SegNorm",
        "ginkgo/$donor/SegStats",
        "ginkgo/$donor/data",
        "ginkgo/$donor/results.txt",
        "ginkgo/$donor/status.xml",
        temp("ginkgo/$donor/ploidyDummy.txt"),
        temp("ginkgo/$donor/refDummy.bed_mapped"),
        plots=protected(expand('ginkgo/$donor/{sample}_{plottype}.jpeg',
            sample=list(bams['$donor']['single_cell'].keys()) + list(bams['$donor']['bulk'].keys()),
            plottype=[ 'CN', 'GC', 'SoS', 'counts', 'dist', 'hist', 'lorenz' ]))
    params:
        workdir="ginkgo/$donor"

ruleorder: ginkgo_run_${donor} > ginkgo_run

EOF
done >> $outfile
