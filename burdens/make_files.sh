#!/bin/bash

# Create SNV and indel burden long table with no metadata
cut -d, -f6,12- ../aging_rates/mutation_burdens_long_table.csv > mutburden_long_table.csv

# Create de novo signature burden long table
# (no longer a separate command - was make_sigburden_longtable.sh)
../snakemake/scripts/annotate_burdens_with_exposures.R \
    sigburden_long_table.csv \
    ../aging_rates/mutation_burdens_long_table.csv \
    denovo,SBS96,snv\;denovo,ID83,indel\;denovo,ID83_corrected,indel \
    ../mutsigs/sigprofilerextractor/all/SBS96/Most_Stab_Sigs/exposures.csv \
    ../mutsigs/sigprofilerextractor/all/ID83/Most_Stab_Sigs/exposures.csv \
    ../mutsigs/sigprofilerextractor/all/ID83_corrected/Most_Stab_Sigs/exposures.csv \
    # can optionally add COSMIC fits from SigProfilerExtractor, though they are not used
    #denovo,SBS96,snv\;denovo,ID83,indel\;denovo,ID83_corrected,indel\;cosmic,SBS96,snv\;cosmic,ID83,indel\;cosmic,ID83_corrected,indel \
    #../mutsigs/sigprofilerextractor/all/SBS96/Most_Stab_Sigs/cosmic_exposures.csv \
    #../mutsigs/sigprofilerextractor/all/ID83/Most_Stab_Sigs/cosmic_exposures.csv \
    #../mutsigs/sigprofilerextractor/all/ID83_corrected/Most_Stab_Sigs/cosmic_exposures.csv 


# Concatenate the two tables
../snakemake/scripts/combine_uniform_tables.R combined_burden_long_table.csv mutburden_long_table.csv sigburden_long_table.csv

# Compute expected burdens from controls
./model_expected_burdens.R combined_residualized_burdens.csv ../metadata/sample_metadata.csv combined_burden_long_table.csv 
