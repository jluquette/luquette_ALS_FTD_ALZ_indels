#### Read and write SEA-AD DLPFC snRNA-seq data for Seurat ####

import anndata as ad

# Load h5ad
adata = ad.read_h5ad("SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.h5ad")
# Clean up sample meta data column names 
new_columns = ['sample_id', 'Neurotypical_reference', 'Donor_ID', 'Organism', 'Brain_Region', 'Sex', 'Gender', 'Age_at_Death', 'Race_choice_White', 'Race_choice_Black_African_American', 'Race_choice_Asian', 'Race_choice_American_Indian_Alaska_Native', 'Race_choice_Native Hawaiian_or_Pacific_Islander', 'Race_choice_Unknown_or_unreported', 'Race_choice_Other', 'specify other race', 'Hispanic/Latino', 'Highest level of education', 'Years of education', 'PMI', 'Fresh Brain Weight', 'Brain pH', 'Overall AD neuropathological Change', 'Thal', 'Braak', 'CERAD score', 'Overall CAA Score', 'Highest Lewy Body Disease', 'Total Microinfarcts_not observed grossly', 'Total microinfarcts in screening sections', 'Atherosclerosis', 'Arteriolosclerosis', 'LATE', 'Cognitive Status', 'Last CASI Score', 'Interval from last CASI in months', 'Last MMSE Score', 'Interval from last MMSE in months', 'Last MOCA Score', 'Interval from last MOCA in months', 'APOE Genotype', 'Primary Study Name', 'Secondary Study Name', 'cell_prep_type', 'facs_population_plan', 'rna_amplification', 'sample_name', 'sample_quantity_count', 'expc_cell_capture', 'method', 'pcr_cycles', 'percent_cdna_longer_than_400bp', 'rna_amplification_pass_fail', 'amplified_quantity_ng', 'load_name', 'library_prep', 'library_input_ng', 'r1_index', 'avg_size_bp', 'quantification_fmol', 'library_prep_pass_fail', 'exp_component_vendor_name', 'batch_vendor_name', 'experiment_component_failed', 'alignment', 'Genome', 'ar_id', 'bc', 'GEX_Estimated_number_of_cells', 'GEX_number_of_reads', 'GEX_sequencing_saturation', 'GEX_Mean_raw_reads_per_cell', 'GEX_Q30_bases_in_barcode', 'GEX_Q30_bases_in_read_2', 'GEX_Q30_bases_in_UMI', 'GEX_Percent_duplicates', 'GEX_Q30_bases_in_sample_index_i1', 'GEX_Q30_bases_in_sample_index_i2', 'GEX_Reads_with_TSO', 'GEX_Sequenced_read_pairs', 'GEX_Valid_UMIs', 'GEX_Valid_barcodes', 'GEX_Reads_mapped_to_genome', 'GEX_Reads_mapped_confidently_to_genome', 'GEX_Reads_mapped_confidently_to_intergenic_regions', 'GEX_Reads_mapped_confidently_to_intronic_regions', 'GEX_Reads_mapped_confidently_to_exonic_regions', 'GEX_Reads_mapped_confidently_to_transcriptome', 'GEX_Reads_mapped_antisense_to_gene', 'GEX_Fraction_of_transcriptomic_reads_in_cells', 'GEX_Total_genes_detected', 'GEX_Median_UMI_counts_per_cell', 'GEX_Median_genes_per_cell', 'Multiome_Feature_linkages_detected', 'Multiome_Linked_genes', 'Multiome_Linked_peaks', 'ATAC_Confidently_mapped_read_pairs', 'ATAC_Fraction_of_genome_in_peaks', 'ATAC_Fraction_of_high_quality_fragments_in_cells', 'ATAC_Fraction_of_high_quality_fragments_overlapping_TSS', 'ATAC_Fraction_of_high_quality_fragments_overlapping_peaks', 'ATAC_Fraction_of_transposition_events_in_peaks_in_cells', 'ATAC_Mean_raw_read_pairs_per_cell', 'ATAC_Median_high_quality_fragments_per_cell', 'ATAC_Non-nuclear_read_pairs', 'ATAC_Number_of_peaks', 'ATAC_Percent_duplicates', 'ATAC_Q30_bases_in_barcode', 'ATAC_Q30_bases_in_read_1', 'ATAC_Q30_bases_in_read_2', 'ATAC_Q30_bases_in_sample_index_i1', 'ATAC_Sequenced_read_pairs', 'ATAC_TSS_enrichment_score', 'ATAC_Unmapped_read_pairs', 'ATAC_Valid_barcodes', 'Number of mapped reads', 'Number of unmapped reads', 'Number of multimapped reads', 'Number of reads', 'Number of UMIs', 'Genes detected', 'Doublet score', 'Fraction mitochondrial UMIs', 'Used in analysis', 'Class confidence', 'Class', 'Subclass confidence', 'Subclass', 'Supertype confidence', 'Supertype_non-expanded', 'Supertype', 'Severely Affected Donor']
new_columns = [sub.replace(' ', '_') for sub in new_columns]
new_columns = [sub.replace('/', '_') for sub in new_columns]
adata.obs.columns = new_columns
# Subset by cell type to solve memory issue
import anndata as ad
from scipy.io import mmwrite
import pandas as pd
import os
subclasses = pd.unique(adata.obs.Subclass)
for subclass in subclasses:
	subclass_out = subclass.replace('/','_').replace(' ','_') 
	out_path = "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass/SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean." + subclass_out + ".h5ad"
	print(out_path)
	# Write h5ad (only X is read by Schard into Seurat)
	adata[adata.obs.Subclass == subclass].write(out_path)
# Write layer UMIs (raw count matrix) as MM for pseudobulk
file_list = os.listdir("SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass")
file_list = [file for file in file_list if file.endswith("h5ad")]
# for file in file_list:
for i in range(3,len(file_list)):
	file = file_list[i]
	subclass_out = file.split(".")[3]
	adata = ad.read_h5ad("SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass/" + file)
	out_path = "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass/SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean." + subclass_out + ".UMIs.mtx"
	print(out_path)
	mmwrite(out_path, adata.layers["UMIs"])

# Try subsetting L2/3 IT by supertypes and do pseudobulk DEG in each supertype. 
import anndata as ad
from scipy.io import mmwrite
import pandas as pd
import os
file = "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.L2_3_IT.h5ad"
adata = ad.read_h5ad("SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass/" + file)
supertypes = pd.unique(adata.obs.Supertype)
for supertype in supertypes:
	supertype_out = supertype.replace('/','_').replace(' ','_')
	print(supertype_out)
	out_path = "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass/SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean." + supertype_out + ".h5ad"
	# Write h5ad (only X is read by Schard into Seurat)
	adata[adata.obs.Supertype == supertype].write(out_path)
	# Write raw UMI count matrix
	out_path = "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass/SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean." + supertype_out + ".UMIs.mtx"
	mmwrite(out_path, adata[adata.obs.Supertype == supertype].layers["UMIs"])
# Subset L2/3 IT to DEG groups due to memory issue with Seurat
amp_meta = pd.read_csv("AD_microglia_Table_S4_mCA_sample_list_individual_metadata.csv")
control_donors = amp_meta.loc[amp_meta['Diagnosis']=='Healthy control', 'Individual ID']
ad_donors = amp_meta.loc[amp_meta['Diagnosis']=="Alzheimer's disease", 'Individual ID']
adata.obs['DEG_group'] = "Unknown"
adata.obs.loc[adata.obs['Donor_ID'].isin(control_donors), 'DEG_group'] = 'Control'
adata.obs.loc[adata.obs['Donor_ID'].isin(ad_donors), 'DEG_group'] = 'AD'
len(adata.obs.loc[adata.obs['DEG_group']=="Unknown",:]) 
adata.obs.loc[:,['Donor_ID','Overall_AD_neuropathological_Change','DEG_group']].drop_duplicates().reset_index(drop=True)
len(adata.obs.loc[adata.obs['DEG_group']=="Control",'Donor_ID'].drop_duplicates())
len(adata.obs.loc[adata.obs['DEG_group']=="AD",'Donor_ID'].drop_duplicates())
# Subset adata
sub = adata[adata.obs['DEG_group'] != "Unknown", :]
# Write h5ad and raw UMI
out_path = "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass/SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.L2_3_IT.DEG_AMPAD.h5ad"
# Write h5ad (only X is read by Schard into Seurat)
sub.write(out_path)
# Write raw UMI count matrix
out_path = "SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.subclass/SEAAD_DLPFC_RNAseq_final-nuclei.2024-02-13.clean.L2_3_IT.DEG_AMPAD.UMIs.mtx"
mmwrite(out_path, sub.layers["UMIs"])

