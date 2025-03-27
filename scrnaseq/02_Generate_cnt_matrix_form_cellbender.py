#from cell bender to cnt_matrix

import os
import re
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import scanpy as sc
from scipy.io import mmwrite
from matplotlib.backends.backend_pdf import PdfPages

from dir_util import dir_data, dir_out, create_dir

class OurGroupProcessor:
    def __init__(self):
        self.data_dir = os.path.join(dir_data, "our_group")
        self.cellbender_dir = os.path.join(dir_out, "02_cellbender_our_group")
        self.df = pd.read_csv(os.path.join(self.data_dir, "sample_name.tsv"), sep='\t')

    def filter_by_sample_name(self, sample_name):
        return self.df[self.df["our_sample_name"] == sample_name]

    def find_cellbender_output_path(self, sample_name):
        return os.path.join(self.cellbender_dir, sample_name, "output_filtered.h5")

    def filter_by_region_and_group(self, region, group_contains):
        region_bool = self.df["brain_region"].str.contains(region)
        group_bool = self.df["case_status"].str.contains(group_contains)
        return self.df[region_bool & group_bool]

    def generate_cellbender_output_list(self, filtered_df):
        sample_list = filtered_df["our_sample_name"].tolist()
        return [self.find_cellbender_output_path(sample) for sample in sample_list]


class SingleCellProcessor:
    def __init__(self, h5_file_path):
        self.h5_file_path = h5_file_path
        self.sample_name = self.extract_sample_name(h5_file_path)
        self.data = self.load_data()

    def extract_sample_name(self, file_path):
        match = re.search(r'/([^/]+)/output_filtered\\.h5$', file_path)
        return match.group(1) if match else None

    def load_data(self):
        data = sc.read_10x_h5(self.h5_file_path)
        data.var_names_make_unique()
        return data

    def perform_qc(self, data):
        data.var["mt"] = data.var_names.str.startswith("MT-")
        data.var["ribo"] = data.var_names.str.startswith(("RPS", "RPL"))
        data.var["hb"] = data.var_names.str.contains("^HB[^(P)]")
        sc.pp.calculate_qc_metrics(data, qc_vars=["mt", "ribo", "hb"], inplace=True, percent_top=[20], log1p=True)
        return data

    def filter_cells(self, data, min_genes=700, max_mt_pct=20):
        outliers = (data.obs['n_genes_by_counts'] < min_genes) | (data.obs['pct_counts_mt'] > max_mt_pct)
        return data[~outliers].copy()

    def normalize_the_sample(self, data):
        data_copy = data.copy()
        sc.pp.normalize_total(data_copy, target_sum=1e4)
        sc.pp.log1p(data_copy)
        data.obs['normalized_log1p_total_counts'] = np.sum(data_copy.X, axis=1)

    def process_data(self):
        self.data = self.perform_qc(self.data)
        self.data_filtered = self.filter_cells(self.data)
        self.normalize_the_sample(self.data)
        self.normalize_the_sample(self.data_filtered)
        return self.data, self.data_filtered

    def plot_qc_metrics(self, adata, sample_name):
        n_obs, n_vars = adata.shape
        fig, axs = plt.subplots(2, 2, figsize=(10, 10))
        fig.suptitle(f"Quality Control Metrics for {sample_name}\n{n_obs} cells × {n_vars} genes", fontsize=16)

        sns.histplot(adata.obs["total_counts"], bins=100, kde=False, ax=axs[0, 0])
        axs[0, 0].set_title("Total Counts Distribution")

        sns.histplot(adata.obs["normalized_log1p_total_counts"], bins=100, kde=False, ax=axs[0, 1])
        axs[0, 1].set_title("Total Counts Distribution (log1p)")

        sc.pl.violin(adata, "pct_counts_mt", ax=axs[1, 0], show=False)
        axs[1, 0].set_title("Percentage of Mitochondrial Counts")

        scatter = axs[1, 1].scatter(adata.obs["total_counts"], adata.obs["n_genes_by_counts"],
                                    c=adata.obs["pct_counts_mt"], cmap='viridis', s=10)
        axs[1, 1].set_xlabel("Total Counts")
        axs[1, 1].set_ylabel("Number of Genes by Counts")
        axs[1, 1].set_title("Total Counts vs Number of Genes")

        fig.colorbar(scatter, ax=axs[1, 1], orientation='vertical').set_label("% Mitochondrial Counts")
        plt.tight_layout()
        plt.close(fig)
        return fig

    def save_counts_matrix_mm(self, adata, matrix_file_path, barcodes_file_path, genes_file_path):
        mmwrite(matrix_file_path, adata.X)
        with open(barcodes_file_path, "w") as f:
            f.writelines(f"{bc}\n" for bc in adata.obs_names)
        with open(genes_file_path, "w") as f:
            f.writelines(f"{gene}\n" for gene in adata.var_names)


def generate_pdf_and_save_anndata(path_list, pdf_output_path, matrixmarket_output_dir):
    with PdfPages(pdf_output_path) as pdf:
        for path in path_list:
            processor = SingleCellProcessor(path)
            a_data, a_data_filtered = processor.process_data()

            #pdf.savefig(processor.plot_qc_metrics(a_data, f"{processor.sample_name} - Original"), bbox_inches='tight')
            #pdf.savefig(processor.plot_qc_metrics(a_data_filtered, f"{processor.sample_name} - Filtered"), bbox_inches='tight')

            sample_out_dir = os.path.join(matrixmarket_output_dir, processor.sample_name)
            os.makedirs(sample_out_dir, exist_ok=True)

            processor.save_counts_matrix_mm(a_data,
                                            os.path.join(sample_out_dir, f"{processor.sample_name}_original.mtx"),
                                            os.path.join(sample_out_dir, f"{processor.sample_name}_original_barcodes.tsv"),
                                            os.path.join(sample_out_dir, f"{processor.sample_name}_original_genes.tsv"))

            processor.save_counts_matrix_mm(a_data_filtered,
                                            os.path.join(sample_out_dir, f"{processor.sample_name}_filtered.mtx"),
                                            os.path.join(sample_out_dir, f"{processor.sample_name}_filtered_barcodes.tsv"),
                                            os.path.join(sample_out_dir, f"{processor.sample_name}_filtered_genes.tsv"))


our_processor = OurGroupProcessor()
our_cb_list = our_processor.generate_cellbender_output_list(our_processor.df)
our_group_QC_pdf = os.path.join(dir_out,"11_cellbender_QC_sparse_matrix", "our_group_QC_metrics.pdf")
our_group_spam_dir = os.path.join(dir_out,"11_cellbender_QC_sparse_matrix","our_group")
generate_pdf_and_save_anndata(our_cb_list,our_group_QC_pdf,our_group_spam_dir)
