## This script performs initial quality control steps. These are based on sc best practices by Fabian Theis
## Ref: https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html#
## 1- minimum genes seen in the cell
## 2- minimum cells with genes
## 3- mitochondria contaminations
## 4- doublet cell removal
## 5- Normalisation

import numpy as np
import scanpy as sc
import seaborn as sns
from scipy.stats import median_abs_deviation
import matplotlib.pyplot as plt
import os
import json
from typing import Literal
from config import _seed, RIBO_GENESET_PATH
from scipy.sparse import csr_matrix
import functools
import operator

class scRNAPreProcessor:
    def __init__(self,
        output_h5,
        output_dir,
        anndata_file = None,
        anndata:sc.AnnData=None, 
        geneset_ids:set=None, 
        mito_genes_ids:set=None, 
        mad_outlier_threshold=5, 
        mito_percentage_threshold=25, 
        genes_exp_in_min_cells_threshold=10, 
        cells_with_min_genes_threshold=100,
        gene_ids_of_interest=None,
        normalise_counts = True,
        filter_doublets = False,
        n_highly_variable_genes=5000,
        high_variable_gene_flavor="seurat_v3_paper",
        n_top_expr_genes=50,

        ):
        self.anndata = anndata if anndata else self.read_ann_data(anndata_file)
        self.geneset_ids = geneset_ids
        self.mito_genes_ids = mito_genes_ids
        self.mad_outlier_threshold = mad_outlier_threshold
        self.mito_percentage_threshold = mito_percentage_threshold
        self.genes_exp_in_min_cells_threshold = genes_exp_in_min_cells_threshold
        self.cells_with_min_genes_threshold = cells_with_min_genes_threshold
        self.gene_ids_of_interest = gene_ids_of_interest
        self.normalise_counts = normalise_counts
        self.drop_doublets = filter_doublets
        self.n_highly_variable_genes = n_highly_variable_genes
        self.high_variable_gene_flavor = high_variable_gene_flavor
        self.n_top_expr_genes = n_top_expr_genes

        self.output_dir = os.path.abspath(output_dir)
        self.output_h5 = output_h5
        self.initial_gene_counts = None
        self.initial_cell_counts = None
        self.mito_gene_percentage = None
        self.anndata_mito_genes = None
        self.all_genes_percentage = None
        self.ribosome_geneset = None
        self.statistics_data = {}
        self._seed = _seed
    
    def process(self):
        self.process_index()
        self.anndata.raw = self.anndata.copy()
        self.preliminary_counts()
        self.qc_process()
        self.plot_counts("initial")
        self.process_outliers()
        self.preliminary_filtering()
        self.populate_stats_data()
        self.process_doublets(filter_doublets=self.drop_doublets)
        self.plot_counts("filtered")
        if self.normalise_counts:
            self.process_normalisation()
        self.process_high_expression_genes()
        self.process_high_variable_genes()
        #self.get_soupX_groups()
        #self.save_data()
        return self.anndata

    def process_index(self):
        self.anndata.var = self.anndata.var.reset_index().rename({'index' : 'gene_symbols'}, axis=1).set_index('gene_ids')
    
    def set_index(self, index_name:str):
        self.anndata.var = self.anndata.var.reset_index().set_index(index_name)
    
    def process_normalisation(self, normal_type:Literal["shifted, pearson"]='shifted'):
        if normal_type == "shifted":
            layer='log1p_norm'
            self.delta_normatlisation(layer)
        elif normal_type == "pearson":
            layer='analytic_pearson_residuals'
            self.pearson_residual_normalisation(layer)
        else:
            raise ValueError("Normalisation algorithm can either be shifted or pearson")
        
        self.plot_normalised_counts(layer=layer, title=normal_type)
    
    def plot_normalised_counts(self, layer , title=None, kde=False):
        title = title if title else layer

        fig, axes = plt.subplots(1, 2, figsize=(10, 5))
        p1 = sns.histplot(self.anndata.obs["total_counts"], bins=100, kde=kde, ax=axes[0])
        axes[0].set_title("Total counts")
        p2 = sns.histplot(self.anndata.layers[layer].sum(1), bins=100, kde=kde, ax=axes[1])
        axes[1].set_title(title)
        fig.savefig(os.path.join(self.output_dir,f"{title}_normalied_counts"))
    
    def preliminary_counts(self):
        self.initial_cell_counts = self.anndata.n_obs
        self.initial_gene_counts = self.anndata.n_vars
        self.initial_mito_genes = self.get_mito_genes().copy()
        self.initial_ribo_genes = self.get_ribo_genes().copy()

    def preliminary_filtering(self):
        sc.pp.filter_genes(self.anndata, min_cells=self.genes_exp_in_min_cells_threshold, inplace=True)
        sc.pp.filter_cells(self.anndata, min_genes=self.cells_with_min_genes_threshold, inplace=True)
    
    def check_index(self):
        ## checks if index is unique or not
        return self.anndata.var.index.is_unique

    def get_mito_genes(self):
         return self.anndata.var['gene_symbols'].str.startswith("MT-")
    
    def get_ribo_genes(self):
        if self.ribosome_geneset is None:
            self.ribosome_geneset = self._get_ribosome_geneset()
        return self.anndata.var['gene_symbols'].isin(self.ribosome_geneset)

    def qc_process(self):
        ## MITO analysis
        if self.mito_genes_ids:
            not_expressed_mito_genes = set(self.get_mito_genes()).difference(self.mito_genes_ids)
            self.mito_gene_percentage = 1 - not_expressed_mito_genes/len(self.mito_genes_ids)
        
        self.anndata.var['mt'] = self.get_mito_genes().copy()
        ### RIBOSOME
        self.anndata.var["ribo"] = self.get_ribo_genes().copy()
        ##TODO: include gene_ids_of_interest to qc_metrics

        self.calculate_qc_metrics(['mt', 'ribo'])
    
    def calculate_qc_metrics(self, varnames:list, percent_top=20):
        sc.pp.calculate_qc_metrics(
            self.anndata, qc_vars=varnames, inplace=True, percent_top=[percent_top], log1p=True
        )
    
    def is_outlier(self, metric: str,):
        ## Taken from https://www.sc-best-practices.org/preprocessing_visualization/quality_control.html
        M = self.anndata.obs[metric]
        outlier = (M < np.median(M) - self.mad_outlier_threshold * median_abs_deviation(M)) | (
            np.median(M) + self.mad_outlier_threshold * median_abs_deviation(M) < M
        )
        return outlier
    
    def process_outliers(self):
        self.anndata.obs["total_counts_outlier"] = (
            self.is_outlier("log1p_total_counts") | self.is_outlier("log1p_n_genes_by_counts") | self.is_outlier("pct_counts_in_top_20_genes")
        )
        self.anndata.obs["mito_counts_outlier"] = (
            self.is_outlier("pct_counts_mt") | self.anndata.obs["pct_counts_mt"] > self.mito_percentage_threshold
            )
        ##in place filtering
        self.anndata = self.anndata[(~self.anndata.obs["total_counts_outlier"]) & (~self.anndata.obs["mito_counts_outlier"])]


    def populate_stats_data(self):
        self.statistics_data.update({
            "initial_gene_counts" : self.initial_gene_counts,
            "initial_cell_counts" : self.initial_cell_counts,
            "initial_mito_gene_counts" : len(self.anndata.var[self.initial_mito_genes]),
            "initial_mito_genes_perc" : self.mito_gene_percentage,
            "initial_genes_perc" : self.initial_cell_counts/self.geneset_ids if self.geneset_ids else None,
            "initial_ribo_gene_counts" : len(self.anndata.var[self.initial_ribo_genes]),
            "post_gene_counts" : self.anndata.n_vars,
            "post_cell_counts": self.anndata.n_obs,
            "post_mito_gene_counts": len(self.anndata.var[self.anndata.var['mt']]),
            "post_mito_genes_perc": 1 - set(self.anndata.var['mt']).difference(self.mito_genes_ids)/len(self.anndata.var['mt']) if self.mito_genes_ids else None,
            "post_genes_perc": self.anndata.n_vars/self.geneset_ids if self.geneset_ids else None,
            "post_ribo_gene_counts": len(self.anndata.var[self.anndata.var['ribo']]),
        })
    
    
    def plot_counts(self, save_name):
        sc.pl.violin(
            self.anndata, 
            ['n_genes_by_counts', 'total_counts', 'pct_counts_mt', 'pct_counts_ribo'],
            show=False,
            save=f"_{save_name}.png",
            jitter=0.4, 
            multi_panel=True
            )
        os.rename(f"./figures/violin_{save_name}.png", os.path.join(self.output_dir, f"violin_{save_name}.png"))
        
    
    def delta_normatlisation(self, layer="log1p_norm", target_sum=1e4):
        ## shifted logarithm
        ## https://www.sc-best-practices.org/preprocessing_visualization/normalization.html#shifted-logarithm
        scales_counts = sc.pp.normalize_total(
            adata=self.anndata,
            target_sum=target_sum,
            inplace=False
        )
        #log1p transform
        self.anndata.layers[layer] = sc.pp.log1p(scales_counts["X"], copy=True)
        return self.anndata
    
    def pearson_residual_normalisation(self, layer="analytic_pearson_residuals"):
        ### https://www.sc-best-practices.org/preprocessing_visualization/normalization.html#analytic-pearson-residuals
        analytic_pearson = sc.experimental.pp.normalize_pearson_residuals(self.anndata, inplace=False)
        self.anndata.layers[layer] = csr_matrix(analytic_pearson["X"])
        return self.anndata
    
    def process_doublets(self, filter_doublets=True):
        self.estimate_doublets(simulate=True, inplace=True)
        doublets_fltr = self.anndata.obs['predicted_doublet']
        self.statistics_data.update({
            'n_predicted_doublets': len(self.anndata.obs[doublets_fltr])
        })
        if filter_doublets:
            self.anndata = self.anndata[~doublets_fltr, :]

    def estimate_doublets(self, simulate=True, inplace=False):
        sim_data = sc.pp.scrublet_simulate_doublets(self.anndata,random_seed=self._seed) if simulate else None
        return sc.pp.scrublet(self.anndata, adata_sim=sim_data, random_state=self._seed, copy=bool(not inplace))
    
    def process_high_expression_genes(self, no_ribo=True, no_mito=True):
        fltr_lst = []
        if no_mito:
            fltr_lst.append(self.anndata.var["mt"])
        if no_ribo:
            fltr_lst.append(self.anndata.var["ribo"])
        
        if fltr_lst:
            combined_filter = functools.reduce(operator.or_, fltr_lst)
            not_ribo_not_mito = self.anndata.var[~combined_filter]
            data_2be_plotted = self.anndata[:,not_ribo_not_mito.index]
        else:
            data_2be_plotted = self.anndata
        
        suffix_out = self.output_h5.split(".")[0]
        sc.pl.highest_expr_genes(data_2be_plotted, gene_symbols='gene_symbols', n_top=self.n_top_expr_genes, log=True, show=False, save=f"_{suffix_out}.png")
        os.rename(f"./figures/highest_expr_genes_{suffix_out}.png", os.path.join(self.output_dir, f"high_expr_genes_{suffix_out}.png"))

    
    def process_high_variable_genes(self):
        self.get_highly_variable_genes(n_genes=self.n_highly_variable_genes, flavor=self.high_variable_gene_flavor, inplace=True)
        self.plot_highly_variable_genes()

    def plot_highly_variable_genes(self, save_name=None, dispersion_data=None):
        suffix_out = self.output_h5.split(".")[0] if save_name is None else save_name
        data2plot = dispersion_data if dispersion_data else self.anndata.var
        sc.pl.highly_variable_genes(data2plot, show=False, save=f"_{suffix_out}.png")
        os.rename(f"./figures/filter_genes_dispersion_{suffix_out}.png", os.path.join(self.output_dir, f"higly_var_genes_{suffix_out}.png"))


    def get_highly_variable_genes(self, n_genes=5000, layer="log1p_norm", flavor='seurat_v3_paper', inplace=False):
        return sc.pp.highly_variable_genes(self.anndata, n_top_genes=n_genes, layer=layer, flavor=flavor, inplace=inplace)

    def save_data(self, reset_index=True):
        with open(os.path.join(self.output_dir, "stats.json"), "w") as fh:
            json.dump(self.statistics_data, fh, indent=4)

        ### Reset the index
        if reset_index:
            self.set_index('gene_symbols') ## Since many tools accept gene names instead of gene ids
        
        self.anndata.write(os.path.join(self.output_dir, self.output_h5))
    
    @staticmethod
    def read_ann_data(file_path):
        return sc.read_10x_h5(filename=file_path)
    
    @staticmethod
    def _get_ribosome_geneset():
        with open(RIBO_GENESET_PATH) as fh:
            return {ribo_gene.strip() for ribo_gene in fh if not ribo_gene.startswith(("#", " "))}
        


class scRNAQCProcessor:
    ...