import numpy as np
import scanpy as sc
import seaborn as sns
import anndata as ad
import matplotlib.pyplot as plt
from scipy.sparse import csr_matrix
from typing import Literal, List
import os
from qc_preprocess import scRNAQCProcessor


class scRNAseqData:
    def __init__(self, input_file:str, output_dir:str, **kwargs):
        self.input_file = input_file
        self.output_dir = output_dir
        self.get_anndata()
        self.stats_data = {}
        self.qc_processor = scRNAQCProcessor(
            anndata=self.anndata,
            output_dir=self.output_dir
            **kwargs
        )

    def get_anndata(self):
        self.anndata = ad.read_h5ad(self.input_file)
        return self.anndata
    
    def process_index(self):
        self.anndata.var = self.anndata.var.reset_index().rename({'index' : 'gene_symbols'}, axis=1).set_index('gene_ids')
    
    def set_index(self, index_name:str):
        self.anndata.var = self.anndata.var.reset_index().set_index(index_name)


class scRNAseqProcess(scRNAseqData):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)


    def plot_projections(self):
        ## 1- plot dispersion
        ## 2- plot pca
        ## 3- plot tsne
        ## 4- plot umap
        ...

    def process_initial_filtering(self):
        self.qc_processor.process()

    def process_normalisation(self):
        ...

    def process_cell_annotation(self):
        ...

    def process_gene_regulatory_networks(self):
        ...

    def process_differential_expression(self):
        ...

    def process_clustering(self):
        ...

    def process_trajectory(self):
        ...
    
    def get_raw_counts(self):
       self.anndata.raw = self.anndata.X
    
    def get_highly_variable_genes(self, layer=None, n_top_genes=2000, flavor="seurat_v3_paper", inplace=False):
        return sc.pp.highly_variable_genes(self, layer=layer, flavor=flavor, n_top_genes=n_top_genes, inplace=inplace)
