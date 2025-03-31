import os
import json
from qc_preprocess import scRNAPreProcessor


def read_gene_ids(file_path):
    with open(file_path, 'r') as f:
        return {gene_id.strip() for gene_id in f}

def execute(args):

    if not os.path.exists(args.output_dir):
        os.makedirs(args.output_dir)

    preprocessor = scRNAPreProcessor(
        anndata_file=args.input_counts_file,
        output_h5=args.output_counts_file,
        output_dir=args.output_dir,
        geneset_ids=read_gene_ids(args.geneset_ids_file) if args.geneset_ids_file else None,
        mad_outlier_threshold=args.mad_outlier_threshold,
        mito_percentage_threshold=args.mito_percentage_threshold,
        genes_exp_in_min_cells_threshold=args.genes_exp_in_min_cells_threshold,
        cells_with_min_genes_threshold=args.cells_with_min_genes_threshold,
        gene_ids_of_interest=read_gene_ids(args.gene_ids_of_interest_file) if args.gene_ids_of_interest_file else None,
        filter_doublets=args.drop_doublets,
        n_highly_variable_genes=args.n_highly_variable_genes,
        high_variable_gene_flavor=args.dispersion_flavor,
        n_top_expr_genes=args.n_top_expr_genes,
    )
    preprocessor.process()
    preprocessor.save_data()

    with open(os.path.join(args.output_dir, "parameters.json"), 'w') as fh:
        json.dump(vars(args), fh, indent=4)
    

if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description='Initites pre-process (pre-filtering and normalisation) step for scRNA data')
    
    parser.add_argument('--input_counts_file', required=True, help='scRNA Raw counts input')
    parser.add_argument('--output_counts_file', required=True, help='scRNA filtered counts output')
    parser.add_argument('--output_dir', required=False, help='Output folder for plots and stats', default=os.getcwd())
    parser.add_argument('--geneset_ids_file', required=False, help='All the geneset')
    parser.add_argument('--mad_outlier_threshold', required=False, default=5, type=int, help="Median Absolute Deviations threshold to filter outliers")
    parser.add_argument('--mito_percentage_threshold', required=False, default=8, type=int, help="Mito percentage upper limit to filter mitochondria contamination")   
    parser.add_argument('--genes_exp_in_min_cells_threshold', required=False, default=0, type=int, help="Drops genes if not seen in at least ... number of cells")
    parser.add_argument('--cells_with_min_genes_threshold', required=False, default=0, type=int, help="Drops cells if have less than ... number of genes")
    parser.add_argument('--n_top_expr_genes', required=False, default=50, type=int, help="Plots number of most expressed genes (ribo and mito genes not included)")
    parser.add_argument('--n_highly_variable_genes', required=False, default=5000, type=int, help="Rnaks n number of highly variable genes")
    parser.add_argument('--dispersion_flavor', required=False, default='seurat_v3_paper', type=str, choices=['seurat_v3_paper', 'seurat_v3', 'seurat', 'cell_ranger'], help="Dispersion flavour defaults 'seurat_v3_paper'")
    parser.add_argument('--gene_ids_of_interest_file', required=False, help="populates additional quality metrics to the genes of interest") 
    parser.add_argument('--normalisation_method', required=False, default='shifted', type=str, help="Normalise the counts either by 'shifted' or 'pearson'")
    parser.add_argument('--drop_doublets', required=False, action='store_true', help="Filters out the predicted doublets inplace")

    args = parser.parse_args()
    execute(args)