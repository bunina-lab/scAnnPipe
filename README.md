# Single-Cell RNA Analysis Pipeline

This pipeline provides tools for preprocessing and quality control analysis of single-cell RNA sequencing data. It performs initial quality filtering, generates diagnostic plots, and annotates genes of interest.

## Features

- Quality control filtering based on:
  - Mitochondrial gene percentage
  - Number of genes per cell
  - Number of cells expressing each gene
  - Outlier detection using Median Absolute Deviations (MAD)
- Doublet detection and filtering
- Highly variable gene selection
- Gene expression visualization
- Custom gene set analysis
- Comprehensive quality metrics and plots

## Installation

1. Clone the repository:
```bash
git clone https://github.com/yourusername/pipe_scanalysis.git
cd pipe_scanalysis
```

2. Create a virtual environment (recommended):
```bash
python -m venv venv
source venv/bin/activate  # On Unix/macOS
# or
.\venv\Scripts\activate  # On Windows
```

3. Install dependencies:
```bash
pip install -r requirements.txt
```

## Usage

The pipeline can be run using the `run_qc_process.py` script. Here's the basic command structure:

```bash
python run_qc_process.py --input_counts_file <input.h5> --output_counts_file <output.h5ad> [options]
```

### Required Parameters

- `--input_counts_file`: Path to the input AnnData file containing raw counts
- `--output_counts_file`: Path where the processed AnnData file will be saved

### Optional Parameters

- `--output_dir`: Output directory for plots and statistics (default: current directory)
- `--geneset_ids_file`: File containing a list of gene IDs to analyze
- `--mad_outlier_threshold`: MAD threshold for outlier filtering (default: 5)
- `--mito_percentage_threshold`: Maximum allowed mitochondrial gene percentage (default: 8)
- `--genes_exp_in_min_cells_threshold`: Minimum number of cells expressing each gene (default: 0)
- `--cells_with_min_genes_threshold`: Minimum number of genes per cell (default: 0)
- `--n_top_expr_genes`: Number of top expressed genes to plot (default: 50)
- `--n_highly_variable_genes`: Number of highly variable genes to select (default: 5000)
- `--dispersion_flavor`: Method for highly variable gene selection (default: 'seurat_v3_paper')
  - Options: 'seurat_v3_paper', 'seurat_v3', 'seurat', 'cell_ranger'
- `--gene_ids_of_interest_file`: File containing gene IDs for additional quality metrics
- `--normalisation_method`: Normalization method to use (default: 'shifted')
  - Options: 'shifted', 'pearson'
- `--drop_doublets`: Flag to filter out predicted doublets

### Example Usage

```bash
python run_qc_process.py \
    --input_counts_file unfiltered_gex_counts.h5 \
    --output_counts_file filtered_counts.h5ad \
    --output_dir results \
    --mito_percentage_threshold 10 \
    --cells_with_min_genes_threshold 200 \
    --n_highly_variable_genes 3000 \
    --drop_doublets
```

## Output

The pipeline generates:
1. A filtered and annotated AnnData file with processed data
2. Quality control plots and statistics in the specified output directory

## Dependencies

The pipeline requires Python 3.x and the following main packages:
- scanpy
- anndata
- numpy
- pandas
- matplotlib
- scikit-learn
- umap-learn

For a complete list of dependencies and their versions, see `requirements.txt`.

## License

See [MIT License](https://opensource.org/license/mit)

## Contributing

[Add contribution guidelines here] 