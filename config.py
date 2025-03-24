import scanpy as sc
import os

### Plotting settings
sc.settings.verbosity = 0
sc.settings.set_figure_params(
    dpi=300,
    facecolor="white",
    frameon=False,
)

### Seed
_seed = 42

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")

RIBO_GENESET_PATH = os.path.join(DATA_DIR, "KEGG_RIBOSOME_GENESET.txt")

