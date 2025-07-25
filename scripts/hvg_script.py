# %%
import logging
import sys

# Configure logging to output to both stdout and stderr
# Create a custom formatter
formatter = logging.Formatter("%(asctime)s - %(levelname)s - %(message)s")

# Create stdout handler for INFO and DEBUG messages
stdout_handler = logging.StreamHandler(sys.stdout)
stdout_handler.setLevel(logging.INFO)
stdout_handler.setFormatter(formatter)

# Create stderr handler for WARNING, ERROR, and CRITICAL messages
stderr_handler = logging.StreamHandler(sys.stderr)
stderr_handler.setLevel(logging.WARNING)
stderr_handler.setFormatter(formatter)

# Get the root logger and configure it
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.addHandler(stdout_handler)
logger.addHandler(stderr_handler)

logging.info("Starting the HVG script...")
logging.info("Importing necessary libraries...")

# Test the different log levels to show where they go:
logging.info("This INFO message goes to stdout")
logging.warning("This WARNING message goes to stderr")
logging.error("This ERROR message goes to stderr")

import gc
import os
import sys
import warnings
import copy
from pathlib import Path

# %%
import anndata as ad
import scanpy as sc
import scipy
import pickle
import data.dataset_keys
import numpy as np
import pandas as pd
import src.preprocessing
import src.utils
from src.harmonize_features import HarmonizeFeatures
import tqdm
from config import config as config

# %%
sc.settings.verbosity = 3

# import matplotlib.pyplot as plt
# import seaborn as sns
# plt.rcParams["font.family"] = "DeJavu Serif"
# plt.rcParams["font.serif"] = ["Times New Roman"]

# %%
gc.enable()

# %%
logging.info("Reading the anndata file...")
# adata_file_path = "/lustre/scratch127/cellgen/cellgeni/tickets/tic-4101/data/subset_100k.h5ad"
adata_file_path = "/lustre/scratch126/cellgen/haniffa/ar32/combine_IA_WE_with_fix/raw_outer_IA_WE_data_20250704.h5ad"
logging.info(f"Loading data from: {adata_file_path}")
adata = sc.read_h5ad(adata_file_path)
logging.info(f"Loaded data with shape: {adata.shape}")
logging.info(f"Data overview: {adata}")

# %%
logging.info("=== STEP 1: Loading gene intersection data ===")
genes_intersect = pd.read_csv(
    "/lustre/scratch126/cellgen/haniffa/ar32/combine_IA_WE_with_fix/gene_var_cut_new_20250704.csv"
)
logging.info(f"Loaded {len(genes_intersect)} genes for intersection")
adata = adata[:, genes_intersect.iloc[:, 0]].copy()
logging.info(f"After gene intersection filtering: {adata}")

# %%
logging.info("=== STEP 2: Removing specific gene types ===")
exclude_genes = {
    "cell_cycle": True,
    "bcr": True,
    "tcr": True,
    "mt": False,
    "ribo": False,
}
logging.info(f"Excluding specific genes: {exclude_genes}")
logging.info(f"Initial anndata shape: {adata.shape}")
logging.info("Starting gene filtering...")
adata, removed_genes_as_specific_genes = src.preprocessing.remove_specific_genes(
    adata=adata,
    exclude_cell_cycle=exclude_genes["cell_cycle"],
    exclude_bcr=exclude_genes["bcr"],
    exclude_tcr=exclude_genes["tcr"],
    exclude_mt=exclude_genes["mt"],
    exclude_ribo=exclude_genes["ribo"],
)
gc.collect()
logging.info(f"After gene filtering, shape: {adata.shape}")
logging.info(f"Resulting anndata shape: {adata.shape}")

# %% [markdown]
# ## Filter genes based on min count

# %%
# logging.info("=== STEP 3: Cell filtering based on minimum gene counts ===")
# logging.info("Filtering genes based on minimum counts...")
# _, scanpy_filter_cells_min_genes = sc.pp.filter_cells(adata, min_genes=1, inplace=False)
# logging.info(f"Filter cells min genes result: {scanpy_filter_cells_min_genes}")

# # %%
# scanpy_filter_cells_min_genes_threshold = 200
# scanpy_filter_cells_min_genes_bool = (
#     scanpy_filter_cells_min_genes < scanpy_filter_cells_min_genes_threshold
# )
# adata = adata[~scanpy_filter_cells_min_genes_bool, :].copy()
# logging.info(f"Shape after filtering: {adata.shape}, is_view: {adata.is_view}")
# gc.collect()

# # %%
# logging.info("=== STEP 4: Cell filtering based on minimum total counts ===")
# _, scanpy_filter_cells_min_counts = sc.pp.filter_cells(
#     adata, min_counts=1, inplace=False
# )

# logging.info(f"Cells with counts < 400: {np.sum(scanpy_filter_cells_min_counts < 400)}")
# cell_counts = adata[scanpy_filter_cells_min_counts < 400, :].obs.value_counts(
#     ["handle_anndata", "cell_type_lvl5"]
# )
# logging.info(f"Cell type counts for low count cells: {cell_counts}")

# %% [markdown]
# ## No filtering ???????!!!!!

# %% [markdown]
# the fuck is happening here

# %%
logging.info("=== STEP 5: Gene filtering based on minimum counts across cells ===")
scanpy_filter_genes_threshold = 5
logging.info(f"Using gene filter threshold: {scanpy_filter_genes_threshold}")
cells_subset, number_per_cell = sc.pp.filter_genes(
    adata, min_counts=scanpy_filter_genes_threshold, inplace=False
)
# Most of them are actually `Braun dataset` only genes.
gene_stats = (
    (
        adata.var[number_per_cell < scanpy_filter_genes_threshold][["hgnc_Braun-0"]]
        != "NA"
    )
    .sum(axis=1)
    .value_counts()
)
logging.info(f"Gene statistics for filtered genes: {gene_stats}")

# %%
removed_genes_as_filtered = set(adata.var.index) - set(adata.var.index[cells_subset])
adata = adata[:, cells_subset].copy()
logging.info(f"Shape after gene filtering: {adata.shape}, is_view: {adata.is_view}")
gc.collect()

# %% [markdown]
# ## Calculate HVG and Define min and max mean values

# %%
logging.info("=== STEP 6: Setting up batch definitions for HVG calculation ===")
logging.info("Creating library batch definitions...")
batch_definitions = dict(
    concatenated_batch_key_columns_for_hvg_calculation=[
        "integration_donor",
        "integration_biological_unit",
        "integration_sample_status",
        "integration_library_platform_coarse",
    ]
)

batch_key = "concatenated_integration_covariates"
logging.info(f"Creating batch key: {batch_key}")
src.preprocessing.library_batch_generator(
    adata=adata,
    batch="concatenated_batch_key_columns_for_hvg_calculation",
    definitions=batch_definitions,
    new_col_name=batch_key,
)
logging.info("Batch key generation completed")

# %%
logging.info("=== STEP 7: Creating batch statistics and intersection copy ===")
adata_obs_with_new_column = adata.obs.copy()

df_concatenated_batch = (
    adata.obs.copy()
    .groupby([batch_key, "handle_anndata"])
    .size()
    .reset_index(name="Counts")
)
df_concatenated_batch = df_concatenated_batch[df_concatenated_batch["Counts"] != 0]
df_concatenated_batch_value_counts = df_concatenated_batch[
    "handle_anndata"
].value_counts()
logging.info(f"Created batch statistics with {len(df_concatenated_batch)} entries")

# genes_intersect = pd.read_csv("/lustre/scratch126/cellgen/haniffa/ar32/combine_IA_WE_with_fix/gene_var_cut_new_20250704.csv")
# adata_intersection = adata[:,genes_intersect.iloc[:,0]]
adata_intersection = adata.copy()
logging.info(f"Created intersection copy: {adata_intersection}")

# %%
logging.info("=== STEP 8: Normalizing and log-transforming data ===")
# as it is for each data row separately, batch key is not needed.
logging.info("Normalizing total counts to 1e6...")
sc.pp.normalize_total(adata_intersection, target_sum=1e6)
logging.info("Applying log1p transformation...")
sc.pp.log1p(adata_intersection)

# %%
logging.info(
    "=== STEP 9: Computing highly variable genes (global, no batch correction) ==="
)
# Calculate without batch keys for testing
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    logging.info("Computing HVG without batch key...")
    global_hvg_adata_var = src.preprocessing.highly_variable_genes(
        adata=adata_intersection,
        layer=None,
        batch_key=None,
        subset=False,
        inplace=False,
    )
global_hvg_adata_var.index = adata_intersection.var.index.copy()
logging.info("Global HVG calculation completed")

# adata_intersection = adata_intersection[adata_intersection.obs[batch_key] != 'HT-231_Yu_et_al_Cell_Fresh_3GEX'].copy()
gc.collect()

# %%
logging.info("=== STEP 10: Computing HVG with batch correction ===")
# Calculate dispersion and mean with scanpy function to determine max and min means.
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    logging.info("Computing HVG with batch key for dispersion and mean calculation...")
    src.preprocessing.highly_variable_genes(
        adata=adata_intersection,
        layer=None,
        batch_key=batch_key,
        subset=False,
        inplace=True,
    )
    adata_intersection.var.drop(columns=["highly_variable"], inplace=True)

adata_var_before = adata_intersection.var.copy()
logging.info("Batch-corrected HVG calculation completed")
gc.collect()

# %%
logging.info("=== STEP 11: HVG selection with thresholds ===")
logging.info(
    f"Highly variable intersection: `{adata_intersection.var['highly_variable_intersection'].sum()}`"
)
logging.info(
    f"Number of batch defined: `{len(adata_intersection.obs[batch_key].unique())}`"
)

hvg_dict = dict(hvg_number=2**14 - 2**12, min_mean=0.05, max_mean=6)
logging.info(f"HVG parameters: {hvg_dict}")

# %%
logging.info("=== STEP 12: Selecting final HVG set ===")
calculated_top_genes = []
logging.info("Running HVG selection algorithm...")
passed_thresholding, final_hvg_selection, adata_var_sorted_index = (
    src.preprocessing.select_hvgs(
        adata_var=adata_intersection.var.copy(),
        top_gene_number=hvg_dict["hvg_number"],
        min_mean=hvg_dict["min_mean"],
        max_mean=hvg_dict["max_mean"],
    )
)
removed_genes_as_mean_threshold = set(adata_var_before.index) - passed_thresholding
genes_not_selected_as_hvg = passed_thresholding - final_hvg_selection
logging.info(f"Passed mean thresholding: {len(passed_thresholding)} genes")
logging.info(f"Final HVG selection: {len(final_hvg_selection)} genes")
calculated_top_genes.append(
    ["main_hvg", hvg_dict["hvg_number"], "rank", adata_var_sorted_index.copy()]
)
calculated_top_genes.append(
    ["main_hvg", hvg_dict["hvg_number"], "hvg", final_hvg_selection.copy()]
)

logging.info(f"Resulting anndata shape: {adata_intersection.shape}")

final_adata_var = adata_intersection.var.copy()

# %%
logging.info("=== STEP 13: Cleaning up and reloading original data ===")
# Keep the lists, sets, dataframes, and dictionaries to report back on the original anndata object.
# This also aims to remove log normalized values, as their sole purpose was HVG calculation
logging.info("Deleting processed data copies to free memory...")
del adata_intersection, adata
gc.collect()
logging.info("Reloading original data...")
adata = ad.read_h5ad(adata_file_path)
adata = adata[adata_obs_with_new_column.index, :].copy()
logging.info(f"Reloaded and filtered data: {adata}")

# %%
logging.info("=== STEP 14: Creating HVG annotation columns ===")
df = adata.var.copy()
columns_to_move = [i for i in final_adata_var.columns if i not in adata.var.columns]
df = pd.merge(
    df, final_adata_var[columns_to_move], left_index=True, right_index=True, how="left"
)
logging.info(f"Merged {len(columns_to_move)} new columns into gene metadata")

df_c1 = "highly_variable_omitted_reason"
logging.info("Creating highly_variable_omitted_reason column...")

df[df_c1] = src.utils.NA_CELL_TYPE_PLACEHOLDER
for specific_genes in removed_genes_as_specific_genes.keys():
    specific_genes_index = df.index.isin(
        removed_genes_as_specific_genes[specific_genes]
    )
    df.loc[specific_genes_index, df_c1] = specific_genes
df.loc[df.index.isin(removed_genes_as_filtered), df_c1] = "scanpy_filter_genes"

# df.loc[adata.var[(adata.var["intersection"] != 1) & (df[df_c1] == "NA")].index, df_c1] = "not_intersection"
df[df_c1] = df[df_c1].astype("category")
omitted_counts = df[df_c1].value_counts()
logging.info(f"Highly variable omitted reason counts: {omitted_counts}")

df_c2 = "highly_variable_reason"
logging.info("Creating highly_variable_reason column...")

df[df_c2] = "highly_variable_omitted"
df.loc[df.index.isin(removed_genes_as_mean_threshold), df_c2] = "thresholded_mean_range"
df.loc[df.index.isin(genes_not_selected_as_hvg), df_c2] = "thresholded_top_genes"
df.loc[df.index.isin(final_hvg_selection), df_c2] = "hvg"
df[df_c2] = df[df_c2].astype("category")

reason_counts = df[df_c2].value_counts()
logging.info(f"Highly variable reason counts: {reason_counts}")

df_c3 = "highly_variable"
logging.info("Creating final highly_variable boolean column...")
df[df_c3] = df[df_c2] == "hvg"
hvg_counts = df[df_c3].value_counts()
logging.info(f"Highly variable counts: {hvg_counts}")

hvg_report = {
    "exclude_genes": exclude_genes,
    "scanpy_filter_genes_threshold": scanpy_filter_genes_threshold,
    "highly_variable_parameters": hvg_dict,
}

# %%
logging.info("=== STEP 15: Creating gene ranking columns ===")
logging.info("Computing gene rankings...")
v1 = np.array(
    [
        (
            np.where(adata_var_sorted_index == i)[0][0]
            if i in adata_var_sorted_index
            else -1
        )
        for ind, i in enumerate(adata.var.index)
    ]
)
v2 = np.array(
    [
        (
            np.where(adata_var_sorted_index == i)[0][0]
            if i in adata_var_sorted_index
            else -1
        )
        for ind, i in enumerate(adata.var.index)
    ]
)
df_c4 = "highly_variable_ranking"
df_c5 = "highly_variable_ranking_global"
df[df_c4] = v1
df[df_c5] = v2
df[df_c4] = df[df_c4].astype(int)
df[df_c5] = df[df_c5].astype(int)
logging.info("Gene ranking columns created")

# %%
logging.info("=== STEP 16: Finalizing gene metadata ===")
logging.info("Processing highly_variable_intersection column...")
df["highly_variable_intersection"] = df["highly_variable_intersection"].astype(str)
df["highly_variable_intersection"].replace(
    src.utils.MISSING_VALUES, src.utils.NA_CELL_TYPE_PLACEHOLDER, inplace=True
)
df["highly_variable_intersection"] = df["highly_variable_intersection"].astype(
    "category"
)
df["highly_variable"] = df["highly_variable"].astype(int)
logging.info("Gene metadata finalization completed")

# %%
logging.info("=== STEP 17: Applying HVG selection and preparing for export ===")
adata.var = df.copy()
adata.obs = adata_obs_with_new_column.copy()
adata.uns["highly_variable"] = hvg_report

initial_shape = adata.shape
adata = adata[:, adata.var["highly_variable"] == 1].copy()
logging.info(f"Applied HVG filter: {initial_shape} -> {adata.shape}")

# %%
logging.info("=== STEP 18: Saving results ===")
split_path = os.path.splitext(adata_file_path)
# write_name = "results/subset_100k_processed.h5ad"
# write_name_pickle = "results/subset_100k_processed.pickle"
write_name = "results/raw_outer_IA_WE_data_20250704.h5ad"
write_name_pickle = "results/raw_outer_IA_WE_data_20250704.pickle"
# %%
logging.info(f"Saving pickle file to: {write_name_pickle}")
with open(write_name_pickle, "wb") as file:
    pickle.dump(calculated_top_genes, file)
logging.info(f"Saving processed data to: {write_name}")
adata.write_h5ad(write_name)
logging.info("File saving completed")

# %%
logging.info("=== ANALYSIS COMPLETED ===")
logging.info(f"Calculated top genes: {calculated_top_genes}")
logging.info("HVG script completed successfully!")
logging.info("FINISHED")
