#!/usr/bin/env python3
"""
Combined preprocessing and downsampling script.

This script performs normalization, HVG detection, PCA calculation,
and then downsampling on single-cell RNA-seq data stored in AnnData format.
It saves both the PCA matrix and downsampled cell indices.
"""

import argparse
import json
import logging
import sys

import anndata as ad
import scanpy as sc
import numpy as np
from celltypist.samples import downsample_adata


def setup_logging(log_level: str = "INFO", log_file: str = "preprocess_downsample.log"):
    """Set up logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler(log_file),
        ],
    )
    return logging.getLogger(__name__)


def load_adata(file_path: str) -> ad.AnnData:
    """Load AnnData object from file."""
    logger = logging.getLogger(__name__)
    logger.info(f"Loading AnnData object from: {file_path}")

    try:
        adata = ad.read_h5ad(file_path)
        logger.info(f"Successfully loaded AnnData with shape: {adata.shape}")
        logger.info(f"Available obs columns: {list(adata.obs.columns)}")
        return adata
    except Exception as e:
        logger.error(f"Failed to load AnnData: {e}")
        raise


def normalize_and_log(adata: ad.AnnData) -> ad.AnnData:
    """Normalize and log-transform the data."""
    logger = logging.getLogger(__name__)
    logger.info("Starting normalization and log transformation")

    # Make a copy to avoid modifying original
    adata_processed = adata.copy()

    # Store raw counts
    adata_processed.raw = adata_processed

    # Normalize to 10,000 reads per cell
    sc.pp.normalize_total(adata_processed, target_sum=1e4)
    logger.info("Normalized to 10,000 reads per cell")

    # Log transform
    sc.pp.log1p(adata_processed)
    logger.info("Applied log1p transformation")

    return adata_processed


def find_hvg(adata: ad.AnnData, n_top_genes: int = 2000) -> ad.AnnData:
    """Find highly variable genes."""
    logger = logging.getLogger(__name__)
    logger.info(f"Finding highly variable genes (top {n_top_genes})")

    # Find highly variable genes
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_genes, subset=False)

    n_hvg = adata.var["highly_variable"].sum()
    logger.info(f"Found {n_hvg} highly variable genes")

    # Subset to HVG
    adata_hvg = adata[:, adata.var["highly_variable"]].copy()
    logger.info(f"Subset to HVG shape: {adata_hvg.shape}")

    return adata_hvg


def calculate_pca(
    adata: ad.AnnData, n_comps: int = 50, use_highly_variable: bool = None
) -> ad.AnnData:
    """Calculate PCA."""
    logger = logging.getLogger(__name__)
    logger.info(f"Calculating PCA with {n_comps} components")

    # Calculate PCA
    sc.tl.pca(
        adata,
        n_comps=n_comps,
        svd_solver="arpack",
        use_highly_variable=use_highly_variable,
    )
    logger.info(f"Calculated PCA with {adata.obsm['X_pca'].shape[1]} components")

    return adata


def perform_downsampling_with_indices(
    adata: ad.AnnData, celltype_key: str, fraction: float, filter_nan: bool = True
) -> tuple[ad.AnnData, list, list]:
    """Perform downsampling using celltypist and return indices."""
    logger = logging.getLogger(__name__)

    # Store reference to original AnnData for index mapping
    original_n_obs = adata.n_obs

    # Check if celltype_key exists
    if celltype_key not in adata.obs.columns:
        available_keys = list(adata.obs.columns)
        logger.error(f"Celltype key '{celltype_key}' not found in adata.obs")
        logger.error(f"Available keys: {available_keys}")
        raise ValueError(f"Celltype key '{celltype_key}' not found in adata.obs")

    # Filter out NaN cells if requested
    if filter_nan:
        logger.info("Filtering out cells with NaN values in celltype annotation")

        # Get non-NaN mask
        non_nan_mask = adata.obs[celltype_key].notna()
        n_nan_cells = (~non_nan_mask).sum()
        n_valid_cells = non_nan_mask.sum()

        logger.info(f"Total cells in original dataset: {original_n_obs}")
        logger.info(f"Cells with NaN in '{celltype_key}': {n_nan_cells}")
        logger.info(f"Valid cells for downsampling: {n_valid_cells}")

        if n_valid_cells == 0:
            raise ValueError(
                f"No valid (non-NaN) cells found in '{celltype_key}' column"
            )

        # Get indices of non-NaN cells in original dataset
        original_valid_indices = np.where(non_nan_mask)[0]

        # Create temporary AnnData view for valid cells only (no copying)
        adata_valid = adata[non_nan_mask]
        logger.info(f"Using {adata_valid.n_obs} valid cells for downsampling")
    else:
        logger.info("Downsampling all cells (including NaN values)")
        original_valid_indices = np.arange(original_n_obs)
        adata_valid = adata
        logger.info(f"Using all {adata_valid.n_obs} cells for downsampling")

    # Calculate target number of cells based on fraction
    total_valid_cells = adata_valid.n_obs
    n_cells = int(total_valid_cells * fraction)

    logger.info(
        f"Starting downsampling with celltype_key='{celltype_key}', fraction={fraction}"
    )
    logger.info(
        f"Valid cells for downsampling: {total_valid_cells}, Target cells: {n_cells}"
    )

    # Validate fraction
    if not 0 < fraction <= 1:
        raise ValueError(f"Fraction must be between 0 and 1, got {fraction}")

    # Log cell type distribution before downsampling (valid cells only)
    celltype_counts = adata_valid.obs[celltype_key].value_counts()
    logger.info("Cell type distribution before downsampling (valid cells only):")
    for celltype, count in celltype_counts.items():
        logger.info(f"  {celltype}: {count} cells")

    # Store original cell barcodes
    original_barcodes = list(adata.obs.index)

    # Perform downsampling on valid cells with return_index=True to get indices
    selected_indices_in_valid = downsample_adata(
        adata=adata_valid,
        mode="total",
        n_cells=n_cells,
        by=celltype_key,
        random_state=4,
        return_index=True,
    )

    # Map back to original dataset indices
    selected_indices_in_original = original_valid_indices[selected_indices_in_valid]

    # Convert numpy array to Python list for JSON serialization
    selected_indices_list = selected_indices_in_original.tolist()

    # Create downsampled adata using the original indices (only copy what we need)
    adata_downsampled = adata[selected_indices_in_original].copy()

    logger.info(f"Downsampling completed. New shape: {adata_downsampled.shape}")
    logger.info(
        f"Selected {len(selected_indices_in_original)} cell indices from original dataset"
    )

    # Log cell type distribution after downsampling
    celltype_counts_after = adata_downsampled.obs[celltype_key].value_counts()
    logger.info("Cell type distribution after downsampling:")
    for celltype, count in celltype_counts_after.items():
        logger.info(f"  {celltype}: {count} cells")

    # Get the selected cell barcodes from original dataset
    selected_barcodes = [original_barcodes[i] for i in selected_indices_in_original]

    return adata_downsampled, selected_indices_list, selected_barcodes


def save_pca_matrix(adata: ad.AnnData, output_path: str):
    """Save PCA matrix as numpy array."""
    logger = logging.getLogger(__name__)
    logger.info(f"Saving PCA matrix to: {output_path}")

    try:
        # Get PCA matrix
        pca_matrix = adata.obsm["X_pca"]

        # Save as numpy array
        np.save(output_path, pca_matrix)
        logger.info(f"Successfully saved PCA matrix with shape: {pca_matrix.shape}")

        # Also save cell barcodes/indices for reference
        barcodes_path = output_path.replace(".npy", "_barcodes.txt")
        with open(barcodes_path, "w", encoding="utf-8") as f:
            for barcode in adata.obs.index:
                f.write(f"{barcode}\n")
        logger.info(f"Saved cell barcodes to: {barcodes_path}")

    except Exception as e:
        logger.error(f"Failed to save PCA matrix: {e}")
        raise


def save_indices(
    indices: list, barcodes: list, output_path: str, format_type: str = "json"
):
    """Save cell indices in specified format."""
    logger = logging.getLogger(__name__)
    logger.info(f"Saving cell indices to: {output_path} in {format_type} format")

    try:
        if format_type.lower() == "json":
            # Save as JSON with both indices and barcodes
            data = {
                "selected_indices": indices,
                "selected_barcodes": barcodes,
                "total_selected": len(indices),
                "metadata": {
                    "description": "Downsampled cell indices and barcodes",
                    "format": "JSON",
                },
            }
            with open(output_path, "w", encoding="utf-8") as f:
                json.dump(data, f, indent=2)

        elif format_type.lower() == "txt":
            # Save as simple text file with barcodes
            with open(output_path, "w", encoding="utf-8") as f:
                for barcode in barcodes:
                    f.write(f"{barcode}\n")

        elif format_type.lower() == "csv":
            # Save as CSV with index and barcode columns
            import csv

            with open(output_path, "w", newline="", encoding="utf-8") as f:
                writer = csv.writer(f)
                writer.writerow(["original_index", "barcode"])
                for idx, barcode in zip(indices, barcodes):
                    writer.writerow([idx, barcode])

        else:
            raise ValueError(
                f"Unsupported format: {format_type}. Use 'json', 'txt', or 'csv'"
            )

        logger.info(f"Successfully saved {len(indices)} cell indices")

    except Exception as e:
        logger.error(f"Failed to save indices: {e}")
        raise


def save_processed_adata(adata: ad.AnnData, output_path: str):
    """Save processed AnnData object (optional)."""
    logger = logging.getLogger(__name__)
    logger.info(f"Saving processed AnnData to: {output_path}")

    try:
        adata.write_h5ad(output_path)
        logger.info(f"Successfully saved processed AnnData with shape: {adata.shape}")
    except Exception as e:
        logger.error(f"Failed to save processed AnnData: {e}")
        raise


def remove_x_matrix(adata: ad.AnnData) -> ad.AnnData:
    """Remove .X matrix from AnnData."""
    logger = logging.getLogger(__name__)
    logger.info(f"Removing .X matrix from AnnData {adata}")

    # Create a copy and remove X matrix
    adata_no_x = adata.copy()
    del adata_no_x.X
    adata_no_x.X = None

    logger.info(f"Successfully removed .X matrix {adata_no_x}")
    return adata_no_x


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Preprocess AnnData, calculate PCA, and downsample cells",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "input_file", type=str, help="Path to input AnnData file (.h5ad)"
    )

    parser.add_argument(
        "--celltype_key",
        type=str,
        required=True,
        help="Column name in adata.obs containing cell type annotations",
    )

    parser.add_argument(
        "--fraction",
        type=float,
        required=True,
        help="Fraction of the dataset to downsample to (between 0 and 1)",
    )

    parser.add_argument(
        "--lognorm",
        action="store_true",
        help="Apply normalization and log1p transformation",
    )

    parser.add_argument(
        "--hvg", action="store_true", help="Find highly variable genes before PCA"
    )

    parser.add_argument(
        "--n_hvg",
        type=int,
        default=2000,
        help="Number of highly variable genes to select",
    )

    parser.add_argument(
        "--n_pca", type=int, default=50, help="Number of PCA components to calculate"
    )

    parser.add_argument(
        "--output_pca",
        type=str,
        default="pca.npy",
        help="Output path for PCA matrix (.npy file)",
    )

    parser.add_argument(
        "--output_indices",
        type=str,
        default="downsampled_indices.json",
        help="Output file path for cell indices",
    )

    parser.add_argument(
        "--filter_nan",
        action="store_true",
        default=True,
        help="Filter out cells with NaN values in celltype annotation before downsampling (default: True)",
    )

    parser.add_argument(
        "--include_nan",
        action="store_true",
        help="Include cells with NaN values in downsampling (overrides --filter_nan)",
    )

    parser.add_argument(
        "--indices_format",
        type=str,
        default="json",
        choices=["json", "txt", "csv"],
        help="Format for saving indices (json, txt, or csv)",
    )

    parser.add_argument(
        "--output_adata",
        type=str,
        help="Optional: Save processed AnnData object to this path",
    )

    parser.add_argument(
        "--log_level",
        type=str,
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR"],
        help="Logging level",
    )

    args = parser.parse_args()

    # Setup logging
    logger = setup_logging(args.log_level)

    # Handle NaN filtering logic
    filter_nan = args.filter_nan and not args.include_nan

    logger.info("=" * 60)
    logger.info("Starting combined preprocessing and downsampling pipeline")
    logger.info("=" * 60)
    logger.info(f"Input file: {args.input_file}")
    logger.info(f"Apply log normalization: {args.lognorm}")
    logger.info(f"Use HVG: {args.hvg}")
    logger.info(f"Number of HVG: {args.n_hvg}")
    logger.info(f"Number of PCA components: {args.n_pca}")
    logger.info(f"Celltype key: {args.celltype_key}")
    logger.info(f"Downsample fraction: {args.fraction}")
    logger.info(f"Filter NaN cells: {filter_nan}")
    logger.info(f"Output PCA file: {args.output_pca}")
    logger.info(f"Output indices file: {args.output_indices}")
    logger.info(f"Indices format: {args.indices_format}")

    # Step 1: Load data
    adata = load_adata(args.input_file)

    # Step 2: Normalize and log transform (optional)
    if args.lognorm:
        adata = normalize_and_log(adata)
    else:
        logger.info("Skipping normalization and log transformation")

    # Step 3: Find HVG (optional)
    if args.hvg:
        adata = find_hvg(adata, n_top_genes=args.n_hvg)

    # Step 4: Calculate PCA
    adata = calculate_pca(adata, n_comps=args.n_pca, use_highly_variable=args.hvg)

    # Step 5: Remove .X matrix
    adata = remove_x_matrix(adata)

    # Step 6: Save PCA matrix
    save_pca_matrix(adata, args.output_pca)

    # Step 7: Perform downsampling and get indices
    adata_downsampled, selected_indices, selected_barcodes = (
        perform_downsampling_with_indices(
            adata, args.celltype_key, args.fraction, filter_nan
        )
    )

    # Step 8: Save indices
    save_indices(
        selected_indices, selected_barcodes, args.output_indices, args.indices_format
    )

    # Step 9: Optionally save processed AnnData
    if args.output_adata:
        save_processed_adata(adata_downsampled, args.output_adata)

    logger.info("=" * 60)
    logger.info("Combined pipeline completed successfully!")
    logger.info(f"PCA matrix saved to: {args.output_pca}")
    logger.info(f"Downsampled indices saved to: {args.output_indices}")
    logger.info(
        f"Selected {len(selected_indices)} cells out of {adata.n_obs} total cells"
    )
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
