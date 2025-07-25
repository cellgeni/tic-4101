#!/usr/bin/env python3
"""
Downsample AnnData object with preprocessing steps.

This script performs normalization, HVG detection, PCA, and downsampling
on single-cell RNA-seq data stored in AnnData format.
"""

import argparse
import logging
import sys
from pathlib import Path

import anndata as ad
import scanpy as sc
import celltypist
from celltypist.samples import downsample_adata


def setup_logging(log_level: str = "INFO"):
    """Set up logging configuration."""
    logging.basicConfig(
        level=getattr(logging, log_level.upper()),
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(sys.stdout),
            logging.FileHandler("downsample.log"),
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

    # Scale data
    # sc.pp.scale(adata, max_value=10)
    # logger.info("Scaled data (max_value=10)")

    # Calculate PCA
    sc.tl.pca(
        adata,
        n_comps=n_comps,
        svd_solver="arpack",
        use_highly_variable=use_highly_variable,
    )
    logger.info(f"Calculated PCA with {adata.obsm['X_pca'].shape[1]} components")

    return adata


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


def perform_downsampling(
    adata: ad.AnnData, celltype_key: str, fraction: float
) -> ad.AnnData:
    """Perform downsampling using celltypist."""
    logger = logging.getLogger(__name__)

    # Calculate target number of cells based on fraction
    total_cells = adata.n_obs
    n_cells = int(total_cells * fraction)

    logger.info(
        f"Starting downsampling with celltype_key='{celltype_key}', fraction={fraction}"
    )
    logger.info(f"Total cells: {total_cells}, Target cells: {n_cells}")

    # Check if celltype_key exists
    if celltype_key not in adata.obs.columns:
        available_keys = list(adata.obs.columns)
        logger.error(f"Celltype key '{celltype_key}' not found in adata.obs")
        logger.error(f"Available keys: {available_keys}")
        raise ValueError(f"Celltype key '{celltype_key}' not found in adata.obs")

    # Validate fraction
    if not 0 < fraction <= 1:
        raise ValueError(f"Fraction must be between 0 and 1, got {fraction}")

    # Log cell type distribution before downsampling
    celltype_counts = adata.obs[celltype_key].value_counts()
    logger.info(f"Cell type distribution before downsampling:")
    for celltype, count in celltype_counts.items():
        logger.info(f"  {celltype}: {count} cells")

    # Perform downsampling
    adata_downsampled = downsample_adata(
        adata=adata,
        mode="total",
        n_cells=n_cells,
        by=celltype_key,
        random_state=4,
        return_index=False,
    )

    logger.info(f"Downsampling completed. New shape: {adata_downsampled.shape}")

    # Log cell type distribution after downsampling
    celltype_counts_after = adata_downsampled.obs[celltype_key].value_counts()
    logger.info(f"Cell type distribution after downsampling:")
    for celltype, count in celltype_counts_after.items():
        logger.info(f"  {celltype}: {count} cells")

    return adata_downsampled


def save_adata(adata: ad.AnnData, output_path: str):
    """Save AnnData object."""
    logger = logging.getLogger(__name__)
    logger.info(f"Saving downsampled AnnData to: {output_path}")

    try:
        adata.write_h5ad(output_path)
        logger.info(f"Successfully saved AnnData with shape: {adata.shape}")
    except Exception as e:
        logger.error(f"Failed to save AnnData: {e}")
        raise


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Downsample AnnData object with preprocessing",
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
        "--hvg", action="store_true", help="Find highly variable genes before PCA"
    )

    parser.add_argument(
        "--lognorm",
        action="store_true",
        help="Apply normalization and log1p transformation",
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
        "--output", type=str, default="downsampled.h5ad", help="Output file path"
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

    logger.info("=" * 60)
    logger.info("Starting AnnData downsampling pipeline")
    logger.info("=" * 60)
    logger.info(f"Input file: {args.input_file}")
    logger.info(f"Celltype key: {args.celltype_key}")
    logger.info(f"Downsample fraction: {args.fraction}")
    logger.info(f"Use HVG: {args.hvg}")
    logger.info(f"Apply log normalization: {args.lognorm}")
    logger.info(f"Output file: {args.output}")

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

    # Step 6: Perform downsampling
    adata_downsampled = perform_downsampling(adata, args.celltype_key, args.fraction)

    # Step 7: Save result
    save_adata(adata_downsampled, args.output)

    logger.info("=" * 60)
    logger.info("Pipeline completed successfully!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
