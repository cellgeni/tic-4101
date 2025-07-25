#!/usr/bin/env python3
"""
Preprocess AnnData object and generate PCA matrix.

This script performs normalization, HVG detection, and PCA
on single-cell RNA-seq data stored in AnnData format, and saves
the PCA matrix as a numpy array.
"""

import argparse
import logging
import sys

import anndata as ad
import scanpy as sc
import numpy as np


def setup_logging(log_level: str = "INFO", log_file: str = "preprocess_pca.log"):
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
        with open(barcodes_path, "w") as f:
            for barcode in adata.obs.index:
                f.write(f"{barcode}\n")
        logger.info(f"Saved cell barcodes to: {barcodes_path}")

    except Exception as e:
        logger.error(f"Failed to save PCA matrix: {e}")
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


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Preprocess AnnData and generate PCA matrix",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "input_file", type=str, help="Path to input AnnData file (.h5ad)"
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

    logger.info("=" * 60)
    logger.info("Starting AnnData preprocessing pipeline")
    logger.info("=" * 60)
    logger.info(f"Input file: {args.input_file}")
    logger.info(f"Apply log normalization: {args.lognorm}")
    logger.info(f"Use HVG: {args.hvg}")
    logger.info(f"Number of HVG: {args.n_hvg}")
    logger.info(f"Number of PCA components: {args.n_pca}")
    logger.info(f"Output PCA file: {args.output_pca}")

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

    # Step 5: Save PCA matrix
    save_pca_matrix(adata, args.output_pca)

    # Step 6: Optionally save processed AnnData
    if args.output_adata:
        save_processed_adata(adata, args.output_adata)

    logger.info("=" * 60)
    logger.info("Preprocessing pipeline completed successfully!")
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
