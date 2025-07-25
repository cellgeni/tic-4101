#!/usr/bin/env python3
"""
Downsample cells based on cell type annotations.

This script performs downsampling on single-cell RNA-seq data
and returns the indices of selected cells in JSON format.
"""

import argparse
import json
import logging
import sys

import anndata as ad
from celltypist.samples import downsample_adata


def setup_logging(log_level: str = "INFO", log_file: str = "downsample_cells.log"):
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


def perform_downsampling_with_indices(
    adata: ad.AnnData, celltype_key: str, fraction: float
) -> tuple[ad.AnnData, list, list]:
    """Perform downsampling using celltypist and return indices."""
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
    logger.info("Cell type distribution before downsampling:")
    for celltype, count in celltype_counts.items():
        logger.info(f"  {celltype}: {count} cells")

    # Store original cell indices (barcodes)
    original_barcodes = list(adata.obs.index)

    # Perform downsampling with return_index=True to get indices
    selected_indices = downsample_adata(
        adata=adata,
        mode="total",
        n_cells=n_cells,
        by=celltype_key,
        random_state=4,
        return_index=True,
    )

    # Create downsampled adata using the indices
    adata_downsampled = adata[selected_indices].copy()

    logger.info(f"Downsampling completed. New shape: {adata_downsampled.shape}")
    logger.info(f"Selected {len(selected_indices)} cell indices")

    # Log cell type distribution after downsampling
    celltype_counts_after = adata_downsampled.obs[celltype_key].value_counts()
    logger.info("Cell type distribution after downsampling:")
    for celltype, count in celltype_counts_after.items():
        logger.info(f"  {celltype}: {count} cells")

    # Convert numpy array to Python list for JSON serialization
    selected_indices_list = selected_indices.tolist()

    # Get the selected cell barcodes
    selected_barcodes = [original_barcodes[i] for i in selected_indices]

    return adata_downsampled, selected_indices_list, selected_barcodes


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
            with open(output_path, "w") as f:
                json.dump(data, f, indent=2)

        elif format_type.lower() == "txt":
            # Save as simple text file with barcodes
            with open(output_path, "w") as f:
                for barcode in barcodes:
                    f.write(f"{barcode}\n")

        elif format_type.lower() == "csv":
            # Save as CSV with index and barcode columns
            import csv

            with open(output_path, "w", newline="") as f:
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


def save_downsampled_adata(adata: ad.AnnData, output_path: str):
    """Save downsampled AnnData object (optional)."""
    logger = logging.getLogger(__name__)
    logger.info(f"Saving downsampled AnnData to: {output_path}")

    try:
        adata.write_h5ad(output_path)
        logger.info(f"Successfully saved downsampled AnnData with shape: {adata.shape}")
    except Exception as e:
        logger.error(f"Failed to save downsampled AnnData: {e}")
        raise


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Downsample cells and save indices",
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
        "--output_indices",
        type=str,
        default="downsampled_indices.json",
        help="Output file path for cell indices",
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
        help="Optional: Save downsampled AnnData object to this path",
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
    logger.info("Starting cell downsampling pipeline")
    logger.info("=" * 60)
    logger.info(f"Input file: {args.input_file}")
    logger.info(f"Celltype key: {args.celltype_key}")
    logger.info(f"Downsample fraction: {args.fraction}")
    logger.info(f"Output indices file: {args.output_indices}")
    logger.info(f"Indices format: {args.indices_format}")

    # Step 1: Load data
    adata = load_adata(args.input_file)

    # Step 2: Perform downsampling and get indices
    adata_downsampled, selected_indices, selected_barcodes = (
        perform_downsampling_with_indices(adata, args.celltype_key, args.fraction)
    )

    # Step 3: Save indices
    save_indices(
        selected_indices, selected_barcodes, args.output_indices, args.indices_format
    )

    # Step 4: Optionally save downsampled AnnData
    if args.output_adata:
        save_downsampled_adata(adata_downsampled, args.output_adata)

    logger.info("=" * 60)
    logger.info("Downsampling pipeline completed successfully!")
    logger.info(
        f"Selected {len(selected_indices)} cells out of {adata.n_obs} total cells"
    )
    logger.info("=" * 60)


if __name__ == "__main__":
    main()
