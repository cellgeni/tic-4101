#!/usr/bin/env python3
"""
Run scib benchmark on downsampled AnnData with representation.

This script loads downsampled AnnData (with PCA), representation matrix, and
downsampled cell indices, attaches the representation to the AnnData object,
and runs scib benchmarking with pre-integrated PCA embedding.
"""

import argparse
import json
import logging
import sys

import anndata as ad
import numpy as np
import pandas as pd
from scib_metrics.benchmark import Benchmarker
from scib_metrics.benchmark import BioConservation, BatchCorrection


def setup_logging(log_level: str = "INFO", log_file: str = "benchmark.log"):
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


def load_representation(file_path: str) -> np.ndarray:
    """Load representation matrix from numpy file."""
    logger = logging.getLogger(__name__)
    logger.info(f"Loading representation matrix from: {file_path}")

    try:
        representation = np.load(file_path)
        logger.info(
            f"Successfully loaded representation with shape: {representation.shape}"
        )
        return representation
    except Exception as e:
        logger.error(f"Failed to load representation: {e}")
        raise


def load_downsampled_indices(file_path: str) -> tuple[list, list]:
    """Load downsampled cell indices from JSON file."""
    logger = logging.getLogger(__name__)
    logger.info(f"Loading downsampled indices from: {file_path}")

    try:
        with open(file_path, "r") as f:
            data = json.load(f)

        selected_indices = data["selected_indices"]
        selected_barcodes = data["selected_barcodes"]

        logger.info(f"Successfully loaded {len(selected_indices)} downsampled indices")
        logger.info(
            f"Total selected cells: {data.get('total_selected', len(selected_indices))}"
        )

        return selected_indices, selected_barcodes
    except Exception as e:
        logger.error(f"Failed to load downsampled indices: {e}")
        raise


def attach_representation_to_adata(
    adata: ad.AnnData,
    representation: np.ndarray,
    selected_indices: list,
    representation_key: str = "X_benchmark",
) -> ad.AnnData:
    """Attach representation matrix to downsampled AnnData object."""
    logger = logging.getLogger(__name__)
    logger.info("Attaching representation to downsampled AnnData object")

    # Verify that the representation matches the original dataset size
    # (before downsampling) and subset it to match the downsampled cells
    logger.info(f"Original representation shape: {representation.shape}")
    logger.info(f"Downsampled AnnData shape: {adata.shape}")
    logger.info(f"Number of selected indices: {len(selected_indices)}")

    # Subset representation to match downsampled cells
    representation_subset = representation[selected_indices]
    logger.info(f"Subset representation shape: {representation_subset.shape}")

    # Verify shapes match
    if representation_subset.shape[0] != adata.n_obs:
        logger.error(
            f"Representation subset shape {representation_subset.shape} doesn't match downsampled AnnData {adata.shape}"
        )
        raise ValueError(
            "Representation subset and downsampled AnnData have incompatible shapes"
        )

    # Add representation to obsm
    adata.obsm[representation_key] = representation_subset
    logger.info(f"Added representation to adata.obsm['{representation_key}']")

    logger.info(f"Final AnnData shape: {adata.shape}")
    logger.info(f"Available obsm keys: {list(adata.obsm.keys())}")

    return adata


def run_benchmark(
    adata: ad.AnnData,
    batch_key: str,
    celltype_key: str,
    embedding_keys: list,
    pre_integrated_key: str = "X_pca",
    n_jobs: int = 1,
    output_file: str = None,
) -> pd.DataFrame:
    """Run scib benchmark with pre-integrated embedding."""
    logger = logging.getLogger(__name__)
    logger.info("Starting scib benchmark")

    # Verify required keys exist
    if batch_key not in adata.obs.columns:
        available_keys = list(adata.obs.columns)
        logger.error(f"Batch key '{batch_key}' not found in adata.obs")
        logger.error(f"Available keys: {available_keys}")
        raise ValueError(f"Batch key '{batch_key}' not found in adata.obs")

    if celltype_key not in adata.obs.columns:
        available_keys = list(adata.obs.columns)
        logger.error(f"Celltype key '{celltype_key}' not found in adata.obs")
        logger.error(f"Available keys: {available_keys}")
        raise ValueError(f"Celltype key '{celltype_key}' not found in adata.obs")

    # Verify embedding keys exist
    for key in embedding_keys:
        if key not in adata.obsm.keys():
            available_keys = list(adata.obsm.keys())
            logger.error(f"Embedding key '{key}' not found in adata.obsm")
            logger.error(f"Available keys: {available_keys}")
            raise ValueError(f"Embedding key '{key}' not found in adata.obsm")

    # Verify pre-integrated key exists
    if pre_integrated_key not in adata.obsm.keys():
        available_keys = list(adata.obsm.keys())
        logger.error(
            f"Pre-integrated key '{pre_integrated_key}' not found in adata.obsm"
        )
        logger.error(f"Available keys: {available_keys}")
        raise ValueError(
            f"Pre-integrated key '{pre_integrated_key}' not found in adata.obsm"
        )

    logger.info(f"Batch key: {batch_key}")
    logger.info(f"Celltype key: {celltype_key}")
    logger.info(f"Embedding keys: {embedding_keys}")
    logger.info(f"Pre-integrated key: {pre_integrated_key}")
    logger.info(f"Number of jobs: {n_jobs}")

    # Log data distribution
    logger.info("Batch distribution:")
    batch_counts = adata.obs[batch_key].value_counts()
    for batch, count in batch_counts.items():
        logger.info(f"  {batch}: {count} cells")

    logger.info("Celltype distribution:")
    celltype_counts = adata.obs[celltype_key].value_counts()
    for celltype, count in celltype_counts.items():
        logger.info(f"  {celltype}: {count} cells")

    try:
        # Create benchmarker with pre-integrated embedding
        bm = Benchmarker(
            adata,
            batch_key=batch_key,
            label_key=celltype_key,
            bio_conservation_metrics=BioConservation(),
            batch_correction_metrics=BatchCorrection(),
            embedding_obsm_keys=embedding_keys,
            pre_integrated_embedding_obsm_key=pre_integrated_key,
            n_jobs=n_jobs,
        )

        logger.info("Running benchmark...")
        bm.benchmark()
        results = bm._results

        logger.info("Benchmark completed successfully!")
        logger.info(f"Results shape: {results.shape}")
        logger.info(f"Metrics: {list(results.columns)}")

        # Save results if output file specified
        if output_file:
            results.to_csv(output_file, index=True)
            logger.info(f"Saved benchmark results to: {output_file}")

        return results

    except Exception as e:
        logger.error(f"Benchmark failed: {e}")
        raise


def main():
    """Main function."""
    parser = argparse.ArgumentParser(
        description="Run scib benchmark on downsampled data with representation",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "adata_file", type=str, help="Path to downsampled AnnData file (.h5ad) with PCA"
    )

    parser.add_argument(
        "representation_file",
        type=str,
        help="Path to full representation matrix (.npy) from original dataset",
    )

    parser.add_argument(
        "downsampled_file", type=str, help="Path to downsampled indices (.json)"
    )

    parser.add_argument(
        "--batch_key",
        type=str,
        required=True,
        help="Column name in adata.obs containing batch information",
    )

    parser.add_argument(
        "--celltype_key",
        type=str,
        required=True,
        help="Column name in adata.obs containing cell type annotations",
    )

    parser.add_argument(
        "--representation_key",
        type=str,
        default="X_benchmark",
        help="Key to store representation in adata.obsm",
    )

    parser.add_argument(
        "--embedding_keys",
        type=str,
        nargs="+",
        help="List of embedding keys to benchmark (if not provided, will use representation_key)",
    )

    parser.add_argument(
        "--n_jobs",
        type=int,
        default=1,
        help="Number of parallel jobs for benchmarking",
    )

    parser.add_argument(
        "--output",
        type=str,
        default="benchmark_results.csv",
        help="Output file for benchmark results",
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
    logger.info("Starting scib benchmark pipeline")
    logger.info("=" * 60)
    logger.info(f"AnnData file: {args.adata_file}")
    logger.info(f"Representation file: {args.representation_file}")
    logger.info(f"Downsampled file: {args.downsampled_file}")
    logger.info(f"Batch key: {args.batch_key}")
    logger.info(f"Celltype key: {args.celltype_key}")
    logger.info(f"Representation key: {args.representation_key}")
    logger.info(f"Output file: {args.output}")

    # Step 1: Load all inputs
    adata = load_adata(args.adata_file)
    representation = load_representation(args.representation_file)
    selected_indices, selected_barcodes = load_downsampled_indices(
        args.downsampled_file
    )

    # Step 2: Attach representation to downsampled AnnData
    adata_benchmark = attach_representation_to_adata(
        adata,
        representation,
        selected_indices,
        args.representation_key,
    )

    # Step 3: Determine embedding keys
    if args.embedding_keys:
        embedding_keys = args.embedding_keys
    else:
        embedding_keys = [args.representation_key]

    logger.info(f"Using embedding keys: {embedding_keys}")

    # Step 4: Run benchmark with pre-integrated PCA
    results = run_benchmark(
        adata_benchmark,
        args.batch_key,
        args.celltype_key,
        embedding_keys,
        pre_integrated_key="X_pca",
        n_jobs=args.n_jobs,
        output_file=args.output,
    )

    logger.info("=" * 60)
    logger.info("Benchmark pipeline completed successfully!")
    logger.info(f"Results saved to: {args.output}")
    logger.info("=" * 60)

    # Print summary of results
    logger.info("Benchmark Results Summary:")
    logger.info("-" * 40)
    for metric in results.columns:
        for embedding in results.index:
            value = results.loc[embedding, metric]
            logger.info(f"{embedding} - {metric}: {value}")


if __name__ == "__main__":
    main()
