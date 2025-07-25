import os
import sys
import glob
import argparse
import logging
import numpy as np
import scanpy as sc

from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

# setup logging to stdout
logging.basicConfig(
    stream=sys.stdout,
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s",
)


def init_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="A script for data integration using scVI"
    )
    parser.add_argument("--adata", type=str, help="path to procesesd adata file")
    parser.add_argument(
        "--representation_list", type=str, help="path to representation files"
    )
    parser.add_argument("--output_dir", type=str, help="path to output file")
    parser.add_argument("--n_jobs", type=int, default=1, help="number of jobs")
    return parser


def main():
    parser = init_parser()
    args = parser.parse_args()
    logging.info(f"Reading adata from {args.adata}")
    adata = sc.read_h5ad(args.adata)

    logging.info(f"Loading representations from {args.representation_list}")
    with open(args.representation_list, "r") as f:
        representation_paths = [p.strip() for p in f.readlines()]
        representation_paths_dict = {
            os.path.basename(os.path.dirname(path)): path
            for path in representation_paths
        }

    for key, path in representation_paths_dict.items():
        logging.info(f"Adding representation {key} from {path}")
        adata.obsm[key] = np.load(path)

    logging.info("Calculating PCA")
    sc.tl.pca(adata)

    scib_max_obs = 500000000
    if adata.shape[0] > scib_max_obs:
        sc.pp.subsample(adata, n_obs=scib_max_obs)

    benchmark_list = ["X_pca"] + [key for key in adata.obsm.keys() if "filtered" in key]
    logging.info(f"Benchmark embeddings: {benchmark_list}")

    bm = Benchmarker(
        adata,
        batch_key="donor_batch",
        label_key="broad_annotation",
        bio_conservation_metrics=BioConservation(),
        batch_correction_metrics=BatchCorrection(),
        embedding_obsm_keys=benchmark_list,
        n_jobs=args.n_jobs,
    )
    logging.info("Running benchmark")
    bm.benchmark()

    logging.info(f"Saving benchmark plots and results to {args.output_dir}")
    os.makedirs(args.output_dir, exist_ok=True)
    bm.plot_results_table(min_max_scale=False, save_dir=args.output_dir)
    bm.get_results(min_max_scale=False, clean_names=False).to_csv(
        os.path.join(args.output_dir, "benchmark_results.csv")
    )
    adata.write_h5ad(os.path.join(args.output_dir, "adata.h5ad"))
    logging.info("Done")


if __name__ == "__main__":
    main()
