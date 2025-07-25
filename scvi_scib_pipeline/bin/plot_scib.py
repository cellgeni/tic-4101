#!/usr/bin/env python3

import os
import glob
import argparse
import pandas as pd
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection

parser = argparse.ArgumentParser(description="Plot SCIB results.")
parser.add_argument(
    "--filedir", type=str, help="Path to directory containing result files"
)
args = parser.parse_args()

result_list = []
files = sorted(glob.glob(f"{args.filedir}/*.csv"))
param_strings = list(map(lambda x: os.path.basename(x).replace(".csv", ""), files))
obsm_keys = [f"params{i}" for i in range(len(files))]
for i, file in enumerate(files):
    df = pd.read_csv(file, index_col=0)
    df.columns = [f"params{i}", df.columns.tolist()[1]]
    result_list.append(df.iloc[:, :-1])
result_list.append(df.iloc[:, -1:])

bm = Benchmarker(
    None,
    None,
    None,
    [None],
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
)

bm._results = pd.concat(result_list, axis=1)
bm._embedding_obsm_keys = obsm_keys
bm._benchmarked = True

plot_table = bm.plot_results_table(show=False)
plot_table.figure.savefig("scib_results_scaled.svg")

plot_table = bm.plot_results_table(min_max_scale=False, show=False)
plot_table.figure.savefig("scib_results.svg")

bm._results.to_csv("scib_results.csv")
pd.DataFrame({"params": param_strings, "obsm_keys": obsm_keys}).to_csv(
    "scib_params.csv", index=False
)
