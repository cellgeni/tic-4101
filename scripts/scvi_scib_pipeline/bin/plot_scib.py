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

# Calculate space needed for parameter text
num_params = len(param_strings)
text_height_per_line = 0.03  # Height per line in figure coordinates
total_text_height = num_params * text_height_per_line + 0.02  # Extra padding

# Adjust the plot to make room for text at the bottom
plot_table.figure.subplots_adjust(bottom=total_text_height + 0.1)

# Add parameter strings as text at the bottom of the figure
for i, param_string in enumerate(param_strings):
    text = f"params{i}: {param_string}"
    # Position text at the bottom of the figure, starting from the bottom up
    y_pos = (num_params - i - 1) * text_height_per_line + 0.02
    plot_table.figure.text(
        0.02,  # Left margin
        y_pos,  # Bottom positioning
        text,
        transform=plot_table.figure.transFigure,
        ha="left",
        va="bottom",
        fontsize=10,
    )

plot_table.figure.savefig("scib_results_scaled.svg")

# Do the same for the unscaled plot
unscaled_plot = bm.plot_results_table(min_max_scale=False, show=False)

# Adjust the unscaled plot to make room for text at the bottom
unscaled_plot.figure.subplots_adjust(bottom=total_text_height + 0.1)

# Add parameter strings as text at the bottom of the unscaled figure
for i, param_string in enumerate(param_strings):
    text = f"params{i}: {param_string}"
    # Position text at the bottom of the figure, starting from the bottom up
    y_pos = (num_params - i - 1) * text_height_per_line + 0.02
    unscaled_plot.figure.text(
        0.02,  # Left margin
        y_pos,  # Bottom positioning
        text,
        transform=unscaled_plot.figure.transFigure,
        ha="left",
        va="bottom",
        fontsize=10,
    )

unscaled_plot.figure.savefig("scib_results.svg")
