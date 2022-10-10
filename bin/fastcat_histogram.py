#!/usr/bin/env python
"""Histogram-json."""

import argparse
import json

from aplanat.util import Limiter
import numpy as np
import pandas as pd


def histogram_counts(data, xlim_min=0, bin_width=100):
    """Histogram bins and counts."""
    x_lim = Limiter()
    x_lim.accumulate(data)
    bins = np.arange(xlim_min, x_lim.max + 1e-10, bin_width)
    counts, _ = np.histogram(data, bins=bins)
    return bins.tolist(), counts.tolist()


def get_stats(seq_summary):
    """Get Stats Json."""
    stats_json = {}
    len_data = seq_summary['read_length']
    len_bins, len_counts = histogram_counts(len_data, 0, 50)
    stats_json["len"] = dict(list(zip(len_bins, len_counts)))
    qual_data = seq_summary['mean_quality']
    qual_bins, qual_counts = histogram_counts(qual_data, 5, 0.2)
    stats_json["qual"] = dict(list(zip(qual_bins, qual_counts)))
    stats_json["total_reads"] = len(seq_summary)
    return stats_json


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--file", nargs='+', required=True,
        help="Read summary file.")
    parser.add_argument(
        "--sample_id", nargs=1
    )
    args = parser.parse_args()
    df = pd.read_csv(args.file[0], sep="\t")
    stats = get_stats(df)
    final = {args.sample_id[0]: stats}
    file_name = args.sample_id[0] + "stats.json"
    with open(file_name, 'w') as fp:
        json.dump(final, fp)


if __name__ == "__main__":
    main()
