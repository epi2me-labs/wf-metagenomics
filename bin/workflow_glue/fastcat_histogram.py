#!/usr/bin/env python
"""Histogram-json."""

import json
import sys

import numpy as np
import pandas as pd
from .util import wf_parser  # noqa: ABS101


def histogram_counts(data, dmin=0, bin_width=100):
    """Histogram bins and counts."""
    bins = np.arange(dmin, max(data) + bin_width, bin_width)
    counts, _ = np.histogram(data, bins=bins)
    return bins.tolist(), counts.tolist()


def get_stats(seq_summary):
    """Get Stats Json."""
    stats_json = {
        "total_reads": len(seq_summary)}
    if len(seq_summary) > 1:
        len_data = seq_summary['read_length']
        len_bins, len_counts = histogram_counts(
            len_data, dmin=0, bin_width=50)
        stats_json["len"] = dict(list(zip(len_bins, len_counts)))

        qual_data = seq_summary['mean_quality']
        qual_bins, qual_counts = histogram_counts(
            qual_data, dmin=0, bin_width=0.2)
        stats_json["qual"] = dict(list(zip(qual_bins, qual_counts)))
    else:
        sys.stderr.write("WARNING: summary file was empty.\n")
        stats_json["len"] = dict()
        stats_json["qual"] = dict()
    return stats_json


def main(args):
    """Run the entry point."""
    df = pd.read_csv(
        args.input, sep="\t",
        usecols=['read_length', 'mean_quality'],
        dtype={'read_length': int, 'mean_quality': float})
    final = {args.sample_id: get_stats(df)}
    with open(args.output, 'w') as fp:
        json.dump(final, fp)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("fastcat_histogram")
    parser.add_argument(
        "input", help="Read summary file.")
    parser.add_argument(
        "output", help="Output summary JSON.")
    parser.add_argument(
        "--sample_id", help="Sample name.")
    return parser
