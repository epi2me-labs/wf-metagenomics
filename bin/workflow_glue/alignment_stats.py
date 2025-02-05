#!/usr/bin/env python
"""Create tables for the report."""
from math import floor
from pathlib import Path

import numpy as np
import pandas as pd

from .util import wf_parser  # noqa: ABS101


def alignment_metrics(depth, stats):
    """Provide some alignment metrics for each reference.

    :param depth_tsv (DataFrame): Sequencing depth per position.
    :param stats (int): Coverage and number of reads for each reference.
    :return reference_stats (DataFrame): a df with all the stats and the mean +
        st. deviation in coverage for each reference.
    """
    # Input depth data may not cover all positions with zero coverage,
    # but we're interested in generating statistics that relate to the
    # entire reference. This code previously allocated memory to cover
    # all possible positions and then filled in the positions from the
    # loaded depth data. This not only double stored the input data, but
    # needlessly kept gigabytes of sparse depth in memory for the purpose
    # of calculating summary metrics.
    # Let's avoid all that. We can trivially calculate the mean as we
    # know the denominator. We can calculate the stddev with a pass over
    # the known data and adjust for the difference in mean after.
    ref_lens = dict(zip(stats.index, stats["endpos"]))

    grouped_stats = depth.groupby('ref')['depth'].agg(["count", "sum"])
    grouped_stats["ref_len"] = grouped_stats.index.map(ref_lens)
    grouped_stats["mean"] = grouped_stats["sum"] / grouped_stats["ref_len"]
    ref_means = dict(zip(grouped_stats.index, grouped_stats["mean"]))
    ref_sum_of_sqdifferences = {ref: 0 for ref in stats.index}

    # calculate the sample variance using the "direct" method
    # this is more robust than sum of squares, if we wanted to be fancy
    # we might write this again to just do a one-pass with Welford
    # but that'll need some rejigging in the caller scope.
    #
    # i'd written this to use kahan summation which was a bit more
    # involved but ultimately unnecessary - np.allclose was happy
    # and this method should suffice for the use case
    #   see https://www.johndcook.com/blog/2008/09/26/comparing-three-methods-of-computing-standard-deviation/  # noqa:E501
    for row in depth.itertuples():
        ref_sum_of_sqdifferences[row.ref] += (row.depth - ref_means[row.ref]) ** 2

    # account for missing zero depth positions
    # "count" is the number of positions observed, anything else is zero
    ref_skipped = dict(
        zip(grouped_stats.index, grouped_stats["ref_len"] - grouped_stats["count"])
    )
    for row in grouped_stats.itertuples():
        n_missing_zero = ref_skipped[row.Index]
        ref_sum_of_sqdifferences[row.Index] += (row.mean ** 2) * n_missing_zero

    grouped_stats["sd"] = \
        np.sqrt(
            grouped_stats.index.map(ref_sum_of_sqdifferences) /
            grouped_stats["ref_len"]  # over N for pop sd
        )
    grouped_stats["Coefficient of Variance"] = \
        grouped_stats["sd"] / grouped_stats["mean"]
    grouped_stats.drop(["count", "sum", "ref_len"], axis=1, inplace=True)

    # Add this df to reference_tsv as before
    reference_stats = pd.concat([stats, grouped_stats], axis=1)
    return reference_stats.reset_index()


def depth2heatmap(depth, reference, min_cov=1, windows=100):
    """
    Calculate depth by windows for those references with a sequencing depth.

    :param depth (Pandas DataFrame): Sequencing depth per position for each reference.
    :param reference (Pandas DataFrame): Output from samtools coverage with the taxonomy
        assigned to each reference.
    :param min_cov (float, optional): Minimum average depth required to plot
        the reference in the heatmap. Defaults to 1x.

    :return (DataFrame): heatmap data structure.
    """
    if min_cov < 0:
        raise ValueError("min_cov must be non-negative.")

    n_seqs = len(reference.index)
    # keep a map of reference names to position in our heatmap matrix
    ref_ids = {}
    # keep explicit array of ref_lens to avoid any assumptions on order
    # use an np array to support broadcasting division later
    ref_lens = np.zeros(n_seqs, dtype=np.uint64)

    # heatmap matrix - declared as float type to support div later
    ref_heatmap = np.zeros((n_seqs, windows), dtype=np.float64)
    for i, ref in enumerate(reference.itertuples()):
        ref_ids[ref.Index] = i
        ref_lens[i] = ref.endpos  # ref_lens are 1 based endpos

    # This previously used pd.qcut to discretise values, which unnecessarily
    # allocated a Series to assign each depth position a corresponding
    # window, and consumed the resulting dataframe to calculate some metrics
    # for each of the window groups. This used a lot of memory.
    # Let's avoid all that noise and simply iterate the positions; summing
    # the depth and later dividing by possible observations (to account for
    # missing zeros) to get the heatmap.
    for row in depth.itertuples():
        this_ref_id = ref_ids[row.ref]
        this_ref_len = ref_lens[this_ref_id]
        # convert row.pos to 0 based to ensure no window can be (window)
        this_window = floor((row.pos - 1) / this_ref_len * windows)
        ref_heatmap[this_ref_id, this_window] += row.depth

    # calculate the average over all windows for the ref
    # to determine the reference mask
    ref_mean_cov = ref_heatmap.sum(axis=1) / ref_lens
    ref_mask = ref_mean_cov >= min_cov

    # now convert window count cells to averages for plotting
    ref_heatmap /= (ref_lens // windows)[:, None]

    # apply the mask to remove refs that do not meet the threshold
    ref_heatmap = ref_heatmap[ref_mask]

    # convert to a dataframe with some nice labels to match the old api
    return pd.DataFrame(
        ref_heatmap,
        index=reference.index[ref_mask]
    )


def main(args):
    """Run the entry point."""
    cov_tsv = Path(args.coverage)
    depth_tsv = Path(args.depth)
    # the depth file is empty if there aren't classfied reads, e.g. barcode03
    if not cov_tsv.exists() or not depth_tsv.exists():
        pass
    else:
        # read the files and set the index
        stats = pd.read_csv(cov_tsv, sep='\t').rename(
            columns={'#rname': 'ref'}).set_index('ref')
        depth = pd.read_csv(
            depth_tsv, sep='\t', header=None, names=["ref", "pos", "depth"])
        # Second check to make sure the dataframe contains something
        if stats.empty:
            raise pd.errors.EmptyDataError(f"File is empty: {cov_tsv}")
        if depth.empty:
            raise pd.errors.EmptyDataError(f"File is empty: {depth_tsv}")
        # generate the table with some statistics per reference
        align_df = alignment_metrics(depth, stats)
        align_df['pcreads'] = (
            align_df['numreads']
            / align_df['numreads'].sum()
            * 100).round(2)
        align_df.rename(
            columns={
                'numreads': 'number of reads',
                'coverage': '% coverage',
                'ref': 'reference',
                'endpos': 'ref length',
                }, inplace=True)
        # Prepare also heatmap data
        heatmap_matrix = depth2heatmap(depth, stats)
        # Write both tables to be imported later in report.py
        align_df.to_csv(args.output, sep='\t', index=False)
        heatmap_matrix.to_csv(args.output_heatmap, sep='\t')


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("alignment_table")
    parser.add_argument(
        "--output", help="Name of the output table.",
        type=Path
    )
    parser.add_argument(
        "--output_heatmap", help="Name of the heatmap table.",
        type=Path
    )
    parser.add_argument(
        "--coverage", help="Output file of samtools coverage.")
    parser.add_argument(
        "--depth", help="Output file of samtools depth")
    return parser
