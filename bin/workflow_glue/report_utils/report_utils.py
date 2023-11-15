#!/usr/bin/env python
"""Create tables for the report."""
import collections
import json
import os
from pathlib import Path

import ezcharts as ezc
from ezcharts.plots.distribution import histplot
import numpy as np
import pandas as pd
import seaborn as sns
import workflow_glue.diversity as diversity

RANK_ORDER = {
    'superkingdom': 1, 'kingdom': 2, 'phylum': 3, 'class': 4, 'order': 5, 'family': 6,
    'genus': 7, 'species': 8}

RANKS = list(RANK_ORDER.keys())

RANKS_NO_SK_K = list(RANK_ORDER.keys())[2:]

RANKS_ABB = {
    "D": "superkingdom", "K": "kingdom", "P": "phylum", "C": "class",
    "O": "order", "F": "family", "G": "genus", "S": "species"
}

# READ INPUT DATA


def parse_lineages(lineages):
    """Join lineage counts into a df for a particular taxonomic rank.

    :param lineages (str): String indicating the path of the dir with the json.
    :return (dict): Taxa counts in json structure.
    """
    path_to_json = lineages
    json_files = [
        pos_json for pos_json in os.listdir(
            path_to_json) if pos_json.endswith('.json')]
    all_json = {}
    for i in json_files:
        with open(os.path.join(lineages, i)) as json_file:
            # process sankey and table counts
            all_json.update(json.load(json_file))
    return all_json


def prepare_data_to_sunburst(lineages_sample, new_lineages=[], position='outside'):
    """Join lineage counts into a df for a particular taxonomic rank.

    :param lineages_sample (str): Dictionary with lineages of one sample.
    :return (dict): Taxa counts in json structure ready for sunburst.
    """
    for taxon, taxon_data in lineages_sample.items():
        # replace single quotes to double quotes to be escaped in ezcharts
        # (as some species names in NCBI taxdmp files contain single quotes,
        # e.g. "'Sphingomonas aeria' Park et al. 2015 non Xue et al. 2018")
        new_lineages.append(dict(
            name=taxon.replace("'", '"'), value=taxon_data["count"]
        ))
        if bool(taxon_data["children"]):
            new_lineages[-1].update(children=[])
            prepare_data_to_sunburst(
                taxon_data["children"], new_lineages[-1]["children"])
        else:
            # add position to make the labels be outside of the outer circle in the plot
            new_lineages[-1].update(label=dict(position=position))
    return new_lineages

# PREPARE INPUT DATA


def split_taxonomy_string(counts_per_taxa_df, set_index=False):
    """Divide the taxonomy string into columns for each different taxonomic rank.

    :param counts_per_taxa_df (DataFrame): Abundance dataframe with taxa in
            rows and samples in columns.
    :param set_index (bool): Set the idx to the last rank for tables in the report.
    :return (DataFrame): Return a new dataFrame with counts per specific rank with
            a new column per specific rank.
    """
    taxonomy_ranks = counts_per_taxa_df.tax.str.split(';', expand=True)
    # Convert numbers to normal rank names
    taxonomy_ranks.columns = [RANKS[i] for i in range(taxonomy_ranks.shape[1])]
    df_ranks = pd.concat([counts_per_taxa_df, taxonomy_ranks], axis=1)
    if set_index:
        ranks = taxonomy_ranks.columns
        df_ranks = df_ranks.set_index(ranks[-1])
    return df_ranks


def most_abundant_table(counts_per_taxa_df, n=10, percent=False):
    """Select most abundant taxa and group the rest into 'Other' taxa.

    :param counts_per_taxa_df (DataFrame): DataFrame with counts per specific rank.
    :param samples (set): Set of samples.
    :param n (int, optional): Number of abundant taxa  to be chosen. Defaults to 10.
    :param percent (bool, optional): Make the percentage by column. Defaults to False.
    :return (DataFrame): DataFrame with most abundant taxa (in the total of all the\
            samples) and the rest group in Others.
    """
    # sort table - just in case
    counts_per_taxa_df = counts_per_taxa_df.sort_values(by=['total'], ascending=False)
    # Extract n most abundant taxa & group less abundant in Others
    most_abundant_taxa = counts_per_taxa_df[:n]
    less_abundant_taxa = counts_per_taxa_df[n:].sum().to_dict()
    less_abundant_taxa['tax'] = 'Other' + ';Other' * (
        len(most_abundant_taxa['tax'].iloc[0].split(';'))-1)
    if less_abundant_taxa['total'] > 1:  # no "other" taxa
        other = pd.DataFrame(less_abundant_taxa,  index=['Other'])
        # Concat other to most abundant, this is data for plotting
        d2plot = pd.concat([most_abundant_taxa, other])
    else:
        d2plot = most_abundant_taxa

    if percent:  # calculate percent
        d2plot = d2plot.set_index('tax')
        d2plot = d2plot.apply(
            lambda x: 100 * x / x.sum()).round(2)
        d2plot = d2plot.reset_index()
    return d2plot


def per_sample_stats(stats):
    """Read from fastcat stats and return total reads per sample.

    :param stats (list): List with fastcat statistics.
    :return (DataFrame): Reads per sample.
    """
    dfs = collections.Counter()
    for fastcat_file in stats:
        df = pd.read_csv(
            fastcat_file,
            sep="\t",
            header=0
            )
        # save total number of reads per sample
        dfs += collections.Counter(df['sample_name'])
    df_allstats = pd.DataFrame.from_dict(
        dfs, orient='index').reset_index()
    df_allstats.columns = ['sample_name', 'Number of reads']
    return df_allstats.sort_values(by=['sample_name'])


def calculate_diversity_metrics(counts_per_taxa_df):
    """Generate diversity indexes from abundance table.

    :param counts_per_taxa_df (DataFrame): Counts per taxa and sample.
    :return  (DataFrame): Diversity metrics for a specific rank
    """
    div = diversity.alpha_diversity(counts_per_taxa_df).round(2).sort_index(
        ).reset_index(level=0)
    div.fillna('None', inplace=True)
    div.rename(columns={
        'index': 'Sample', 'S': 'Richness', 'H': 'Shannon Diversity Index (H)'
        }, inplace=True)
    return div


def filter_by_abundance(
        df, column_to_filter, abundance_threshold=0, column_to_group=None):
    """Given a df, return a filtered dataframe after applying a threshold of abundances.

    :param df (DataFrame): Dataframe with counts.
    :param column_to_filter (string): name of the column in which apply the filter.
    :param abundance_threshold (int or float, optional): If it is a natural number it
        remove those rows with less counts. Decimal between 0-1 are treated as
        percentages. Defaults to 0.
    :param column_to_group (string): name of the column to make groups before filtering.

    :return (DataFrame): Dataframe with rows that satisfy the threshold.
    """
    if abundance_threshold < 1:
        abundance_threshold = abundance_threshold * round(
            df[column_to_filter].sum())
    # Group & filter them
    if column_to_group:
        # Subset just columns that are going to be used
        # In case the df contains character/factor columns
        interesting_cols = [
            colname for colname in [
                column_to_filter, column_to_group
                ] if colname in df.columns]
        mini_df = df[interesting_cols]
        mini_df = mini_df.groupby(column_to_group).sum()
        df_filtered = df[
            df[column_to_group].isin(  # list of species that satisfy the threshold
                mini_df.loc[mini_df[column_to_filter] > abundance_threshold].index
            )]
    else:
        df_filtered = df[df[column_to_filter] > abundance_threshold]
    df_filtered = df_filtered.dropna()
    return df_filtered

# PLOTS

# The SeqSummary from ezcharts.components.fastcat cannot be used.
# It groups data into bins, but from the real time analysis output
# the input data is already grouped into bins.


def read_quality_plot(seq_summary, min_qual=4, max_qual=30, title='Read quality'):
    """Create read quality summary plot."""
    df = pd.DataFrame.from_dict(seq_summary['qual'].items())
    df.columns = ['mean_quality', 'counts']
    df['mean_quality'] = df['mean_quality'].astype('float')
    plt = histplot(
        data=df['mean_quality'],
        bins=len(df),
        weights=list(df['counts'])
        )
    plt.title = dict(text=title)
    plt.xAxis.name = 'Quality score'
    plt.xAxis.min, plt.xAxis.max = min_qual, max_qual
    plt.yAxis.name = 'Number of reads'
    return plt


def read_length_plot(seq_summary, title='Read length'):
    """Create a read length plot."""
    df = pd.DataFrame.from_dict(seq_summary['len'].items())
    df.columns = ['read_length', 'counts']
    df['read_length'] = df['read_length'].astype('uint64')
    df['read_length'] = df['read_length'] / 1000
    plt = histplot(
        data=df['read_length'],
        bins=len(df),
        weights=list(df['counts']))
    plt.title = dict(text=title)
    plt.xAxis.name = 'Read length / kb'
    plt.yAxis.name = 'Number of reads'
    return plt


def parse_amr(amr_dir):
    """Join sample jsons together.

    :param amr_dir (str): String indicating the path of the dir with amr json results.
    :return (dict): Grouped results in json structure.
    """
    all_output = dict()
    for file in os.listdir(amr_dir):
        with open(f"{amr_dir}/{file}") as fh:
            data = json.load(fh)
        all_output.update(data)
    return all_output


def nreads_all_positions(depth_ref, ref_len):
    """For each reference, provide mean and standard deviation of the sequencing depth.

    :param depth (Pandas DataFrame): Sequencing depth per position per reference.
    :param ref_len (int): length in bp of the reference.
    :return all_depth_pos, DataFrame: contains the depth of each position for reference.
    """
    # Initiate an empty dataframe with all the positions of the reference.
    # Samtools depth returns just positions that at least had 1 read mapped,
    # exclude those positions whose depth is 0 to avoid huge files.
    # To make percentiles, we need to count these empty positions,
    # but better do it now by each reference than having a file with all pos for all the
    # references.
    all_depth_pos = pd.DataFrame(0, index=range(1, ref_len + 1), columns=["depth"])
    all_depth_pos.index.name = "pos"
    # Fill the df with the sequencing depth of each position.
    all_depth_pos.loc[depth_ref["pos"], 'depth'] = list(depth_ref["depth"])
    # total_pos contains now the depth of each position, even those whose depth=0
    return all_depth_pos


def coverage_dispersion(depth_ref, ref_len):
    """For each reference, provide mean and standard deviation of the sequencing depth.

    :param depth (Pandas DataFrame): Sequencing depth per position per reference.
    :param ref_len (int): length in bp of the reference.
    :return average_depth, st. deviation, deviation/ref_len (tuple): Average of the
        sequencing depth, deviation from this average, coefficient of variation.
    """
    # Get sequencing depth at each positions
    all_depth_pos = nreads_all_positions(depth_ref, ref_len)
    # Calculate the mean of the reads distribution
    m_depth = all_depth_pos.depth.mean()
    # Calculate the deviation from this mean
    s_depth = all_depth_pos.depth.std()
    return m_depth, s_depth, s_depth/m_depth


def alignment_metrics(depth, stats):
    """Provide some alignment metrics for each reference.

    :param depth_tsv (DataFrame): Sequencing depth per position.
    :param stats (int): Coverage and number of reads for each reference.
    :return reference_stats (DataFrame): a df with all the stats and the mean +
        st. deviation in coverage for each reference.
    """
    # For each row/reference, get the sequence length and positional depth.
    # Then, calculate the depth and coverage dispersion.
    # Each row is a reference, as in the samtools coverage output ->
    # apply 1 or ‘columns’: apply function to each row.
    metrics = stats.apply(lambda row: coverage_dispersion(
        depth.query("ref == @row.name"), stats.loc[
            row.name, 'endpos']), axis=1)
    metrics = pd.DataFrame(metrics.tolist(), index=metrics.index)
    metrics.columns = ['mean', 'sd', 'Coefficient of Variance']

    # Add this df to reference_tsv
    reference_stats = pd.concat([stats, metrics], axis=1)
    return reference_stats.reset_index()


def depth_windows(depth_ref, ref_len, nwindows=100):
    """For each reference, return mean depth by windows.

    :param depth (Pandas DataFrame): Sequencing depth per position for each reference.
    :param ref_len (int): length in bp of the reference.
    :param nwindows (int): Number of windows to make intervals in the length of the
        reference.
    """
    # Get sequencing depth of all the positions
    all_depth_pos = nreads_all_positions(depth_ref, ref_len)
    # There is likely a very large variability in reference lengths: virus + Eukaryota
    # Calculate sequencing depth for each n quantile (nwindows) -> break in intervals
    all_depth_pos['windows'] = pd.qcut(all_depth_pos.index, nwindows)
    # Sequencing depth in each window
    windows = all_depth_pos.groupby(['windows']).mean()
    # nwindows is the same for all the references regardless each length.
    # use these nwindows as the axis to plot.
    windows.index = list(range(len(windows)))
    # return the series (it will be added to a df in the depth2heatmap function to make
    # the df to be plotted).
    return pd.Series(windows['depth'])


def depth2heatmap(depth, reference, min_nreads=0.001):
    """
    Calculate depth by windows for those references with a sequencing depth.

    :param depth (Pandas DataFrame): Sequencing depth per position for each reference.
    :param reference (Pandas DataFrame): Output from samtools coverage with the taxonomy
        assigned to each reference.
    :param min_nreads (float, optional): Minimum percentage of depth required to plot
        the reference in the heatmap. Defaults to 0.001.

    :return plot.
    """
    # For each row/reference, get the sequence length and positional depth.
    # Then, calculate the depth  in sliding windows (percentiles).
    # Each row is a reference, as in the samtools coverage output ->
    # apply 1 or ‘columns’: apply function to each row.
    metrics = reference.apply(lambda row: depth_windows(
        depth.query("ref == @row.name"), reference.loc[
            row.name, 'endpos']), axis=1)
    # Apply an abundance cutoff for the heatmap.
    # Get references with depth greater than cutoff.
    cutoff = min_nreads*(metrics.sum(axis=1).max())
    mask = metrics.mean(axis=1) > cutoff
    metrics = metrics[mask.values]
    return metrics


def load_alignment_data(align_stats, sample, rank='species'):
    """Load alignment data and return components for the report.

    :param align_stats (string): Path to the files
    :param sample (string): sample to be analyzed

    :return [table, scatter, heatmap]
    """
    cov_tsv = Path(f"{align_stats}/{sample}.reference.tsv.gz")
    depth_tsv = Path(f"{align_stats}/{sample}.depth.tsv.gz")
    # the depth file is empty if there aren't classfied reads, e.g. barcode03
    if not cov_tsv.exists() or not depth_tsv.exists():
        return None
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
    align_df.rename(columns={'numreads': 'number of reads'}, inplace=True)
    # create the scatter plot.
    # x-axis is the number of reads, y-axis is the coverage.
    # Depth: number of reads that map in X position
    # Coverage: number of bases in the references that are covered by the reads.
    # Color: by default the species of each reference to help to identify different
    # references of the same species. If 'species doesn't exist (e.g. silva db),
    # use the previous rank.
    if 'species' not in align_df.columns:
        rank = 'genus'
    plt_scatter = ezc.scatterplot(
        data=align_df.reset_index(),
        x='number of reads', y='coverage', hue=rank)
    # move to ezcharts package!!!!
    # Add custom label tooltips: nreads, cov, species name, reference (ref)
    plt_scatter.dataset[0].source = np.hstack((
        plt_scatter.dataset[0].source,
        align_df['ref'].values[:, None]))
    plt_scatter.dataset[0].dimensions = list(
        plt_scatter.dataset[0].dimensions) + ['ref']
    # Add them and custom format-color by species
    for series in plt_scatter.series:
        series.encode = dict(
            x="number of reads",
            y="coverage",
            seriesName=None,
            itemName="ref",
            tooltip=['x', 'y'],
        )
        plt_scatter.tooltip = {"trigger": "item"}
    # Add heatmap in blues to visualize the coverage of the reference and
    # the sequencing depth by quartiles (depth2heatmap).
    cmap = sns.color_palette(
        "Blues",  n_colors=1000, as_cmap=False).as_hex()
    plt_heatmap = ezc.heatmap(
        depth2heatmap(depth, stats),
        annot=False, cmap=cmap)
    plt_heatmap.tooltip = dict(position='top', trigger='item')
    plt_heatmap.yAxis.name = "Relative position"
    # return table, scatter, heatmap
    return [align_df, plt_scatter, plt_heatmap]
