#!/usr/bin/env python
"""Create tables for the report."""
import json
import os
from pathlib import Path

from bokeh.models import HoverTool
import ezcharts as ezc
from ezcharts.plots.distribution import histplot
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
        abundance_threshold = round(
            abundance_threshold * df[column_to_filter].sum())
    # Group & filter them
    if column_to_group:
        # Subset just columns that are going to be used
        # E.g. not use lineages columns
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
    plt._fig.xaxis.axis_label = 'Quality score'
    plt._fig.x_range.start, plt._fig.x_range.end = min_qual, max_qual
    plt._fig.yaxis.axis_label = 'Number of reads'
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
    plt._fig.xaxis.axis_label = 'Read length / kb'
    plt._fig.yaxis.axis_label = 'Number of reads'
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
    m_depth = all_depth_pos.depth.mean().round(2)
    # Calculate the deviation from this mean
    s_depth = all_depth_pos.depth.std().round(2)
    return m_depth, s_depth, (s_depth/m_depth).round(2)


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


def depth2heatmap(depth, reference, min_nreads=0.01):
    """
    Calculate depth by windows for those references with a sequencing depth.

    :param depth (Pandas DataFrame): Sequencing depth per position for each reference.
    :param reference (Pandas DataFrame): Output from samtools coverage with the taxonomy
        assigned to each reference.
    :param min_nreads (float, optional): Minimum percentage of depth required to plot
        the reference in the heatmap. Defaults to 0.01.

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
    cutoff = min_nreads*(metrics.mean(axis=1).max())
    mask = metrics.mean(axis=1) > cutoff
    metrics = metrics[mask.values]
    return metrics


def load_alignment_data(align_stats, sample, rank='species'):
    """Load alignment data and return components for the report.

    :param align_stats (string): Path to the files
    :param sample (string): sample to be analyzed

    :return [table, scatter, heatmap]
    """
    rank = RANKS_ABB[rank]
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
    align_df.rename(
        columns={
            'numreads': 'number of reads',
            'coverage': 'pc coverage',
            'ref': 'reference',
            'endpos': 'ref length',
            }, inplace=True)

    # create the scatter plot.
    # x-axis is the number of reads, y-axis is the coverage.
    # Depth: number of reads that map in X position
    # Coverage: number of bases in the references that are covered by the reads.
    # Color: the rank (eg. species) of each reference to identify different
    # references of the same species. If species doesn't exist (e.g. silva db)
    # use the previous rank.
    if (rank == 'species') and (rank not in align_df.columns):
        rank = 'genus'
    plt_scatter = ezc.scatterplot(
        data=align_df,
        x='number of reads', y='pc coverage', hue=rank)
    plt_scatter._fig.title.text = sample
    # Add ref to the plot to be able to access later in hover
    # There is two glyd per species (rank), which is the hue
    # eg. 10 different families -> 20 glyd
    # Replace with the corresponding info based on the index value
    for gly in plt_scatter._fig.renderers:
        # add more properties for hover later
        gly.data_source.data['reference'] = align_df.loc[
            list(gly.data_source.to_df().index)]['reference']
        gly.data_source.data['rank'] = align_df.loc[
            list(gly.data_source.to_df().index)][rank]

    hover = plt_scatter._fig.select(dict(type=HoverTool))
    hover.tooltips = [
        ("number of reads", "$x"),
        ("coverage (%)", "$y"),
        ("reference", "@reference"),
        ("taxa", "@rank")
    ]
    # Add heatmap in blues to visualize the coverage of the reference and
    # the sequencing depth by quartiles (depth2heatmap).
    cmap = sns.color_palette(
        "Blues",  n_colors=1000, as_cmap=False).as_hex()
    plt_heatmap = ezc.heatmap(
        depth2heatmap(depth, stats),
        annot=False, cmap=cmap)
    plt_heatmap.tooltip = dict(position='top', trigger='item')
    plt_heatmap.yAxis.name = "Relative position"
    plt_heatmap.xAxis.axisLabel = dict(rotate=90)
    plt_heatmap.title = dict(text=sample)
    # return table, scatter, heatmap
    return [align_df, plt_scatter, plt_heatmap]


def n_reads_pass(metadata):
    """Table with number of reads discarded per sample."""
    n_reads = {}
    with open(metadata) as meta:
        meta_list = json.load(meta)
        for sample_meta in meta_list:
            n_reads[sample_meta["alias"]] = {
                "Reads": int(sample_meta.get("n_seqs")),
                "Reads after host depletion": sample_meta.get(
                    "n_seqs_passed_host_depletion"),
                "Unclassified": sample_meta.get("n_unclassified")
            }
    df = pd.DataFrame.from_dict(n_reads).T.sort_index()
    # NaNs are expected in the case of skipping steps
    # e.g. when the exclude_host is not used
    # remove those columns in the report
    df = df.dropna(axis=1, how='all')
    df.rename_axis('Sample alias')
    # Calculate percentages based on the initial number of reads (after fastcat)
    reference_column = 'Reads'
    cols_to_make_pc = df.columns.difference([reference_column])
    # Add pc suffix
    df[cols_to_make_pc + ' (%)'] = df[cols_to_make_pc].apply(
        lambda x: round(x / df[reference_column] * 100, 2))
    return df
