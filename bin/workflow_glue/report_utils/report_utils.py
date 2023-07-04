#!/usr/bin/env python
"""Create tables for the report."""

import json
import os
import re

import anytree
from ezcharts.plots.distribution import histplot
import pandas as pd
import workflow_glue.diversity as diversity

RANK_ORDER = {
    'superkingdom': 1, 'kingdom': 2, 'phylum': 3, 'class': 4, 'order': 5, 'family': 6,
    'genus': 7, 'species': 8}

RANKS = list(RANK_ORDER.keys())

RANKS_NO_SK_K = list(RANK_ORDER.keys())[2:]

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
        if taxon_data.get('count'):
            new_lineages.append(dict(name=taxon, value=taxon_data["count"]))
        if bool(taxon_data["children"]):
            new_lineages[-1].update(children=[])
            prepare_data_to_sunburst(
                taxon_data["children"], new_lineages[-1]["children"])
        else:
            # add position to make the labels be outside of the outer circle in the plot
            new_lineages[-1].update(label=dict(position=position))
    return new_lineages

# PREPARE INPUT DATA


def tax_tree(lineage_trees_dict):
    """From lineages json, create a dictionary with {sample:tree (counts per lineage)}.

    Save the full taxonomy string.

    :param d (dict): Taxa counts in json structure. Nested dictionary.
        {"Sample": {"Taxon name":{"rank":str, "count":int, "children": dict},}}.
    :return (dict): Dictionary with one lineage tree per sample.
    """
    # Root project
    trees = {}
    for sample_id, lineage_sample in lineage_trees_dict.items():
        # Root sample
        rootnode = anytree.Node("root", rank=None, count=0, sample=None)

        def itertaxa(d, parent=None, parent_node=None):
            """From lineages json, construct a tree with the counts per each taxon.

            :param d (dict): Taxa counts in json structure. Nested dictionary:
                {"Taxon":{"rank":str, "count":int, "children": {}}}}.
            :param parent (str): Name of the parent taxon.
            :param parent_node (Node): Node with the info related to the parent taxon.
            """
            for taxon, taxon_data in d.items():
                if parent:
                    # Databases don't always follow a perfect structure.
                    # Check that parent belongs to immediate previous rank
                    # to avoid those cases in which some ranks are missing
                    # (e.g.: Wohlfahrtiimonas)
                    diff_ranks = RANK_ORDER[taxon_data['rank']] - parent_node.rank
                    if diff_ranks == 1:  # Consecutive ranks
                        node = anytree.Node(
                            taxon, parent=parent_node,
                            count=taxon_data['count'],
                            rank=RANK_ORDER[taxon_data['rank']],
                            )
                    else:  # Add missing ranks
                        original_parent_node = parent_node
                        for i in range(1, diff_ranks):
                            # Fill missing ranks with label: Incertae sedis,
                            # i.e. https://en.wikipedia.org/wiki/Wohlfahrtiimonas
                            label = 'Incertae_sedis'
                            if (
                                'Bacteria' in original_parent_node.name or
                                    'Archaea' in original_parent_node.name):
                                # Leave in blank
                                node = anytree.Node(
                                    'none',
                                    parent=parent_node,
                                    count=taxon_data['count'],
                                    rank=(parent_node.rank + 1)
                                    )
                            elif 'Unclassified' not in original_parent_node.name:
                                node = anytree.Node(
                                    f'{parent_node.name}_{label}',
                                    parent=parent_node,
                                    count=taxon_data['count'],
                                    rank=(parent_node.rank + 1)
                                    )
                            else:  # Do not do this for Unclassified
                                node = anytree.Node(
                                    taxon,
                                    parent=parent_node,
                                    count=taxon_data['count'],
                                    rank=(parent_node.rank + 1)
                                    )
                            parent_node = node
                        # Add last identified rank
                        node = anytree.Node(
                            taxon, parent=node,
                            count=taxon_data['count'],
                            rank=RANK_ORDER[taxon_data['rank']]
                            )
                        # Return to original parent node
                        parent_node = original_parent_node
                else:
                    node = anytree.Node(
                        taxon, parent=rootnode,
                        count=taxon_data['count'],
                        rank=RANK_ORDER[taxon_data['rank']]
                        )
                if isinstance(taxon_data['children'], dict):
                    itertaxa(taxon_data['children'], parent=taxon, parent_node=node)

        itertaxa(lineage_sample)
        trees[sample_id] = rootnode
    return trees


def check_counts(taxa_trees):
    """Verify that more general taxonomic ranks contains counts from more specific ones.

    Otherwise, add an 'Unclassified_<taxon>' Node to not miss those counts.

    :param taxa_trees (anytree.Node): Tree structure with a rootnode (the sample)
        and children nodes following the lineage structure.
    :return (dict): Tree structure with a rootnode (the sample)
        and children nodes following the lineage structure.
    """
    # TODO: This shouln't happen from bracken report, but allow this in future?.
    # Total counts per node == sum(counts of children nodes), i.e. [in counts]:
    # Gammaproteobacteria (c) = Enterobacterales (o) + found Order (o)
    # If not:
    # Gammaproteobacteria_unclassified (o) = Gammaproteobacteria (c) - sum(found order)

    for node in anytree.PreOrderIter(taxa_trees, maxlevel=taxa_trees.height):
        if sum([n.count for n in node.children]) != node.count:
            # Add unclassified node
            if not re.match('^Unclassified|^root', node.name):
                unclassified_taxon = f'Unclassified_{node.name}'
                counts = node.count - sum([n.count for n in node.children])
                anytree.Node(
                    unclassified_taxon, parent=node, rank=node.rank + 1, count=counts)
    return taxa_trees


def prepare_taxa_table(taxa_tree, rank):
    """Filter the abundance table for a specific rank.

    :param sample_taxa_tree (anytree.Node): Node with counts for all possible ranks.
    :param rank (str): Index of the taxonomic rank to use as filter.
    :return (DataFrame): DataFrame with counts per specific rank for one sample.
    """
    # Make df from anytree.Node
    d = []
    for node in anytree.search.findall(taxa_tree, lambda node: node.rank == rank):
        d.append(
            {
                'rank': node.rank,
                'count':  node.count,
                'tax': ';'.join([str(node.name) for node in node.path[1:]])
            }
        )
    # If the rank level does not exist, return None
    if d:
        df = pd.DataFrame(d)
        return df
    else:
        return None


def join_abundance_tables(taxa_trees, rank):
    """Take taxonomy trees and create a dataframe from them.

    :param taxa_trees (dict): Dictionary with one lineage tree per sample.
    :param rank (str): Taxonomic rank to subset the df.

    :return (DataFrame): Table with counts per specific rank for all the samples.
    """
    # Trees into tables
    tables = {}
    for s, t in taxa_trees.items():
        tables[s] = prepare_taxa_table(check_counts(t), RANK_ORDER[rank])
    df_all = pd.concat(list(tables.values()), keys=list(tables.keys()))
    # Remove multiindex
    df_all_samples = df_all.reset_index().rename(columns={'level_0': 'sample'})
    # Pivot table. Make each sample a different column
    df = df_all_samples.pivot_table(
        columns=['sample'], values='count', index=['tax'], fill_value=0)
    df['total'] = df.sum(axis=1)
    df = df.sort_values('total', ascending=False).reset_index()
    df.columns.name = None
    return df


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


def most_abundant_table(counts_per_taxa_df, samples, n=10, percent=False):
    """Select most abundant taxa and group the rest into 'Other' taxa.

    :param counts_per_taxa_df (DataFrame): DataFrame with counts per specific rank.
    :param samples (set): Set of samples.
    :param n (int, optional): Number of abundant taxa  to be chosen. Defaults to 10.
    :param percent (bool, optional): Make the percentage by column. Defaults to False.
    :return (DataFrame): DataFrame with most abundant taxa (in the total of all the\
            samples) and the rest group in Others.
    """
    # Extract n most abundant taxa & group less abundant in Others
    most_abundant_taxa = counts_per_taxa_df[:n]
    less_abundant_taxa = counts_per_taxa_df[n:].sum()[samples + ['total']].to_dict()
    less_abundant_taxa['tax'] = 'Other' + ';Other' * len(list(RANK_ORDER.keys())[1:])
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


def per_sample_stats(df_allstats):
    """Read from fastcat stats and return total reads per sample.

    :param df_allstats (DataFrame): Fastcat statistics.
    :return (DataFrame): Reads per sample.
    """
    df = df_allstats.groupby('sample_name').size().to_frame('count').reset_index()
    return df


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
