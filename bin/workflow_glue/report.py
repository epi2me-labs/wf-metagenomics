#!/usr/bin/env python
"""Create workflow report."""
import json

from dominate.tags import em, p
import ezcharts as ezc
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import SeqSummary
from ezcharts.components.reports.labs import LabsReport
from ezcharts.layout.snippets import DataTable
from ezcharts.layout.snippets import Grid
from ezcharts.layout.snippets import Tabs
from ezcharts.plots import util
import pandas as pd
import workflow_glue.diversity as diversity
import workflow_glue.report_utils.report_utils as report_utils
from workflow_glue.report_utils.sankey import sankey_plot

from .util import get_named_logger, wf_parser  # noqa: ABS101

# Setup simple globals
WORKFLOW_NAME = 'wf-metagenomics'
REPORT_TITLE = f'{WORKFLOW_NAME}-report'
THEME = 'epi2melabs'
N_BARPLOT = 8  # number of taxa to plot in the barplot


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = LabsReport(
        "Workflow Metagenomics Sequencing Report", "wf-metagenomics",
        args.params, args.versions)

    #
    # 1. READ SUMMARY
    #
    if args.pipeline == 'minimap':
        if args.stats:
            stats_files = util.read_files(args.stats)
            with report.add_section("Read summary", "Read summary"):
                SeqSummary(stats_files)
                # Add metadata section
            with report.add_section("Samples summary", "Samples summary"):
                tabs = Tabs()
                with tabs.add_tab('Reads'):
                    with Grid(columns=1):
                        # Add barplot with reads per sample
                        sample_reads = ezc.barplot(data=report_utils.per_sample_stats(
                            stats_files), x='sample_name', y='count')
                        sample_reads.title = {"text": "Number of reads per sample."}
                        EZChart(sample_reads, THEME)
    # kraken stats are not in tsv, they came from json files with binned counts
    # to be plotted directly as an histogram
    else:
        with report.add_section("Read summary", "Read summary"):
            with open(args.stats[0]) as f:
                datas = json.load(f)
                tabs = Tabs()
                total_reads = {}
                for sample_id, data in sorted(datas.items()):
                    with tabs.add_tab(sample_id):
                        with Grid(columns=2):
                            EZChart(report_utils.read_quality_plot(data), THEME)
                            EZChart(report_utils.read_length_plot(data), THEME)
                            total_reads[sample_id] = data['total_reads']
                with tabs.add_tab('total'):
                    with Grid(columns=1):  # total read counts per sample
                        df_stats = pd.DataFrame.from_dict(total_reads.items())
                        df_stats.columns = ['Sample_name', 'Number of reads']
                        plt = ezc.barplot(
                            data=df_stats, x='Sample_name', y='Number of reads')
                        plt.title = {"text": "Number of reads per sample."}
                        plt.tooltip = {'trigger': 'axis'}
                        EZChart(plt, THEME)

    #
    # 2. TAXONOMY RESULTS
    #

    # Join taxonomy data
    all_json = report_utils.parse_lineages(args.lineages[0])
    # Extract all possible lineages
    allranks_tree = report_utils.tax_tree(all_json)
    samples = list(allranks_tree.keys())
    # Save all ranks info
    ranks_counts = []
    # 2.1. SANKEY
    sankey_plot(all_json, report)

    for rank in report_utils.RANKS_NO_SK_K:  # avoid superkingdom (SK), kingdom(K)
        counts_per_taxa_df = report_utils.join_abundance_tables(allranks_tree, rank)
        if not counts_per_taxa_df.empty:
            ranks_counts.append(counts_per_taxa_df)
    # Write report
    with report.add_section('Taxonomy', 'Taxonomy'):
        tabs = Tabs()
        # 2.2. BARPLOT
        with tabs.add_dropdown_menu('Rank', change_header=False):
            for i, counts_per_taxa_per_rank_df in enumerate(ranks_counts):
                with tabs.add_dropdown_tab(report_utils.RANKS_NO_SK_K[i]):
                    most_abundant = report_utils.most_abundant_table(
                        counts_per_taxa_per_rank_df, samples, n=N_BARPLOT, percent=True)
                    d2plot = report_utils.split_taxonomy_string(most_abundant)
                    # Long wide format
                    d2plot_melt = d2plot.melt(
                        id_vars=[report_utils.RANKS_NO_SK_K[i]], value_vars=samples,
                        var_name='samples', value_name='counts')
                    # Plot
                    p(f"Barplot of the {N_BARPLOT} most abundant taxa\
                    at the {report_utils.RANKS_NO_SK_K[i]} rank.\
                    Any remaining taxa have been collapsed under the \'Other\' category\
                    to facilitate the visualization.\
                    The y-axis indicates the relative abundance of each taxon\
                    in percentages\
                    for each sample.")
                    plt = ezc.barplot(
                        d2plot_melt, x='samples', y='counts',
                        hue=report_utils.RANKS_NO_SK_K[i], dodge=False)
                    plt.yAxis = dict(name='Relative abundance')
                    plt.legend = {
                        'orient': 'horizontal', 'left': 'center', 'top': 'bottom'}
                    plt.tooltip = {'trigger': 'axis', 'axisPointer': {'type': 'shadow'}}
                    # Not best option to avoid overlap, but echarts hasn't solved it yet
                    # https://github.com/apache/echarts/issues/15654
                    # https://github.com/apache/echarts/issues/14252
                    plt.grid = {'bottom': '18%'}
                    plt.title = {
                        "text": f"{report_utils.RANKS_NO_SK_K[i].capitalize()} rank"}
                    EZChart(plt, THEME)
        # 2.3. ABUNDANCE TABLE
    with report.add_section('Abundances', 'Abundances'):
        tabs = Tabs()
        with tabs.add_dropdown_menu('Abundance tables', change_header=False):
            for i, counts_per_taxa_per_rank_df in enumerate(ranks_counts):
                with tabs.add_dropdown_tab(report_utils.RANKS_NO_SK_K[i]):
                    p(f"Abundance table for the {report_utils.RANKS_NO_SK_K[i]} rank.")
                    export_table = report_utils.split_taxonomy_string(
                        counts_per_taxa_per_rank_df, set_index=True)
                    # Move tax column to end to not spoil visualization
                    temp_cols = export_table.columns.tolist()
                    new_cols = temp_cols[1:] + temp_cols[0:1]
                    export_table = export_table[new_cols]
                    # Table to report with the export option
                    DataTable.from_pandas(
                        export_table,
                        export=True, file_name=f'wf-metagenomics-counts-\
                            {report_utils.RANKS_NO_SK_K[i]}')
        # 2.4. RAREFIED ABUNDANCE TABLE
        with tabs.add_dropdown_menu('Rarefied Abundance tables', change_header=False):
            for i, counts_per_taxa_per_rank_df in enumerate(ranks_counts):
                with tabs.add_dropdown_tab(report_utils.RANKS_NO_SK_K[i]):
                    p(f"Rarefied abundance table for the \
                      {report_utils.RANKS_NO_SK_K[i]} rank.")
                    p("All samples have been randomly subsetted to have the same number\
                       of reads.")
                    # Table to report with the export option
                    rarefied_df = diversity.global_rarefaction(
                        counts_per_taxa_per_rank_df.set_index('tax')).sort_index(
                        axis=1).reset_index()
                    export_table = report_utils.split_taxonomy_string(
                        rarefied_df, set_index=True)
                    # Move tax column to end to not spoil visualization
                    temp_cols = export_table.columns.tolist()
                    new_cols = temp_cols[1:] + temp_cols[0:1]
                    export_table = export_table[new_cols]
                    DataTable.from_pandas(
                        export_table,
                        export=True, file_name=f'wf-metagenomics-rarefied-\
                        {report_utils.RANKS_NO_SK_K[i]}')

    #
    # 3. DIVERSITY
    #

    # Return counts for the last analyzed rank to calculate diversity
    last_analyzed_rank = ranks_counts[-1].set_index('tax')
    # Rarefy step by step
    rarefied_counts = last_analyzed_rank.apply(
        lambda x: diversity.rarefaction_curve(x), axis=0)
    richness_curve = {
        rarefied_counts.index[i]: rarefied_counts[i]
        for i in range(len(rarefied_counts.index))}

    with report.add_section("Alpha Diversity", "Diversity"):
        tabs = Tabs()
        #
        # 3.1. ALPHA DIVERSITY METRICS for the last analyzed level
        #
        with tabs.add_tab('Diversity indices'):
            p("Sample diversity indices. Indices are calculated from the original\
               abundance table.")
            DataTable.from_pandas(report_utils.calculate_diversity_metrics(
                last_analyzed_rank).set_index('Sample'))
            em("Note that the taxon 'Unkown' is considered as a unique taxon.")
        #
        # 3.2. SPECIES RICHNESS CURVES
        #
        with tabs.add_tab('Species richness curves'):
            p("Sample-based rarefaction curves to display observed taxa richness.\
                Sample size shows the number of reads sampled from the total amount\
              of reads analyzed during the real time analysis.\
              The Y-axis indicates the number of unique taxa at the last analyzed\
              taxonomic rank in those subsampled reads.\
            ")
            with Grid(columns=1):
                df_richness = pd.DataFrame.from_dict(
                    richness_curve, orient='index')
                df_richness['Sample'] = list(df_richness.index)
                df_richness_melt = df_richness.melt(
                    id_vars='Sample', value_vars=list(df_richness.columns),
                    var_name='Sample size', value_name='Richness')
                df_richness_melt_sort = df_richness_melt.sort_values(by=['Sample size'])
                plot = ezc.lineplot(
                    data=df_richness_melt_sort,
                    x='Sample size', y='Richness', hue='Sample')
                EZChart(plot, 'epi2melabs')
            em("Note that Unknown taxon is considered as a unique taxon.")

    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--stats", nargs='+', required=False,
        help="Read summary file.")
    parser.add_argument(
        "--lineages", nargs='+', required=True,
        help="Read lineage file.")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.")
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values")
    parser.add_argument(
        "--revision", default='unknown',
        help="git branch/tag of the executed workflow")
    parser.add_argument(
        "--commit", default='unknown',
        help="git commit of the executed workflow")
    parser.add_argument(
        "--pipeline", default='unknown',
        help="kraken or minimap")
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
