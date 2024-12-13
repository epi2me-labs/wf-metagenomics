#!/usr/bin/env python
"""Create workflow report."""

import json

from bokeh.models import HoverTool
from dominate import tags as html_tags
from dominate.tags import em, p
import ezcharts as ezc
from ezcharts.components.ezchart import EZChart
from ezcharts.components.fastcat import SeqSummary
from ezcharts.components.reports.labs import LabsReport
from ezcharts.layout.snippets import DataTable
from ezcharts.layout.snippets import Grid
from ezcharts.layout.snippets import Tabs
import pandas as pd
import workflow_glue.diversity as diversity
import workflow_glue.report_utils.report_utils as report_utils


from .util import get_named_logger, wf_parser  # noqa: ABS101

# Setup simple globals
THEME = 'epi2melabs'


def amr_section(amr_data, html_id):
    """Parse amr JSON for accordion style table.

    params: amr_data (dict): Dict containing amr results for sample.
    params: html_id (str): String used for grouping data correctly in html.
    returns (str): html format of data.
    """
    _div = html_tags.div(cls="accordion-item")
    for i, (gene, data) in enumerate(amr_data.items()):
        _head = html_tags.h2(id=str(i), style="border: 1px solid rgba(0,0,0,.125);\
                            border-collapse: collapse;\
                            padding:0;\
                            margin-bottom:0")
        _button = html_tags.button(
            html_tags.span(html_tags.b(gene)),
            html_tags.span(
                    data["count"],
                    cls="badge rounded-pill bg-danger",
                    style="top:8px"
                    ),
            cls="accordion-button collapsed",
            type="button",
            data_bs_toggle="collapse",
            data_bs_target=f"#collapse{i}",
            aria_expanded="false",
            aria_controls=f"collapse{i}",
            style="display: grid; \
                    align-items: center;\
                    grid-template-columns: 1fr max-content max-content;\
                    grid-gap: 25px"
            )
        _head.add(_button)
        _div.add(_head)
        _div1 = html_tags.div(
            id=f"collapse{i}",
            cls="accordion-collapse collapse",
            fr=str(i),
            aria_labelledby=str(i),
            data_bs_parent=f"#{html_id}")
        _div2 = html_tags.div(cls="accordion body")
        _table = html_tags.table(cls="table table-striped")
        _thead = html_tags.thead()
        _thead.add(
            html_tags.tr(
                html_tags.th("ReadID"),
                html_tags.th("Coverage %"),
                html_tags.th("Identity %"),
                html_tags.th("Resistance")
            )
        )
        _table.add(_thead)
        for hit in data["meta"]:
            _tr = html_tags.tr()
            _tr.add(
                html_tags.td(hit["SEQUENCE"]),
                html_tags.td(hit["%COVERAGE"]),
                html_tags.td(hit["%IDENTITY"]),
                html_tags.td(hit["RESISTANCE"].replace("_", ""))
            )
            _table.add(_tr)
        _div2.add(_table)
        _div1.add(_div2)
        _div.add(_div1)
    return _div


def main(args):
    """Run the entry point."""
    logger = get_named_logger("Report")
    report = LabsReport(
        f"{args.workflow_name} Sequencing Report", args.workflow_name,
        args.params, args.versions, args.wf_version)
    global SELECTED_RANKS
    SELECTED_RANKS = []
    for i in report_utils.RANKS:
        # append each rank
        SELECTED_RANKS.append(i)
        if i == report_utils.RANKS_ABB[args.taxonomic_rank]:
            # stop when the user rank has been added
            break
    ranks_no_sk_k = SELECTED_RANKS[2:]
    abundance_threshold = args.abundance_threshold
    n_taxa_barplot = args.n_taxa_barplot

    #
    # 1. READ SUMMARY
    #
    if args.pipeline != 'real_time':
        # Samples
        with open(args.metadata) as metadata:
            sample_details = [{
                'sample': d['alias'],
                'type': d['type'],
                'barcode': d['barcode']
            } for d in json.load(metadata)]
        if args.read_stats:
            with report.add_section("Read summary", "Read summary"):
                names = tuple(d['sample'] for d in sample_details)
                stats = tuple(args.read_stats)
                if len(stats) == 1:
                    stats = stats[0]
                    names = names[0]
                SeqSummary(stats, sample_names=names)

    # kraken stats for the real time are not in tsv,
    # they came from json files with binned counts
    # to be plotted directly as an histogram
    else:
        with report.add_section("Read summary", "Read summary"):
            with open(args.read_stats[0]) as f:
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
                        plt._fig.title.text = "Number of reads per sample."
                        plt._fig.xaxis.major_label_orientation = 45
                        hover = plt._fig.select(dict(type=HoverTool))
                        hover.tooltips = [("Number of reads", "@top")]
                        EZChart(plt, THEME)
    if args.pipeline != 'real_time':
        with report.add_section("Number of reads", "Reads"):
            p("""
                Number of reads after applying read length and quality filters.
                Read counts will also reflect host depletion and unclassified
                reads (kraken2 approach) and unmapped reads (minimap2 approach).
                Percentages are calculated from reads after the filtering.
            """)
            p("""
                Note that in the minimap2 approach, there are two filters
                based on the mapping identity and coverage that can increase
                the number of unclassified sequences. This table shows the
                unmapped reads before both filters are applied.
              """)
            DataTable.from_pandas(report_utils.n_reads_pass(args.metadata))

    #
    # 2. TAXONOMY RESULTS
    #

    # Join taxonomy data
    all_json = report_utils.parse_lineages(args.lineages[0])
    abundance_table = pd.read_csv(args.abundance_table, sep='\t')
    # need to read samples from the table, as it is in a different process
    # is not updated the same time that lineages: lineages can include 2 samples and
    # at the same moment abundance table just contains 1 of them.
    samples = [i for i in abundance_table.columns if i not in ['tax', 'total']]
    # 2.1. SANKEY
    with report.add_section('Lineages', 'Lineages'):
        ezc.metagenomics_sankey(all_json)
    # 2.2. SUNBURST
    with report.add_section('Sunburst', 'Sunburst'):
        tabs = Tabs()
        with tabs.add_dropdown_menu('Sample', change_header=True):
            for barcode in sorted(all_json.keys()):
                logger.info(f"Sample {barcode}.")
                with tabs.add_dropdown_tab(barcode):
                    p("""
                    This visualization can be useful to interactively explore the
                    taxonomic composition of the sample. It shows hierarchical data:
                    each layer represents a taxonomic rank. The color indicates the
                    abundance of that taxon in comparison with the total number of
                    reads. Zooming into an specific taxon is possible by clicking on it.
                    """)
                    lineages = report_utils.prepare_data_to_sunburst(all_json[barcode])
                    plt = ezc.sunburst(
                        lineages, label_rotate="tangential",
                        label_minAngle=25)
                    EZChart(plt, THEME)
                    lineages.clear()
    # Save all ranks info
    ranks_counts_filtered = []
    for rank in ranks_no_sk_k:  # avoid superkingdom (SK), kingdom(K)
        abundance_table_rank = abundance_table.copy()
        # update full taxonomy string to an string which finishes at the rank asked.
        abundance_table_rank['tax'] = [';'.join(
            i.split(';')[:report_utils.RANK_ORDER[rank]]
            ) for i in abundance_table['tax']]
        counts_per_taxa_df = abundance_table_rank.groupby(['tax']).sum().reset_index()
        if not counts_per_taxa_df.empty:
            # Filter by abundance threshold.
            # Distinguish between natural number cutoff or percentage.
            counts_per_taxa_df_filtered = report_utils.filter_by_abundance(
                counts_per_taxa_df, 'total', abundance_threshold)
            ranks_counts_filtered.append(counts_per_taxa_df_filtered)
            ranks_counts = counts_per_taxa_df  # save last table
    # Write report
    with report.add_section('Taxonomy', 'Taxonomy'):
        tabs = Tabs()
        # 2.3. BARPLOT
        with tabs.add_dropdown_menu('Rank', change_header=True):
            for i, counts_per_taxa_per_rank_df in enumerate(ranks_counts_filtered):
                with tabs.add_dropdown_tab(ranks_no_sk_k[i]):
                    logger.info(f"rank {ranks_no_sk_k[i]}.")
                    most_abundant = report_utils.most_abundant_table(
                        counts_per_taxa_per_rank_df,
                        n=n_taxa_barplot,
                        percent=True)
                    d2plot = report_utils.split_taxonomy_string(most_abundant)
                    # Long wide format
                    d2plot_melt = d2plot.melt(
                        id_vars=[ranks_no_sk_k[i]], value_vars=samples,
                        var_name='samples', value_name='counts')
                    # Plot
                    p(f"""
                    Barplot of the {n_taxa_barplot} most abundant taxa at the
                    {ranks_no_sk_k[i]} rank in all the samples.
                    Any remaining taxa have been collapsed under the \'Other\' category
                    to facilitate the visualization.
                    The y-axis indicates the relative abundance of each taxon
                    in percentages for each sample.
                    """)
                    plt = ezc.barplot(
                        d2plot_melt,
                        x="samples",
                        y="counts",
                        hue=ranks_no_sk_k[i],
                        dodge=False,
                    )
                    plt._fig.yaxis.axis_label = "Relative abundance"
                    plt._fig.xaxis.major_label_orientation = 45
                    # distribute names in 5 columns as per default 10 taxa are shown.
                    # > 10 is difficult to distinguish colors.
                    plt._fig.legend.ncols = 5
                    # mute allows to soften the color of the chosen taxa in the legend.
                    plt._fig.legend.click_policy = "mute"
                    plt._fig.title.text = f"{ranks_no_sk_k[i].capitalize()} rank"
                    plt._fig.title.align = "center"
                    hover = plt._fig.select(dict(type=HoverTool))
                    hover.tooltips = [("Taxon", "$name"), ("Percent", "@$name %")]
                    EZChart(plt, THEME)
        # 2.4. ABUNDANCE TABLE
    with report.add_section('Abundances', 'Abundances'):
        tabs = Tabs()
        with tabs.add_dropdown_menu('Abundance tables', change_header=False):
            for i, counts_per_taxa_per_rank_df in enumerate(ranks_counts_filtered):
                with tabs.add_dropdown_tab(ranks_no_sk_k[i]):
                    p(f"Abundance table for the {ranks_no_sk_k[i]} rank.")
                    p(
                        "Only taxa whose global abundance are above the ",
                        html_tags.code("abundance_threshold"),
                        " parameter will appear in the table."
                    )
                    export_table = report_utils.split_taxonomy_string(
                        counts_per_taxa_per_rank_df, set_index=True)
                    # Move tax column to end to not spoil visualization
                    temp_cols = export_table.columns.tolist()
                    new_cols = temp_cols[1:] + temp_cols[0:1]
                    export_table = export_table[new_cols]
                    # Table to report with the export option
                    DataTable.from_pandas(
                        export_table,
                        export=True,
                        file_name=(
                            f'{args.workflow_name}-counts-{ranks_no_sk_k[i]}'
                        )
                    )
        # 2.5. RAREFIED ABUNDANCE TABLE
        with tabs.add_dropdown_menu('Rarefied Abundance tables', change_header=False):
            for i, counts_per_taxa_per_rank_df in enumerate(ranks_counts_filtered):
                with tabs.add_dropdown_tab(ranks_no_sk_k[i]):
                    p(f"Rarefied abundance table for the \
                      {ranks_no_sk_k[i]} rank.")
                    p("""
                        All samples have been randomly subsetted to have the same number
                        of reads.
                    """)
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
                        export=True,
                        file_name=(
                            f'{args.workflow_name}-rarefied-{ranks_no_sk_k[i]}'
                        )
                    )

    #
    # 3. DIVERSITY
    #

    # Return counts for the last analyzed rank to calculate diversity
    last_analyzed_rank = ranks_counts[samples + ['total']]
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
            p("""
            Sample diversity indices. Indices are calculated from the original
            abundance table.
            """)
            diversity_df = report_utils.calculate_diversity_metrics(
                last_analyzed_rank)
            diversity_df.rename(columns={'Sample': 'Indices'}, inplace=True)
            DataTable.from_pandas(
                diversity_df.set_index('Indices'), export=True,
                file_name=(f'{args.workflow_name}-diversity'))
            em("Note that the taxon 'Unknown' is considered as a unique taxon.\n")

        #
        # 3.2. SPECIES RICHNESS CURVES
        #
        with tabs.add_tab('Species richness curves'):
            p("""Sample-based rarefaction curves to display observed taxa richness.
                Sample size shows the number of reads sampled from the total amount
              of reads analyzed during the real time analysis.
              The Y-axis indicates the number of unique taxa at the last analyzed
              taxonomic rank in those subsampled reads.
            """)
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

        #
        # 3.3. SPECIES ABUNDANCE DISTRIBUTION: SAD
        #
        with tabs.add_dropdown_menu('Taxa abundance distribution', change_header=False):
            for barcode in samples:
                with tabs.add_dropdown_tab(barcode):
                    # Remove Unknown
                    ranks_counts = ranks_counts[
                        ~(ranks_counts['tax'].str.contains('Unknown'))]
                    df_sample_counts = (
                        ranks_counts[barcode]
                        .sort_values(ascending=False)
                        .rename("freq")
                        .to_frame()
                    )
                    df_sample_counts.columns = ["freq"]
                    if df_sample_counts["freq"].sum() > 0:
                        df_sample_counts.index = list(
                            range(1, df_sample_counts.shape[0] + 1)
                        )
                        plt = ezc.barplot(
                            df_sample_counts.reset_index().rename(
                                columns={"index": "rank"}
                            ),
                            x="rank",
                            y="freq",
                            dodge=False,
                            color="#0079a4",
                        )
                        plt._fig.title.text = barcode
                        plt._fig.title.align = "center"
                        plt._fig.xaxis.axis_label = "Rank"
                        plt._fig.xaxis.major_label_text_color = None
                        plt._fig.xaxis.major_tick_line_color = None
                        EZChart(plt, 'epi2melabs')
                        em("""This plot includes all the counts (except Unknown),
                        before applying any filter threshold based on abundances.
                        """)
                    # case of barcode03 with all unclassified.
                    else:
                        em("""There are no taxa to display.""")

    #
    # 4. AMR
    #
    if args.amr:
        with report.add_section("Antimicrobial resistance", "AMR"):
            with open(args.params) as f:
                params = json.load(f)
            amr_db = params["amr_db"].capitalize()
            p(f"""Detection of acquired AMR genes within sample using Abricate
                with the {amr_db} database.
            Please note that SNP-mediated AMR cannot be detected.
            """)
            amr_data = report_utils.parse_amr(args.amr)
            tabs = Tabs()
            for sample_id, data in sorted(amr_data.items()):
                with tabs.add_tab(sample_id):
                    with html_tags.div(cls="accordion", id=f"accordian_{sample_id}"):
                        if data["pass"]:
                            amr_section(data["results"], f"accordion_{sample_id}")
                        else:
                            p("No AMR genes detected in sample")

    #
    # 5. ALIGNMENT STATS
    #
    if args.align_stats:
        heatmap_min_cov = 1
        samples_references = {}
        for s in samples:
            samples_references[s] = report_utils.load_alignment_data(
                args.align_stats,
                s,
                args.taxonomic_rank,
                heatmap_min_cov=heatmap_min_cov,
            )
        # make sure that the samples really have data.
        dataset_results = {k: v for k, v in samples_references.items() if v is not None}
        if len(dataset_results) >= 1:
            with report.add_section("Alignment Statistics", "References"):
                tabs = Tabs()
                # Show table with the detailed output.
                with tabs.add_dropdown_menu("Dataset"):
                    # It will apply on the last analyzed rank.
                    taxa = [  # split the taxonomy string to the last value
                        i.split(';')[-1] for i in counts_per_taxa_df_filtered['tax']]
                    for barcode, metrics in dataset_results.items():
                        # Choose those species above the abundance threshold.
                        align_stats = metrics[0]
                        align_stats_filtered = align_stats[align_stats[rank].isin(taxa)]
                        with tabs.add_dropdown_tab(barcode):
                            p(
                                "Only taxa present in the abundance table above the ",
                                html_tags.code("abundance_threshold"),
                                " parameter appear in the table."
                            )
                            DataTable.from_pandas(
                                    align_stats_filtered,
                                    export=True,
                                    file_name='wf-metagenomics-alignment'
                            )
                # Show reference scatterplot of number of reads by coverage.
                with tabs.add_dropdown_menu("Scatter", change_header=False):
                    for barcode, metrics in dataset_results.items():
                        with tabs.add_dropdown_tab(barcode):
                            EZChart(metrics[1], 'epi2melabs')
                # Show heatmap of relative positional coverage.
                with tabs.add_dropdown_menu("Heatmap", change_header=False):
                    for barcode, metrics in dataset_results.items():
                        with tabs.add_dropdown_tab(barcode):
                            if metrics[2]:
                                p(
                                    "To illustrate consistency of coverage between ",
                                    "the reference sequences, each reference is ",
                                    "divided into 100 evenly sized windows, and the ",
                                    "average depth across all positions in the window ",
                                    "is plotted in a cell in the heatmap. Only ",
                                    f"references with {heatmap_min_cov}% average ",
                                    "coverage across the entire sequence are included ",
                                    "in the heatmap."
                                )
                                EZChart(metrics[2], 'epi2melabs')
                            else:
                                p(
                                    "No taxa present with sufficient coverage for heatmap."  # noqa:E501
                                )
    report.write(args.report)
    logger.info(f"Report written to {args.report}.")


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--workflow_name", required=True,
        help="The name of the workflow.")
    parser.add_argument(
        "--metadata", default='metadata.json', required=False,
        help="sample metadata")
    parser.add_argument(
        "--read_stats",  nargs='+', required=False,
        help="Fastcat per-read stats, ordered as per entries in --metadata",
        type=report_utils.is_not_empty_or_exit)
    parser.add_argument(
        "--lineages", nargs='+', required=True,
        help="Read lineage file.",
        type=report_utils.is_not_empty_or_exit)
    parser.add_argument(
        "--align_stats", required=False,
        help="Folder containing the mapping and depth statistics in TSV format.",
        type=report_utils.is_not_empty_or_exit)
    parser.add_argument(
        "--abundance_table", required=True,
        help="Read abundance tsv file.",
        type=report_utils.is_not_empty_or_exit)
    parser.add_argument(
        '--taxonomic_rank', required=True, choices=["S", "G", "k", "F", "O", "C", "P"],
        help="Taxonomic rank.")
    parser.add_argument(
        "--versions", required=True,
        help="directory containing CSVs containing name,version.",
        type=report_utils.is_not_empty_or_exit)
    parser.add_argument(
        "--params", default=None, required=True,
        help="A JSON file containing the workflow parameter key/values",
        type=report_utils.is_not_empty_or_exit)
    parser.add_argument(
        "--pipeline", default='kraken2', choices=["kraken2", "minimap2", "real_time"],
        help="kraken2, minimap2 or real_time")
    parser.add_argument(
        "--amr", default=None,
        help="Path to combined AMR results",
        type=report_utils.is_not_empty_or_exit)
    parser.add_argument(
        "--abundance_threshold", default=1, type=float,
        help="Remove those taxa whose abundance is below this cut-off.")
    parser.add_argument(
        "--n_taxa_barplot", default=9, type=int,
        help="Number of taxa to be displayed in the barplot.")
    parser.add_argument(
        "--wf_version", default='unknown',
        help="version of the executed workflow")
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
