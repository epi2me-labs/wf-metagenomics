#!/usr/bin/env python
"""Create workflow report."""
import json
import os

from aplanat import bars, lines
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from aplanat.util import Colors
from bokeh.layouts import layout
from bokeh.plotting import figure
from bokeh.resources import INLINE
from jinja2 import Template
import pandas as pd
import workflow_glue.diversity as diversity

from .util import wf_parser  # noqa: ABS101


def read_files(summaries, sep='\t'):
    """Read a set of files and join to single dataframe."""
    dfs = list()
    for fname in sorted(summaries):
        dfs.append(pd.read_csv(fname, sep=sep))
    return pd.concat(dfs)


def plot_hist_data(counts, edges, **kwargs):
    """Create histogram from json."""
    defaults = {
        "output_backend": "webgl",
        "height": 300, "width": 600}
    defaults.update(**kwargs)
    p = figure(**defaults)
    p.quad(
            top=counts, bottom=0, left=edges[:-1], right=edges[1:],
            alpha=0.6)
    return p


def kraken(summaries, section):
    """Kraken quality report section."""
    with open(summaries) as f:
        datas = json.load(f)
    bc_counts = {}
    for sample_id, data in sorted(datas.items()):
        len_hist = list(data["len"].items())
        lbins, lcounts = list(zip(*len_hist))
        len_plot = plot_hist_data(
            lcounts, lbins, title="Read length distribution.",
            x_axis_label='Read Length / bases',
            y_axis_label='Number of reads'
        )
        qual_hist = list(data["qual"].items())
        edges, counts = list(zip(*qual_hist))
        total_reads = data["total_reads"]
        qual_plot = plot_hist_data(
            counts, edges, title="Read quality score.",
            x_axis_label="Quality score",
            y_axis_label="Number of reads")

        section.markdown("###"+sample_id)
        section.plot(
            layout([[len_plot, qual_plot]], sizing_mode="stretch_width"))
        bc_counts[sample_id] = total_reads

    bc_counts_plot = bars.simple_bar(
        list(bc_counts.keys()),
        list(bc_counts.values()),
        colors=[Colors.cerulean]*len(bc_counts),
        title='Number of reads per sample.',
        plot_width=800)
    bc_counts_plot.xaxis.major_label_orientation = 3.14/2
    section.markdown("### Samples")
    section.plot(layout([[bc_counts_plot]], sizing_mode="stretch_width"))


def minimap(summary, report):
    """Create minimap quality section."""
    seq_summary = read_files(summary)
    bc_counts = (
        pd.DataFrame(
            seq_summary['sample_name'].value_counts()
        ).sort_index().reset_index().rename(
            columns={'index': 'sample', 'sample_name': 'count'})
    )

    bc_counts_plot = bars.simple_bar(
        bc_counts['sample'].astype(str),
        bc_counts['count'],
        colors=[Colors.cerulean]*len(bc_counts),
        title='Number of reads per sample.',
        plot_width=800)
    bc_counts_plot.xaxis.major_label_orientation = 3.14/2
    section = report.add_section()
    section.markdown("### Samples")
    section.plot(layout([[bc_counts_plot]], sizing_mode="stretch_width"))
    #
    # Standard read metrics
    #
    stats_summaries = dict(tuple(seq_summary.groupby(['sample_name'])))

    for key, summ in stats_summaries.items():
        summ.to_csv(
            'file_name.csv', mode='w+', index=False, header=True, sep="\t")
        section = report.add_section(
            section=fastcat.full_report(
                'file_name.csv',
                header='#### Read stats: {}'.format(str(key))
            ))


def collector_curves(rarefaction, section):
    """Plot observed species richness curves.

    :param rarefaction: dictionary with richness
        per number of sampled reads.
    :type rarefaction: dict.
    :param section: new section to add the plots.
    """
    # Also called Species accumulation curves
    for sample_id, data in sorted(rarefaction.items()):
        rarefaction_curve = lines.line(
            x_datas=[list(map(float, data.keys()))],
            y_datas=[list(map(float, data.values()))],
            names=[sample_id],
            title="Observed Taxa Richness.",
            x_axis_label='Sample size',
            y_axis_label='S (Number of unique taxa)'
        )
        section.markdown("### " + sample_id)
        section.plot(
            layout([[rarefaction_curve]], sizing_mode="stretch_width"))


def main(args):
    """Run the entry point."""
    report = WFReport(
        "Workflow Metagenomics Report", "wf-metagenomics",
        revision=args.revision, commit=args.commit)

    templ = None
    with open(args.vistempl, "r") as vistempl:
        templ = vistempl.read()

    #
    # Sankey plot
    #
    path_to_json = args.lineages[0]
    json_files = [
        pos_json for pos_json in os.listdir(
            path_to_json) if pos_json.endswith('.json')]
    all_json = {}
    all_json_rank = {}
    all_rarefied_counts = {}
    for i in json_files:
        with open(os.path.join(args.lineages[0], i)) as json_file:
            # process sankey and table counts
            sample_lineages = json.load(json_file)
            sample_name = list(sample_lineages.keys())[0]
            all_json.update(sample_lineages)
            # process counts per specific taxa rank to rarefy it
            rank_counts = diversity.extract_rank(
                list(sample_lineages.values())[0],
                rank_dict={}, rank=args.rank)
            # rarefy step by step
            species_richness = {sample_name: diversity.rarefaction_curve(
                list(rank_counts.values()))}
            all_rarefied_counts.update(species_richness)
            # join all count data in a dict to then do the final rarefaction
            all_json_rank[sample_name] = rank_counts

    templ = templ.replace(
        "replace_me",
        json.dumps(all_json).replace('"', '\\"'))
    report.template = Template(templ)
    bokeh_resources = INLINE.render()
    report.template.render(bokeh_resources=bokeh_resources)
    #
    # Standard read metrics
    #
    section = report.add_section()
    if args.pipeline == "minimap":
        minimap(args.summaries, report)
    else:
        kraken(args.summaries[0], section)
    #
    # Plot richness with different sample sizes.
    #
    section = report.add_section()
    section.markdown(" # Diversity ")
    section.markdown(" Sample-based rarefaction curves to \
                    display observed taxa richness. ")
    section.markdown(" Sample size shows the number of reads \
                        sampled from the total amount of reads analyzed \
                        during the real time analysis. \
                        Y-axis indicates the number of \
                        unique taxa (S) found in those subsampled reads. ")
    section.markdown(" *Note that Unknown taxon is \
                    considered as an unique taxon.* ")
    collector_curves(all_rarefied_counts, section)
    #
    # Provide alpha-diversity indexes
    #
    # Global rarefaction
    df = pd.DataFrame.from_dict(all_json_rank, orient='columns').sort_index(axis=1)
    df = df.fillna(0)
    df.to_csv('wf-metagenomics-counts.tsv', sep="\t")
    rarefied_df = diversity.global_rarefaction(df)
    # rarefied_df = rarefied_df.sort_index(axis=1)
    rarefied_df.to_csv('wf-metagenomics-rarefied.tsv', sep="\t")
    div = diversity.alpha_diversity(df).round(2).sort_index().reset_index(
        level=0)
    div.fillna('None', inplace=True)
    div.rename(columns={'index': 'Sample'}, inplace=True)
    # div.to_csv('wf-metagenomics-alpha.tsv', sep="\t")
    section = report.add_section()
    section.markdown(" # Alpha diversity ")
    section.markdown('''
        Indexes are calculated from \
    the abundance table \
    (see description in:\
    Provided outputs)
    ''')
    section.table(div)
    section = report.add_section()
    section.markdown(" ### Provided outputs ")
    d_output = {
        'wf-metagenomics-report.html': ['this report', 'html', ''],
        'wf-metagenomics-counts.tsv': [
            'counts per taxa', 'tsv',
            'cols: samples, rows: taxa at a specific rank'
        ],
        'wf-metagenomics-rarefied.tsv': [
            'rarefied counts per sample', 'tsv',
            'cols: samples, rows: taxa at a specific rank'
        ]
    }
    df_output = pd.DataFrame.from_dict(d_output).T
    df_output.columns = ['Description', 'Format', 'Notes']
    df_output.reset_index(level=0, inplace=True)
    df_output.rename(columns={'index': 'Output'}, inplace=True)
    section.markdown("This is a description of the output directory")
    section.table(df_output)
    #
    # Standard wf reporting
    #

    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))
    report.write(args.report)


def argparser():
    """Argument parser for entrypoint."""
    parser = wf_parser("report")
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--summaries", nargs='+', required=True,
        help="Read summary file.")
    parser.add_argument(
        "--lineages", nargs='+', required=True,
        help="Read lineage file.")
    parser.add_argument(
        "--rank", default='S',
        help="Taxonomic rank to run divresity.")
    parser.add_argument(
        "--vistempl", required=True)
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
