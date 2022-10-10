#!/usr/bin/env python
"""Create workflow report."""
import argparse
import json
import os

from aplanat import bars
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from aplanat.util import Colors
from bokeh.layouts import layout
from bokeh.plotting import figure
from bokeh.resources import INLINE
from jinja2 import Template
import pandas as pd


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
    for sample_id, data in datas.items():
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
            counts, edges, title="Read quality score",
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
        title='Number of reads per sample',
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
        title='Number of reads per sample',
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


def main():
    """Run the entry point."""
    parser = argparse.ArgumentParser()
    parser.add_argument("report", help="Report output file")
    parser.add_argument(
        "--summaries", nargs='+', required=True,
        help="Read summary file.")
    parser.add_argument(
        "--lineages", nargs='+', required=True,
        help="Read lineage file.")
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
    args = parser.parse_args()

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
    for i in json_files:
        with open(os.path.join(args.lineages[0], i)) as json_file:
            sample_lineages = json.load(json_file)
            all_json.update(sample_lineages)
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
    # Standard wf reporting
    #
    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    report.write(args.report)


if __name__ == "__main__":
    main()
