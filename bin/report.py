#!/usr/bin/env python
"""Create workflow report."""
import argparse
import json

from aplanat import bars
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport
from aplanat.util import Colors
from bokeh.layouts import layout
from jinja2 import Template
import natsort
import pandas as pd


def read_files(summaries, sep='\t'):
    """Read a set of files and join to single dataframe."""
    dfs = list()
    for fname in sorted(summaries):
        dfs.append(pd.read_csv(fname, sep=sep))
    return pd.concat(dfs)


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
    args = parser.parse_args()

    report = WFReport(
        "Workflow Metagenomics Report", "wf-metagenomics",
        revision=args.revision, commit=args.commit)

    templ = None
    with open(args.vistempl, "r") as vistempl:
        templ = vistempl.read()

    sample_lineages = {}
    for lineage in natsort.natsorted(args.lineages):
        lineage_name = lineage.split('.')[0]
        with open(lineage, 'r') as lf:
            sample_lineages[lineage_name] = json.load(lf)

    templ = templ.replace(
        "replace_me",
        json.dumps(sample_lineages).replace('"', '\\"'))
    report.template = Template(templ)

    #
    # Plot read counts per barcode
    #
    seq_summary = read_files(args.summaries)
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
        plot_width=None)
    bc_counts_plot.xaxis.major_label_orientation = 3.14/2

    section = report.add_section()
    section.markdown("### Samples")
    section.plot(layout([[bc_counts_plot]], sizing_mode="stretch_width"))

    #
    # Standard read metrics
    #
    for summ in args.summaries:
        section = report.add_section(
            section=fastcat.full_report(
                [summ],
                header='#### Read stats: {}'.format(str(summ.split('.')[0]))
            ))

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
