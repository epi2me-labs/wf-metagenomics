#!/usr/bin/env python
"""Create workflow report."""

import json
import argparse
import holoviews as hv
from bokeh.layouts import layout
from aplanat.components import fastcat
from aplanat.components import simple as scomponents
from aplanat.report import WFReport


hv.extension('bokeh')


def yield_sankey_items(entries, current, threshold=10):
    for i, j in entries.items():
        if j['count'] < threshold or j['rank'] in ['species']:
            continue
        yield(current, i, j['count'])
        for k in yield_sankey_items(j['children'], i, threshold):
            yield k


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
        "Workflow Ribosomal Survey Report", "wf-ribosomal-survey",
        revision=args.revision, commit=args.commit)

    sample_files = []
    for summ in args.summaries:
        sample_name = summ.split('.')[0]
        for lin in args.lineages:
            if lin.startswith(sample_name):
                sample_files.append([sample_name, summ, lin])
                break

    for sample in sample_files:
        section = report.add_section()
        with open(sample[2]) as lineage_file:
            data = json.load(lineage_file)

            entries = []
            for k, v in data.items():
                for entry in yield_sankey_items(v['children'], k):
                    entries.append(entry)

            sankey = hv.Sankey(entries)
            sankey.opts(height=900, width=1000)
            p = hv.render(sankey, backend='bokeh')

            section.markdown('## Sample: {}'.format(str(sample[0])))
            section.plot(
                layout([[p]], sizing_mode='scale_width')
            )

            section = report.add_section(
                section=fastcat.full_report(
                    [sample[1]],
                    header='### {} read stats'.format(str(sample[0]))
                ))

            section.markdown('---')

    report.add_section(
        section=scomponents.version_table(args.versions))
    report.add_section(
        section=scomponents.params_table(args.params))

    # write report
    report.write(args.report)


if __name__ == "__main__":
    main()
