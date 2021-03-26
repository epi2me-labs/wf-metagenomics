#!/usr/bin/env python3
"""
Convert taxids to names.

positional arguments:
  master_table  Master table from generate_master_table.py

optional arguments:
  -h, --help    show this help message and exit
"""
import argparse
import pathlib

from aplanat import annot, hist
from aplanat.components import fastcat
import aplanat.graphics
import aplanat.report
import ete3
import pandas as pd


parser = argparse.ArgumentParser(description="Convert taxids to names.")
parser.add_argument(
    "master_table", help="Master table from generate_master_table.py"
)
parser.add_argument(
    "fastcat_summary", help="fastcat sequencing summary stats"
)

report = aplanat.report.HTMLReport(
    "Metagenomics report",
    "Results generated through the wf-metagenomics nextflow workflow provided "
    "by Oxford Nanopore Technologies.",
)
exec_summary = aplanat.graphics.InfoGraphItems()

NCBI = ete3.NCBITaxa()


def _find_name(taxid):
    """Find scientific name from taxid using ete3."""
    if taxid == -1:
        return "Unclassified"
    try:
        return NCBI.get_taxid_translator([taxid])[taxid]
    except ValueError:
        return f"Unknown taxid: {taxid}"


def get_read_stats(master_table):
    """Generate read stats from master table."""
    read_stats = {
        "total_bases": 0,
        "mean_quality": 0,
        "n50_length": 0,
        "mean_length": 0,
    }
    if master_table.empty:
        return read_stats

    sorted_lengths = master_table.len.sort_values(ascending=False).reset_index(
        drop=True
    )
    cumulative_length = sorted_lengths.cumsum()
    total_bases = cumulative_length.iloc[-1]
    read_stats["total_bases"] = total_bases
    mean_quality = round(master_table.meanqscore.mean(), 2)
    read_stats["mean_quality"] = mean_quality

    # read length stats
    n50_index = cumulative_length.searchsorted(total_bases / 2)
    n50_length = int(sorted_lengths.iloc[n50_index])
    read_stats["n50_length"] = n50_length
    mean_length = round(total_bases / len(master_table), 2)
    read_stats["mean_length"] = mean_length

    return read_stats


def main(argv=None):
    """Process table and generate report.html."""
    args = parser.parse_args(argv)
    master_table = pathlib.Path(args.master_table)
    master_table = pd.read_csv(master_table)
    fastcat_stats = args.fastcat_summary

    pd.options.display.float_format = "{:.1f}".format

    # Classification numbers
    unclassified_reads = master_table[master_table["label"] == "unclassified"]
    classified_reads = master_table[master_table["label"] != "unclassified"]
    exec_summary.append(
        "Unclassified reads", len(unclassified_reads), "chart-pie"
    )
    exec_summary.append("Classified reads", len(classified_reads), "chart-pie")
    exec_summary.append("Total reads", len(master_table), "calculator")

    # Find genera counts
    df = pd.value_counts(master_table["genus"]).to_frame().reset_index()
    df.columns = ["taxid", "count"]
    df["scientific_name"] = df.apply(lambda row: _find_name(row.taxid), axis=1)
    is_classified = df["taxid"] != -1
    report.markdown("### Genera overview", key="genera-head")
    report.markdown(
        "The following shows the top 10 genera found in the sample sorted by "
        "number of reads classified.  Please note: only organisms present in "
        "the original database can be discovered in the sample!",
        key="genera-desc",
    )

    exec_summary.append(
        "Genera found", str(int(len(df[is_classified]))), "project-diagram"
    )
    # print(df[is_classified])  # genera counts (excluding unclassified)
    top_gen = (
        df[is_classified].sort_values(by="count", ascending=False).head(10)
    )
    top_gen = top_gen.reset_index()
    top_gen.drop(top_gen.columns[0], axis=1, inplace=True)
    report.table(
        top_gen, key="Top 10 Genera detected by read count", index=False
    )

    # Fastq numbers
    report.markdown("### Sequence overview", key="fq-head")
    report.markdown(
        "Read length and quality impact the success of classification.  "
        "Below is a summary of the data.",
        key="fq-desc",
    )
    read_stats = get_read_stats(master_table)
    exec_summary.append("Total bases", read_stats["total_bases"], "dna")
    exec_summary.append("Mean read qscore", read_stats["mean_quality"], "gem")
    exec_summary.append(
        "Read length N50", read_stats["n50_length"], "ruler-horizontal"
    )
    exec_summary.append(
        "Mean length", read_stats["mean_length"], "ruler-horizontal"
    )
    fq_summary = pd.DataFrame.from_dict(
        read_stats, orient="index", columns=["Value"]
    )
    report.table(fq_summary, key="Fastq summary stats")


    #------
    report.add_section(section=fastcat.full_report(fastcat_stats))
    #------


    # # fastq graphs
    # datas = [master_table.len]
    # xlim_max = read_stats["mean_length"] + (master_table["len"].std() * 6)
    # plot = hist.histogram(
    #     datas, bins=400, title="Read length distribution", xlim=[0, xlim_max]
    # )
    # # add vertical lines for mean and N50 read length
    # annot.marker_vline(
    #     plot,
    #     read_stats["mean_length"],
    #     "Mean: {:.0f}".format(read_stats["mean_length"]),
    # )
    # annot.marker_vline(
    #     plot,
    #     read_stats["n50_length"],
    #     "N50: {}".format(read_stats["n50_length"]),
    #     text_baseline="top",
    # )
    # # add axis labels
    # plot.xaxis.axis_label = "Read Length / bases"
    # plot.yaxis.axis_label = "Number of reads"
    # report.plot(plot, key="length-plain")

    report.markdown(
        "### Classified vs unclassified reads", key="classification-debug-head"
    )
    report.markdown(
        "When considering the validity of classifications the quality of the"
        " data and composition of the original database are crucial.  Some"
        " unclassified reads in a classification experiment are expected, the"
        " proportion of unclassified reads will depend on the experiment.  One"
        " way to determine if there are organisms present that are not"
        " classified due to lack of representation in the original database is"
        " to compare the distribution of length and quality of reads that are"
        " unclassified vs classified.  If there are many unclassified reads"
        " but there is not a significant difference between the statistics of"
        " the population of classified and unclassified reads this indicates"
        " that the organism may not be adequately represented in the original"
        " database.  Consequently, as poorer quality data is less likely to"
        " classify as anything specific, adjusting prep to improve sequencing"
        " would be a good avenue to pursue.",
        key="classification-debug-desc",
    )

    report.markdown(
        "The following table shows the key summary statistics of the data "
        "split by classification status",
        key="classification-debug-table-desc",
    )

    # Read stats classified vs unclassified
    read_stats_classified = get_read_stats(classified_reads)
    read_stats_unclassified = get_read_stats(unclassified_reads)
    df_class = pd.DataFrame.from_dict(
        read_stats_classified, orient="index", columns=["Value"]
    )
    df_unclass = pd.DataFrame.from_dict(
        read_stats_unclassified, orient="index", columns=["Value"]
    )
    df_cvu_stats = df_class.join(
        df_unclass, lsuffix="_classified", rsuffix="_unclassified"
    )
    df_cvu_stats.rename(
        columns={
            "Value_classified": "Classified",
            "Value_unclassified": "Unclassified",
        },
        inplace=True,
    )
    report.table(
        df_cvu_stats,
        key="Key read stats for classified and unclassified reads",
    )

    # Executive summary
    report.markdown("## Executive summary", key="exec-head")
    report.markdown(
        "The following summarises the key findings of this workflow.",
        key="exec-desc",
    )
    report.plot(None, "exec-plot")
    exec_plot = aplanat.graphics.infographic(exec_summary.values(), ncols=4)
    report.plot(exec_plot, key="exec-plot")

    # About the report
    report.markdown("### About", key="about")
    report.markdown(
        "**Oxford Nanopore Technologies products are not intended for use for "
        "health assessment or to diagnose, treat, mitigate, cure or prevent "
        "any disease or condition.** This report was produced using the "
        "[epi2me-labs/wf-metagenomics](https://github.com/epi2me-labs/"
        "wf-metagenomics)."
        "The workflow can be run using `nextflow epi2me-labs/wf-metagenomics "
        "--help`",
        key="about-desc",
    )

    fname = "report.html"
    report.write(fname)


if __name__ == "__main__":
    main()
