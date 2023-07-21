#!/usr/bin/env python
"""Tools for data wrangling abricate."""
import json
import os

import pandas as pd

from .util import wf_parser  # noqa: ABS101


def amr_df_to_json(sample, amr_input):
    """Create json of relevant information from abricate df.

    :param sample (str): Name of sample being analysed.
    :param amr_input ([str]): list of string indicating path of abricate tsv results.
    :return (dict): Sample amr results in json structure.
    """
    for i in amr_input:
        amr_df = pd.read_csv(i, sep="\t")
        if amr_df.empty:
            return {sample: {"pass": False}}
        results_dict = dict()
        for product, df in amr_df.groupby("PRODUCT"):
            df = df.fillna("Unknown")
            meta = df[
                ["SEQUENCE", "START", "END", "%COVERAGE", "%IDENTITY", "RESISTANCE"]
                ].to_dict(orient="records")
            results_dict[product] = {"count": df.shape[0], "meta": meta}
    return {sample: {"pass": True, "results": results_dict}}


def combine_jsons(sample, inputs):
    """Update json of sample during watch path.

    :param sample (str): Name of sample being analysed.
    :param inputs ([str]): List of file paths with current results and previous
        results that are to be combined.
    :return (dict): Combined amr results in json structure.
    """
    combined = dict()
    results_pass = False
    for i in inputs:
        if os.stat(i).st_size == 0:
            continue
        with open(i) as f:
            json_raw = json.load(f)
        sample_data = json_raw[sample]
        if sample_data["pass"]:
            results = sample_data["results"]
            for gene, data in results.items():
                if gene not in combined:
                    combined[gene] = data
                else:
                    combined[gene]["count"] += data["count"]
                    combined[gene]["meta"].extend(data["meta"])
            results_pass = True
    return {sample: {"pass": results_pass, "results": combined}}


def main(args):
    """Run the entry point."""
    if not args.combine:
        out_json = amr_df_to_json(args.sampleid, args.input)
    else:
        out_json = combine_jsons(args.sampleid, args.input)
    with open(args.output, "w") as outfile:
        json.dump(out_json, outfile)


def argparser():
    """Create argument parser."""
    parser = wf_parser("abricate_json")
    parser.add_argument("--sampleid", required=True)
    parser.add_argument("--input", required=True, nargs="+")
    parser.add_argument("--output", required=True)
    parser.add_argument(
        "--combine", action="store_true",
        help="Combine jsons for watch_path"
                        )
    return parser
