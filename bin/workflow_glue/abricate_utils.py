#!/usr/bin/env python
"""Tools for data wrangling abricate."""
import json

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


def main(args):
    """Run the entry point."""
    out_json = amr_df_to_json(args.sampleid, args.input)
    with open(args.output, "w") as outfile:
        json.dump(out_json, outfile)


def argparser():
    """Create argument parser."""
    parser = wf_parser("abricate_json")
    parser.add_argument("--sampleid", required=True)
    parser.add_argument("--input", required=True, nargs="+")
    parser.add_argument("--output", required=True)
    return parser
