#!/usr/bin/env python
"""Histogram-json."""

import json
import os

from .util import wf_parser  # noqa: ABS101


def add_dicts(d1, d2):
    """Extend json, sum values."""
    def sum_a(v1, v2):
        if v2 is None:
            return v1
        try:
            if isinstance(v1 + v2, int):
                return v1 + v2
            elif isinstance(v1 + v2, str):
                return v1
        except TypeError:
            return add_dicts(v1, v2)
    result = d2.copy()
    result.update({k: sum_a(v, d2.get(k)) for k, v in d1.items()})
    return result


def main(args):
    """Run the entry point."""
    if os.stat(args.state).st_size == 0:
        state = {}
    else:
        with open(args.state) as json_file:
            state = json.load(json_file)
    with open(args.new_file) as json_file:
        new_file = json.load(json_file)
    combined = add_dicts(state, new_file)
    with open(args.output, "w") as outfile:
        json.dump(combined, outfile)


def argparser():
    """Create argument parser."""
    parser = wf_parser("add_jsons")
    parser.add_argument("new_file")
    parser.add_argument("state")
    parser.add_argument("output")
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
