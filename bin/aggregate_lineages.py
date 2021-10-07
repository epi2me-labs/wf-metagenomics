#!/usr/bin/env python
"""Script to check that sample sheet is well-formatted."""
import sys
import json
import argparse
from operator import itemgetter
from typing import TypedDict, Dict, List, Union, Tuple

UNCLASSIFIED = 'Unclassified'

RANKS = [
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species"
]


def update_or_create_unclassified(entries):
    unclass = entries.get(UNCLASSIFIED)
    if not unclass:
        entries[UNCLASSIFIED] = {
            'rank': UNCLASSIFIED,
            'count': 1,
            'children': {}
        }
    else:
        unclass['count'] += 1
    return entries


def update_or_create_count(entry, entries):
    try:
        _, lineage, ranks = entry.rstrip().split('\t')
    except ValueError:
        if str(entry).strip() == '0':
            return update_or_create_unclassified(entries)
        else:
            print('Error: unknown lineage {}'.format(entry))
            sys.exit(1)

    lineage_split = lineage.split(';')
    ranks_split = ranks.split(';')

    previous = entries
    for [name, rank] in zip(lineage_split, ranks_split):

        if rank not in RANKS:
            continue

        current = previous.get(name)
        if not current:
            new_entry = {
                'rank': rank,
                'count': 1,
                'children': {}
            }
            previous[name] = new_entry
            previous = new_entry['children']
            continue

        current['count'] += 1
        previous = current['children']

    return entries


def yield_entries(entries, total, indent=0):
    for i, j in entries.items():
        perc = "{:.2%}".format(j['count'] / total)
        yield(indent, i, j['count'], perc, j['rank'])
        for k in yield_entries(j['children'], total, indent + 1):
            yield k


def main(
    prefix: str,
    lineages: Union[str, None] = None
) -> None:
    """
    For each alignment received, writes a formatted
    row in TSV format to the designated path which 
    contains the assigned taxid and classification
    status.
    """
    if lineages:
        infile = open(lineages)
    else:
        infile = sys.stdin

    entries = {}
    total = 0
    for line in infile:
        entries = update_or_create_count(line, entries)
        total += 1

    output_report = open('{}.lineages.txt'.format(prefix), 'w')
    output_json = open('{}.lineages.json'.format(prefix), 'w')
    for entry in yield_entries(entries, total):
        [indent, name, count, perc, rank] = entry
        output_report.write(' '.join(
            ['-' * (indent + 1), name, str(count), perc, rank, '\n'])
        )
    output_json.write(json.dumps(entries))


def execute(argv) -> None:
    """
    Parses command line arguments and runs main.
    """
    parser = argparse.ArgumentParser(
        description="Outputs aggregated lineage counts in a kraken2-like format",
    )

    parser.add_argument(
        '-i',
        help=(
            "Lineages .tsv (taxid, lineage)."
        ),
        dest="lineages",
        metavar='',
    )

    parser.add_argument(
        '-p',
        help="Prefix to append to output file names.",
        dest="prefix",
        required=True,
        metavar='',
    )

    args = parser.parse_args(argv)

    main(
        lineages=args.lineages,
        prefix=args.prefix,
    )


if __name__ == "__main__":
    execute(sys.argv[1:])
