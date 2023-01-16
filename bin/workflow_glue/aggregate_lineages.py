#!/usr/bin/env python
"""Script to create an aggregated count from lineage data."""
import json
import logging
import sys
from .util import wf_parser  # noqa: ABS101

UNCLASSIFIED = 'Unclassified'
UNKNOWN = 'Unknown'

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
    """Handle unclassified entries."""
    unclass = entries.get(UNCLASSIFIED)
    if not unclass:
        entries[UNCLASSIFIED] = {
            'rank': RANKS[0],
            'count': 1,
            'children': {
                UNKNOWN: {
                    'rank': "species",
                    'count': 1,
                    'children': {}
                }
            }
        }
    else:
        unclass['count'] += 1
        unclass['children'][UNKNOWN]['count'] += 1
    return entries


def update_or_create_count(entry, entries):
    """Increment lineage counts given entries."""
    try:
        _, lineage, ranks = entry.rstrip().split('\t')
    except ValueError:
        if str(entry).strip() == '0':
            return update_or_create_unclassified(entries)
        else:
            log = logging.getLogger(__name__)
            log.warning('Error: unknown lineage {}'.format(entry))
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
    """Get entries in printable form."""
    for i, j in entries.items():
        perc = "{:.2%}".format(j['count'] / total)
        yield (indent, i, j['count'], perc, j['rank'])
        for k in yield_entries(j['children'], total, indent + 1):
            yield k


def main(args):
    """Run lineage aggregation algorithm."""
    lineages = args.lineages
    prefix = args.prefix
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


def argparser():
    """Create argument parser."""
    parser = wf_parser("aggregate_lineages")
    parser.add_argument(
        '-i',
        help=(
            "Lineages .tsv (taxid, lineage)."
        ),
        dest="lineages",
        metavar='',
        default=None
    )  # not sure if default can be none here

    parser.add_argument(
        '-p',
        help="Prefix to append to output file names.",
        dest="prefix",
        required=True,
        metavar='',
    )
    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
