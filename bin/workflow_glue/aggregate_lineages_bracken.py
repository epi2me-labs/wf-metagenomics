#!/usr/bin/env python
"""Script to create an aggregated count from lineage data."""
import json
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

RANKS_ABB = {
    "D": "superkingdom", "K": "kingdom", "P": "phylum", "C": "class",
    "O": "order", "F": "family", "G": "genus", "S": "species"
}


def update_or_create_unclassified(entries, unclassified_count):
    """Handle unclassified entries."""
    entries[UNCLASSIFIED] = {
        'rank': RANKS[0],
        'count': int(unclassified_count),
        'children': {
            UNKNOWN: {
                'rank': SELECTED_RANKS[-1],
                'count': int(unclassified_count),
                'children': {}
            }
        }
    }
    return entries


def update_or_create_count(entry, entries, bracken_counts):
    """Increment lineage counts given entries."""
    tax_id, lineage, ranks = entry.rstrip().split('\t')
    lineage_split = lineage.split(';')
    ranks_split = ranks.split(';')
    count = int(bracken_counts[tax_id])

    previous = entries
    for [name, rank] in zip(lineage_split, ranks_split):

        # Add this temporarily [https://github.com/DerrickWood/kraken2/issues/739]
        if rank not in SELECTED_RANKS or name in ['Holozoa', 'Nucletmycea']:
            continue

        current = previous.get(name)
        if not current:
            new_entry = {
                'rank': rank,
                'count': count,
                'children': {}
            }
            previous[name] = new_entry
            previous = new_entry['children']
            continue

        current['count'] += count
        previous = current['children']

    return entries


def yield_entries(entries, total, indent=0):
    """Get entries in printable form."""
    for i, j in entries.items():
        if total != 0:
            perc = "{:.2%}".format(j['count'] / total)
        else:
            perc = "0.00"
        yield (indent, i, j['count'], perc, j['rank'])
        for k in yield_entries(j['children'], total, indent + 1):
            yield k


def main(args):
    """Run lineage aggregation algorithm."""
    lineages = args.lineages
    prefix = args.prefix
    bracken = args.bracken
    report = args.report
    taxonomic_rank = args.taxonomic_rank
    bracken_counts = {}
    entries = {}
    total = 0
    global SELECTED_RANKS
    SELECTED_RANKS = []
    for i in RANKS:
        # append each rank
        SELECTED_RANKS.append(i)
        if i == RANKS_ABB[taxonomic_rank]:
            # stop when the user rank has been added
            break
    with open(bracken) as f:
        bracken = f.readlines()
    if len(bracken) > 0:
        for i in bracken:
            bracken_counts[i.split()[0]] = i.split()[1]
        with open(lineages) as f:
            infile = f.readlines()
        for line in infile:
            try:
                entries = update_or_create_count(line, entries, bracken_counts)
                total += 1
            except ValueError:
                sys.stderr.write(
                    """Lineage for tax id {} not found in taxonomy database"""
                    .format(str(line)))
    with open(report) as f:
        report_file = f.readlines()
        for line in report_file:
            # the last field could be something like "Unclassified Marinobacter"; we
            # thus check if the last field is exactly "unclassified"
            if "unclassified" == line.split()[-1]:
                unclassified_count = line.split()[1]
                entries = update_or_create_unclassified(
                    entries, unclassified_count)
                total += int(unclassified_count)
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
    parser = wf_parser("aggregate_lineages_bracken")
    parser.add_argument(
        '-i',
        help=(
            "Lineages .tsv (taxid, lineage)."
        ),
        dest="lineages",
    )

    parser.add_argument(
        '-b',
        help=(
            "Bracken Lineages .tsv (taxid, count)."
        ),
        dest="bracken",
    )

    parser.add_argument(
        '-u',
        help=(
            "full report to get unclassified count"
        ),
        dest="report",
    )

    parser.add_argument(
        '-p',
        help="Prefix to append to output file names.",
        dest="prefix",
        required=True,
    )

    parser.add_argument(
        '-r',
        help=("Taxonomic rank."),
        dest="taxonomic_rank",
        required=True,
    )

    return parser


if __name__ == "__main__":
    args = argparser().parse_args()
    main(args)
