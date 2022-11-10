#!/usr/bin/env python
"""Script to create an aggregated count from lineage data."""
import argparse
import json
import logging
import sys

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


def update_or_create_count(entry, entries, bracken_counts):
    """Increment lineage counts given entries."""
    try:
        tax_id, lineage, ranks = entry.rstrip().split('\t')
    except ValueError:
        if str(entry).strip() == '0':
            return update_or_create_unclassified(entries)
        else:
            log = logging.getLogger(__name__)
            log.warning('Error: unknown lineage {}'.format(entry))
            sys.exit(1)

    lineage_split = lineage.split(';')
    ranks_split = ranks.split(';')
    count = int(bracken_counts[tax_id])

    previous = entries
    for [name, rank] in zip(lineage_split, ranks_split):

        if rank not in RANKS:
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
        perc = "{:.2%}".format(j['count'] / total)
        yield (indent, i, j['count'], perc, j['rank'])
        for k in yield_entries(j['children'], total, indent + 1):
            yield k


def main(prefix, lineages, bracken):
    """Run lineage aggregation algorithm."""
    bracken_counts = {}
    with open(bracken) as f:
        bracken = f.readlines()
    for i in bracken:
        bracken_counts[i.split()[0]] = i.split()[1]
    with open(lineages) as f:
        infile = f.readlines()
    entries = {}
    total = 0
    for line in infile:
        entries = update_or_create_count(line, entries, bracken_counts)
        total += 1
    output_report = open('{}.lineages.txt'.format(prefix), 'w')
    output_json = open('{}.lineages.json'.format(prefix), 'w')
    for entry in yield_entries(entries, total):
        [indent, name, count, perc, rank] = entry
        output_report.write(' '.join(
            ['-' * (indent + 1), name, str(count), perc, rank, '\n'])
        )
    output_json.write(json.dumps(entries))


def execute():
    """Parse command line arguments and run main."""
    parser = argparse.ArgumentParser(
        description="Aggregates lineage counts in a kraken2-like format",
    )
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
        '-p',
        help="Prefix to append to output file names.",
        dest="prefix",
        required=True,
    )

    args = parser.parse_args()

    main(
        lineages=args.lineages,
        prefix=args.prefix,
        bracken=args.bracken,
    )


if __name__ == "__main__":
    execute()
