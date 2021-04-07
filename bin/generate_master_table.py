#!/usr/bin/env python3
"""Script to split fastq by taxid using the output of centrifuge."""
import argparse
from pathlib import Path

import ete3
import pandas as pd


NCBI = ete3.NCBITaxa()
RANKS = [
    "superkingdom",
    "kingdom",
    "phylum",
    "class",
    "order",
    "family",
    "genus",
    "species",
]


parser = argparse.ArgumentParser(
    description=(
        "Script to split fastq by taxid using the output of centrifuge."
    )
)
parser.add_argument(
    "read_classification", type=Path, help="File produced by centrifuge"
)
parser.add_argument(
    "fastcat_output", type=Path, help="File produced by fastcat"
)
parser.add_argument("outdir", type=Path, help="Directory to place results")
split_type = parser.add_mutually_exclusive_group(required=True)
split_type.add_argument(
    "--rank",
    choices=RANKS,
    help=(
        "Group fastq by a specific rank.  Note: all "
        "reads classified above this level will be grouped in one group."
    ),
)
split_type.add_argument(
    "--taxid",
    type=int,
    help=(
        "Group reads classified at this rank and below in one file. "
        "Produces taxid group, unclassified and other."
    ),
)
split_type.add_argument(
    "--classified",
    action="store_true",
    default=False,
    help=(
        "If selected, reads that were successfully classified "
        "will be grouped separately from unclassified reads."
    ),
)


def _get_ranked_lineage(taxid):
    if not taxid:
        return {}
    try:
        lineage = NCBI.get_rank(NCBI.get_lineage(taxid))
        return {r: str(t) for t, r in lineage.items() if r in RANKS}
    except ValueError:
        # This will happen if the ID is not found or the dbs are out of date
        return {}


def _create_lineage_master(taxids):
    return pd.DataFrame(
        (_get_ranked_lineage(taxid) for taxid in taxids),
        index=taxids,
        columns=RANKS,
    )


def _get_lineage(taxid):
    try:
        return NCBI.get_lineage(int(taxid))
    except ValueError:
        return {}


def main():
    """Generate master table."""
    args = parser.parse_args()
    tsv_path = args.read_classification
    fastcat_output = args.fastcat_output
    output_dir = args.outdir
    split_by_rank = args.rank
    split_by_taxid = args.taxid
    split_by_status = args.classified

    output_dir.mkdir(exist_ok=True)

    if not tsv_path.is_file:
        raise FileNotFoundError
    elif split_by_taxid is not None and int(split_by_taxid) <= 0:
        raise Exception(f"Invalid taxid: {split_by_taxid}")
    elif not fastcat_output.is_file:
        raise FileNotFoundError

    read_stats = pd.read_csv(fastcat_output, sep="\t")
    read_stats.rename(columns={'read_id': 'readID'}, inplace=True)
    read_stats = read_stats.set_index('readID')
    classifications = pd.read_csv(tsv_path, sep="\t")

    tax_groups = classifications["taxID"].unique()
    assignments = None
    if split_by_rank:
        assignments = {
            taxid: "unclassified"
            if not taxid
            else _get_ranked_lineage(taxid).get(split_by_rank, "other")
            for taxid in tax_groups
        }
    elif split_by_taxid:
        assignments = {
            taxid: split_by_taxid
            if _get_lineage(int(taxid))
            and split_by_taxid in NCBI.get_lineage(int(taxid))
            else "unclassified"
            if taxid == 0
            else "other"
            for taxid in tax_groups
        }
    elif split_by_status:
        assignments = {
            taxid: "unclassified" if not taxid else "classified"
            for taxid in tax_groups
        }
    if not assignments:
        Exception("Something weird has happened")
    lineage_master_table = _create_lineage_master(tax_groups)
    assignment_df = pd.DataFrame.from_dict(
        assignments, orient="index", columns=["label"]
    )
    assignment_df = lineage_master_table.join(assignment_df)
    read_classification_master = classifications.join(
        assignment_df, on="taxID"
    )
    read_classification_master = read_classification_master.join(
        read_stats, on="readID"
    )
    print(read_classification_master.head(20))
    with open(output_dir / "read_classification_master.tsv", "wb") as f:
        read_classification_master.to_csv(f, na_rep="-1", index=False)

    print(
        read_classification_master[["readID", "label"]]
        .groupby("label")
        .count()
    )


if __name__ == "__main__":
    main()
