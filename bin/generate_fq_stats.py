#!/usr/bin/env python3
"""
Add fastq stats to master table.

positional arguments:
  fastq       Path to fastq/fastq.tar.gz

optional arguments:
  -h, --help  show this help message and exit
  -f F        Path to master centrifuge table with label
"""
import argparse
from functools import lru_cache
from math import log10
from pathlib import Path

import pandas as pd
from pysam import FastxFile


parser = argparse.ArgumentParser(
    description="Add fastq stats to master table."
)
parser.add_argument("fastq", type=Path, help="Path to fastq/fastq.tar.gz")
parser.add_argument(
    "-f",
    type=Path,
    help="Path to master centrifuge table with label",
    default=None,
)


def qstring_to_phred(quality):
    """Compute standard phred scores from a quality string."""
    qscores = [ord(q) - 33 for q in quality]
    return qscores


@lru_cache(100)
def _get_probability(qscore):
    return pow(10, -0.1 * qscore)


def mean_qscore(scores):
    """
    Return the phred score corresponding to the mean of the probabilities.

    :param scores: Iterable of phred scores.

    :returns: Phred score corresponding to the average error rate, as
        estimated from the input phred scores.
    """
    if len(scores) == 0:
        return 0.0
    sum_prob = 0.0
    for val in scores:
        sum_prob += _get_probability(val)
    mean_prob = sum_prob / len(scores)
    return -10.0 * log10(mean_prob)


def _process_record(fq_path):
    with FastxFile(fq_path) as fh:
        for fq_record in fh:
            mean = round(mean_qscore(qstring_to_phred(fq_record.quality)), 2)
            yield [fq_record.name, len(fq_record.sequence), mean]


def main(argv=None):
    """Process master table."""
    args = parser.parse_args(argv)
    fq_path = args.fastq
    master_path = args.f

    master = pd.DataFrame(
        _process_record(fq_path), columns=["readID", "len", "meanqscore"]
    )
    if master_path:
        master = master.set_index("readID")
        classifications = pd.read_csv(master_path, sep=",", index_col=0)
        read_master = classifications.join(master)
    else:
        read_master = master

    # print(read_master)
    with open("read_classification_master.tsv", "wb") as f:
        read_master.to_csv(f, na_rep="-1", index=True)


if __name__ == "__main__":
    main()
