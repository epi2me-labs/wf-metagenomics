"""Test the table generation code.

Test the functions in bin/workflow_glue/report_utils.py that are used to pick
the most abundant taxa and to plot and filter them by abundance thresholds.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from workflow_glue.report_utils.report_utils import (
    filter_by_abundance,
    most_abundant_table,
    RANK_ORDER,
    split_taxonomy_string
)


@pytest.fixture
def test_data(request):
    """Define data location fixture."""
    return Path(request.config.getoption("--test_data")) / "workflow_glue/case03"


def _compare_frames(actual, expected):
    """Compare two dataframes."""
    pd.testing.assert_frame_equal(
        actual,
        expected,
        check_dtype=True,
        check_categorical=True,
        check_exact=True,
        check_datetimelike_compat=True,
    )


def _abundance_table():
    """Create an abundance table with dummy data."""
    species_names = [
        "previous_rank;Sp1",
        "previous_rank;Sp2",
        "previous_rank;Sp3",
        "previous_rank;Sp4",
        "previous_rank;Sp5",
        "previous_rank;Sp6",
        "previous_rank;Sp7",
        "previous_rank;Unclassified",
    ]
    data = {
        "sample_1": [1, 80, 0, 20, 0, 0, 1, 0],
        "sample_2": [1, 50, 25, 10, 5, 1, 0, 9],
        "sample_3": [1, 20, 20, 20, 20, 20, 0, 0],  # equally distributed
        "sample_4": [1, 0, 70, 30, 0, 0, 0, 200],  # more than 100 counts
        "sample_5": [0, 0, 0, 0, 0, 0, 0, 0],  # 0 counts
        "sample_Uncl": [0, 0, 0, 0, 0, 0, 0, 100],
    }
    df = pd.DataFrame(data)
    df["total"] = df.sum(axis=1)
    df["tax"] = species_names
    return df


@pytest.mark.parametrize(
    "threshold, species_to_be_removed", [
        (0,  []),
        (1,  ["previous_rank;Sp7"]),
        (0.01, ["previous_rank;Sp1", "previous_rank;Sp7"]),
        (25, [
            "previous_rank;Sp1", "previous_rank;Sp5",
            "previous_rank;Sp6", "previous_rank;Sp7"]),
    ]
)
def test_001_filter_by_abundance(threshold, species_to_be_removed):
    """Test the filtering based on thresholds."""
    df = _abundance_table()
    actual = filter_by_abundance(
        df, column_to_filter="total", abundance_threshold=threshold
    )
    expected = _abundance_table()
    expected = expected[~expected.tax.isin(species_to_be_removed)]
    _compare_frames(
        actual, expected)


def test_002_most_abundant_absolute():
    """Check most abundant function."""
    df = _abundance_table()
    actual = most_abundant_table(df, n=5, percent=False)
    species_names_expected = [
        "previous_rank;Unclassified",
        "previous_rank;Sp2",
        "previous_rank;Sp3",
        "previous_rank;Sp4",
        "previous_rank;Sp5",
        "Other;Other",
    ]
    expected_d = {  # last value is others (sum of not selected species)
        "sample_1": [0, 80, 0, 20, 0, 2],
        "sample_2": [9, 50, 25, 10, 5, 2],
        "sample_3": [0, 20, 20, 20, 20, 21],  # equally distributed
        "sample_4": [200, 0, 70, 30, 0, 1],  # more than 100 counts
        "sample_5": [0, 0, 0, 0, 0, 0],
        "sample_Uncl": [100, 0, 0, 0, 0, 0],
        "total": [309, 150, 115, 80, 25, 26]
    }
    expected = pd.DataFrame(expected_d)
    expected["tax"] = species_names_expected
    _compare_frames(actual.reset_index(drop=True), expected.reset_index(drop=True))


def test_003_most_abundant_percent():
    """Check most abundant using percentages."""
    df = _abundance_table()
    actual = most_abundant_table(df, n=3, percent=True)
    species_names_expected = [
        "previous_rank;Unclassified",
        "previous_rank;Sp2",
        "previous_rank;Sp3",
        "Other;Other",
    ]
    expected_d = {
        "tax": species_names_expected,
        "sample_1": [0, 78.43, 0.0, 21.57],
        "sample_2": [8.91, 49.50, 24.75, 16.83],
        "sample_3": [0.0, 19.80, 19.80, 60.40],
        "sample_4": [66.45, 0.0, 23.26, 10.30],
        # 0 divided by 0 total reads
        "sample_5": [np.nan, np.nan, np.nan, np.nan],
        "sample_Uncl": [100, 0.0, 0.0, 0.0],
        "total": [43.83, 21.28, 16.31, 18.58]
    }
    expected = pd.DataFrame(expected_d)
    _compare_frames(actual, expected)


@pytest.mark.parametrize(
    "expected_columns",
    [
        (['superkingdom', 'kingdom']),
        (['superkingdom', 'kingdom', 'phylum']),
    ],
)
def test_004_split_taxonomy_string(expected_columns):
    """Test taxonomy string splitting."""
    # Taxonomy strings to split
    df_tax_str = pd.DataFrame({
        "tax": [
            "Bacteria;None;Proteobacteria", "Bacteria;None;Cyanobacteria",
            "Archaea;None;Crenarchaeota", "Eukarya;Viridiplantae",
            "Archaea;None;Crenarchaeota", "Archaea;None;Crenarchaeota",
            "Eukarya;Fungi;Basidiomycota", "Unclassified"
        ],
    })

    expected = pd.DataFrame({
        'superkingdom': [
            'Bacteria', 'Bacteria', 'Archaea', 'Eukarya',
            'Archaea', 'Archaea', 'Eukarya', 'Unclassified'],
        'kingdom': [
            'None', 'None', 'None', 'Viridiplantae',
            'None', 'None', 'Fungi', 'None'],
        'phylum': [
            'Proteobacteria', 'Cyanobacteria', 'Crenarchaeota', 'None',
            'Crenarchaeota', 'Crenarchaeota', 'Basidiomycota', 'None']
    })
    expected = expected[expected_columns]
    actual = split_taxonomy_string(df_tax_str)[expected_columns]
    _compare_frames(expected.astype('str'), actual.astype('str'))


# Compare abundances in the abundance table
# with raw abundances from the JSON file
@pytest.mark.parametrize(
    "threshold",
    [
        40,  # remove Vibrio sinaloensis, Limnococcus fonticola
    ],
)
@pytest.mark.parametrize("rank", ["species"])
def test_101_filter_by_abundance_file(
    test_data, threshold, rank, sample_name="test_sample"
):
    """Run abundance related functions using realistic files."""
    input_table = test_data / "abundance_table.tsv"
    expected_table = test_data / "test_expected/test_101_species_40.tsv"
    abundance_table_rank = pd.read_csv(input_table, sep='\t')
    # Prepare actual results as it is done in report.
    # update full taxonomy string to a string which finishes at the rank asked.
    abundance_table_rank['tax'] = [';'.join(
            i.split(';')[:RANK_ORDER[rank]]) for i in abundance_table_rank['tax']]
    # Group samples by rank
    counts_per_taxa_df = abundance_table_rank.groupby(['tax']).sum().reset_index()
    # Filter by abundance threshold.
    # Distinguish between natural number cutoff or percentage.
    actual_counts = filter_by_abundance(counts_per_taxa_df, 'total', threshold)
    # Select just counts and split taxonomy string to last level
    actual_counts['species'] = actual_counts['tax'].apply(lambda x: x.split(';')[-1])
    actual = actual_counts[['species', sample_name]]
    actual = actual.sort_values(by=['species']).set_index('species')
    # Expected results
    expected_counts = pd.read_csv(expected_table, sep='\t')
    expected_counts = expected_counts.sort_values(by=['species']).set_index('species')
    # Compare if selected taxa are the same
    _compare_frames(
        expected_counts,
        actual
    )
