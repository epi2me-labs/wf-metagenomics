"""Test the diversity module."""

import pandas as pd
from workflow_glue.diversity import global_rarefaction


def _abundance_table():
    """Create a samtools depth table with dummy data."""
    data = {
        'barcode01': [10, 10, 2, 9],
        'barcode02': [100, 0, 0, 0],
        'barcode03': [50, 30, 5, 15],
        'barcode04': [0, 0, 0, 1]
    }
    df = pd.DataFrame(data)
    return df


def test_001_global_rarefaction():
    """Test alignment_metrics function."""
    abundance_table = _abundance_table()
    expected = {
        'barcode02': [95, 0, 0, 0],
        'barcode03': [48, 28, 5, 14],
    }
    expected_df = pd.DataFrame(expected)
    actual = global_rarefaction(df=abundance_table, seed=4)
    pd.testing.assert_frame_equal(
        actual,
        expected_df,
        check_dtype=True,
        check_exact=False,
        rtol=0.001,
    )
