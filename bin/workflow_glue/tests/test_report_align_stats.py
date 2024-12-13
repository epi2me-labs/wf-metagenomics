"""Test the table generation code.

Test the alignment statistic helper functions in bin/workflow_glue/report_utils.py.
"""

from pathlib import Path

import pandas as pd
import pytest
from workflow_glue.report_utils.report_utils import (
    alignment_metrics,
    depth2heatmap,
)


@pytest.fixture
def valid_inputs_dir(request):
    """Define data location fixture."""
    return Path(
        request.config.getoption(
            "--test_data")) / "workflow_glue/case03/alignment_stats"


def _depth_table():
    """Create a samtools depth table with dummy data."""
    data = {
        'ref': ['A', 'A', 'A', 'A', 'A', 'B', 'B', 'B', 'C'],
        'pos': [2, 3, 5, 7, 9, 1, 2, 3, 1],
        'depth': [2, 2, 1, 2, 3, 4, 4, 4, 1]
    }
    df = pd.DataFrame(data)
    return df


def _cov_table():
    """Create a simplified samtools coverage table with dummy data."""
    data = {
        'ref': ['B', 'A', 'C'],
        'endpos': [3, 10, 4],
    }
    df = pd.DataFrame(data)
    return df.set_index('ref')


@pytest.mark.parametrize(
    "expected", [
        {
            'ref': ['B', 'A', 'C'],
            'endpos': [3, 10, 4],
            'mean': [4.0, 1.0, 0.25],
            'sd': [0.0, 1.095, 0.433],
            'Coefficient of Variance': [0.0, 1.095, 1.732],
        }
    ]
)
def test_001_alignment_metrics(expected):
    """Test alignment_metrics function."""
    depth_df = _depth_table()
    cov_df = _cov_table()

    expected_df = pd.DataFrame(expected)
    actual = alignment_metrics(depth=depth_df, stats=cov_df)
    pd.testing.assert_frame_equal(
        actual,
        expected_df,
        check_dtype=True,
        check_categorical=True,
        check_exact=False,
        rtol=0.001,
        check_datetimelike_compat=True,
    )


@pytest.mark.parametrize(
    "windows, expected", [
        # windows = 1, the results are the same as the mean_depth from alignment metrics
        # no percentiles, just 1
        (1, {
            'ref': ['B', 'A'],
            0: [4.0, 1.0],
        }),
        (2, {
            'ref': ['B', 'A'],
            0: [8.0, 1.0],
            1: [4.0, 1.0],
        }),
    ]
)
def test_002_depth2heatmap(windows, expected):
    """Test depth2heatmap function."""
    depth_df = _depth_table()
    cov_df = _cov_table()
    expected_df = pd.DataFrame(expected).set_index('ref')
    expected_df.columns = [int(val) for val in expected_df.columns]
    actual = depth2heatmap(depth=depth_df, reference=cov_df, windows=windows)
    pd.testing.assert_frame_equal(
        actual,
        expected_df,
        check_dtype=True,
        check_categorical=True,
        check_exact=False,
        rtol=0.001,
        check_datetimelike_compat=True,
    )
