"""Test report errors when empty files are input."""

import pytest
from workflow_glue import report


BASE_PARAMS = [
    "wf-metagenomics-report.html",
    "--workflow_name", "wf-metagenomics-test",
    "--pipeline", "minimap2",
    "--taxonomic_rank", "S",
]

INPUT_PARAMS = [
    # `[filename, flag, is_dir, no_test]`
    ["versions", "--versions", True, True],
    ["params", "--params", False, True],
    ["metadata", "--metadata", False, False],
    ["stats", "--read_stats", True, False],
    ["lineages", "--lineages", True, True],
    ["abundance_table_species", "--abundance_table", False, True],
    ["alignment_stats", "--align_stats", True, True],
    ["amr", "--amr", True, True],
]


def test_empty_inputs(tmp_path):
    """Test that the report script properly checks that relevant inputs aren't empty.

    For each file/directory input param that is not allowed to be empty, create
    an argument list with current param pointing to an empty file/dir and all
    other params pointing to valid file/directory.
    tmp_path is a pytest fixture to create tmp dir for testing

    """
    for empty_param, flag, is_empty_dir, no_test in INPUT_PARAMS:
        # Reset args list in each loop, to determine which param
        # will have the empty file/dir
        args = BASE_PARAMS.copy()
        if not no_test:
            # don't need to test this one as it can be empty (e.g. `metadata` in
            # real-time)
            continue
        # use an empty file / dir for `input_name`
        empty_input_dir = tmp_path / f'{empty_param}_empty'
        empty_input_dir.mkdir()
        if is_empty_dir:
            args += [flag, str(empty_input_dir)]
        else:
            # Add empty file
            empty_input_file = empty_input_dir / 'empty.txt'
            empty_input_file.touch()
            args += [flag, str(empty_input_file)]
        # Build args for the remaining non-empty params
        for new_input_name, flag, is_dir, _ in INPUT_PARAMS:
            if new_input_name != empty_param:
                # Again can be a dir with files or just a file
                fname = tmp_path / new_input_name
                if is_dir:
                    # The dir should contain non empty files
                    fname.mkdir(parents=True, exist_ok=True)
                    (fname / f"{new_input_name}.txt").write_text('blo')
                else:
                    fname.write_text('blo')
                # Add the flag and folder/file to args.
                args += [flag, str(fname)]
        # Define pattern for the expected message
        pattern = f"Empty {'directory' if is_empty_dir else 'input file'}"
        with pytest.raises(SystemExit, match=pattern):
            report.argparser().parse_args(args)
