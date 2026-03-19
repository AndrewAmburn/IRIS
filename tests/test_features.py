from pathlib import Path
from unittest.mock import patch

import pytest

from iris.features import (
    get_legacy_scripts_dir,
    get_repo_root,
    run_feature_generation,
    run_legacy_script,
)


def test_get_repo_root_exists():
    repo_root = get_repo_root()
    assert repo_root.exists()
    assert repo_root.is_dir()


def test_get_legacy_scripts_dir_points_to_scripts_folder():
    scripts_dir = get_legacy_scripts_dir()
    assert scripts_dir.name == "scripts"


def test_run_legacy_script_raises_for_missing_script(tmp_path):
    with pytest.raises(FileNotFoundError):
        run_legacy_script(
            script_name="not_a_real_script.py",
            args=[],
            scripts_dir=tmp_path,
        )


def test_run_legacy_script_calls_subprocess(tmp_path):
    script_path = tmp_path / "dummy.py"
    script_path.write_text("print('hello')\n")

    with patch("iris.features.subprocess.run") as mock_run:
        run_legacy_script(
            script_name="dummy.py",
            args=["arg1", "arg2"],
            scripts_dir=tmp_path,
            python_executable="python",
            quiet=True,
            delay_seconds=0.0,
        )

    mock_run.assert_called_once()
    called_args = mock_run.call_args[0][0]

    assert called_args[0] == "python"
    assert called_args[1] == str(script_path)
    assert called_args[2:] == ["arg1", "arg2"]


def test_run_feature_generation_raises_for_invalid_working_directory(tmp_path):
    missing_dir = tmp_path / "missing_dir"
    pymol_path = tmp_path / "pymol"
    plugin_dir = tmp_path / "plugins"

    pymol_path.write_text("")
    plugin_dir.mkdir()

    with pytest.raises(NotADirectoryError):
        run_feature_generation(
            working_directory=missing_dir,
            pymol_path=pymol_path,
            pymol_plugin_path=plugin_dir,
            delay_seconds=0.0,
        )


def test_run_feature_generation_raises_for_missing_pymol(tmp_path):
    working_dir = tmp_path / "complex"
    working_dir.mkdir()
    plugin_dir = tmp_path / "plugins"
    plugin_dir.mkdir()

    with pytest.raises(FileNotFoundError):
        run_feature_generation(
            working_directory=working_dir,
            pymol_path=tmp_path / "missing_pymol",
            pymol_plugin_path=plugin_dir,
            delay_seconds=0.0,
        )


def test_run_feature_generation_raises_for_missing_plugin_dir(tmp_path):
    working_dir = tmp_path / "complex"
    working_dir.mkdir()
    pymol_path = tmp_path / "pymol"
    pymol_path.write_text("")

    with pytest.raises(FileNotFoundError):
        run_feature_generation(
            working_directory=working_dir,
            pymol_path=pymol_path,
            pymol_plugin_path=tmp_path / "missing_plugins",
            delay_seconds=0.0,
        )


def test_run_feature_generation_calls_scripts_in_order(tmp_path):
    working_dir = tmp_path / "complex"
    working_dir.mkdir()

    pymol_path = tmp_path / "pymol"
    pymol_path.write_text("")

    plugin_dir = tmp_path / "plugins"
    plugin_dir.mkdir()

    with patch("iris.features.run_legacy_script") as mock_run_legacy_script:
        run_feature_generation(
            working_directory=working_dir,
            pymol_path=pymol_path,
            pymol_plugin_path=plugin_dir,
            scripts_dir=tmp_path,
            python_executable="python",
            quiet=True,
            delay_seconds=0.0,
        )

    assert mock_run_legacy_script.call_count == 9

    called_script_names = [
        call.kwargs["script_name"] for call in mock_run_legacy_script.call_args_list
    ]

    assert called_script_names == [
        "sd2sdf.py",
        "consolidate_scores.py",
        "lig_sd2pdb.py",
        "consolidate_scores_and_rmsd.py",
        "generate_descriptors.py",
        "rdkit_full.py",
        "cav_volume.py",
        "distances.py",
        "receptor_descriptors.py",
    ]