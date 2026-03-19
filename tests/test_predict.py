from pathlib import Path

import pytest

from iris.predict import load_model, read_scores_file, run_prediction


def test_read_scores_file_parses_multiple_ligands(tmp_path):
    descriptor_file = tmp_path / "descriptors_all.txt"
    descriptor_file.write_text(
        "Ligand 1 Scores:\n"
        "SCORE: -32.5\n"
        "FeatureA: 1.2\n"
        "FeatureB: 3.4\n"
        "\n"
        "Ligand 2 Scores:\n"
        "SCORE: -28.1\n"
        "FeatureA: 5.6\n"
        "FeatureB: 7.8\n"
    )

    df = read_scores_file(descriptor_file)

    assert df.shape[0] == 2
    assert "SCORE" in df.columns
    assert "FeatureA" in df.columns
    assert "FeatureB" in df.columns
    assert df.iloc[0]["SCORE"] == -32.5
    assert df.iloc[1]["FeatureA"] == 5.6


def test_read_scores_file_raises_for_missing_file():
    with pytest.raises(FileNotFoundError):
        read_scores_file("/this/file/does/not/exist.txt")


def test_read_scores_file_raises_for_empty_parse(tmp_path):
    descriptor_file = tmp_path / "descriptors_all.txt"
    descriptor_file.write_text("This is not a valid descriptor file\n")

    with pytest.raises(ValueError):
        read_scores_file(descriptor_file)


def test_load_model_raises_for_missing_file():
    with pytest.raises(FileNotFoundError):
        load_model("/this/model/does/not/exist.pkl")


def test_run_prediction_raises_for_missing_directory():
    with pytest.raises(FileNotFoundError):
        run_prediction("/this/path/does/not/exist")


def test_run_prediction_raises_for_file_instead_of_directory(tmp_path):
    not_a_directory = tmp_path / "example.txt"
    not_a_directory.write_text("test")

    with pytest.raises(NotADirectoryError):
        run_prediction(not_a_directory)