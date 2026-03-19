import sys
from pathlib import Path
from unittest.mock import patch

from iris.cli import main


def test_cli_without_subcommand_prints_help(capsys, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["iris"])
    main()
    captured = capsys.readouterr()
    assert "usage:" in captured.out
    assert "features" in captured.out
    assert "predict" in captured.out


def test_cli_features_subcommand_calls_feature_generation(capsys, monkeypatch, tmp_path):
    input_dir = tmp_path / "complex"
    input_dir.mkdir()

    pymol_path = tmp_path / "pymol"
    pymol_path.write_text("")

    plugin_dir = tmp_path / "plugins"
    plugin_dir.mkdir()

    with patch("iris.cli.run_feature_generation") as mock_run_feature_generation:
        monkeypatch.setattr(
            sys,
            "argv",
            ["iris", "features", str(input_dir), str(pymol_path), str(plugin_dir)],
        )
        main()

    mock_run_feature_generation.assert_called_once_with(
        str(input_dir),
        str(pymol_path),
        str(plugin_dir),
    )

    captured = capsys.readouterr()
    assert "Features generated successfully!" in captured.out


def test_cli_predict_subcommand_calls_prediction(capsys, monkeypatch, tmp_path):
    input_dir = tmp_path / "complex"
    input_dir.mkdir()

    expected_output = input_dir / "complex_IRIS_sorted.sdf"

    with patch("iris.cli.run_prediction", return_value=expected_output) as mock_run_prediction:
        monkeypatch.setattr(
            sys,
            "argv",
            ["iris", "predict", str(input_dir)],
        )
        main()

    mock_run_prediction.assert_called_once_with(str(input_dir), model_path=None)

    captured = capsys.readouterr()
    assert "Processed and wrote sorted SDF file to" in captured.out
    assert str(expected_output) in captured.out


def test_cli_predict_subcommand_accepts_model_path(capsys, monkeypatch, tmp_path):
    input_dir = tmp_path / "complex"
    input_dir.mkdir()

    model_path = tmp_path / "custom_model.pkl"
    model_path.write_text("")

    expected_output = input_dir / "complex_IRIS_sorted.sdf"

    with patch("iris.cli.run_prediction", return_value=expected_output) as mock_run_prediction:
        monkeypatch.setattr(
            sys,
            "argv",
            ["iris", "predict", str(input_dir), "--model-path", str(model_path)],
        )
        main()

    mock_run_prediction.assert_called_once_with(
        str(input_dir),
        model_path=str(model_path),
    )

    captured = capsys.readouterr()
    assert "Processed and wrote sorted SDF file to" in captured.out