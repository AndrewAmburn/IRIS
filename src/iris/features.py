from __future__ import annotations

import subprocess
import time
from pathlib import Path


def get_repo_root() -> Path:
    return Path(__file__).resolve().parents[2]


def get_legacy_scripts_dir() -> Path:
    return get_repo_root() / "scripts"


def run_legacy_script(
    script_name: str,
    args: list[str],
    scripts_dir: str | Path | None = None,
    python_executable: str = "python",
    quiet: bool = True,
    delay_seconds: float = 3.0,
) -> None:
    scripts_path = Path(scripts_dir) if scripts_dir is not None else get_legacy_scripts_dir()
    script_path = scripts_path / script_name

    if not script_path.is_file():
        raise FileNotFoundError(f"Legacy script not found: {script_path}")

    command = [python_executable, str(script_path), *args]

    stdout = subprocess.DEVNULL if quiet else None
    stderr = subprocess.DEVNULL if quiet else None

    try:
        subprocess.run(
            command,
            check=True,
            stdout=stdout,
            stderr=stderr,
            cwd=str(scripts_path),
        )
        if delay_seconds > 0:
            time.sleep(delay_seconds)
    except subprocess.CalledProcessError as exc:
        raise RuntimeError(f"Error running {script_name}: {exc}") from exc


def run_feature_generation(
    working_directory: str | Path,
    pymol_path: str | Path,
    pymol_plugin_path: str | Path,
    scripts_dir: str | Path | None = None,
    python_executable: str = "python",
    quiet: bool = True,
    delay_seconds: float = 3.0,
) -> None:
    working_directory = Path(working_directory).resolve()
    pymol_path = Path(pymol_path).resolve()
    pymol_plugin_path = Path(pymol_plugin_path).resolve()

    if not working_directory.is_dir():
        raise NotADirectoryError(f"Working directory is not valid: {working_directory}")

    if not pymol_path.exists():
        raise FileNotFoundError(f"PyMOL executable not found: {pymol_path}")

    if not pymol_plugin_path.exists():
        raise FileNotFoundError(f"PyMOL plugin directory not found: {pymol_plugin_path}")

    scripts = [
        ("sd2sdf.py", [str(working_directory)]),
        ("consolidate_scores.py", [str(working_directory)]),
        ("lig_sd2pdb.py", [str(working_directory)]),
        ("consolidate_scores_and_rmsd.py", [str(working_directory)]),
        (
            "generate_descriptors.py",
            [str(working_directory), str(pymol_path), str(pymol_plugin_path)],
        ),
        ("rdkit_full.py", [str(working_directory)]),
        ("cav_volume.py", [str(working_directory)]),
        ("distances.py", [str(working_directory)]),
        ("receptor_descriptors.py", [str(working_directory)]),
    ]

    for script_name, args in scripts:
        run_legacy_script(
            script_name=script_name,
            args=args,
            scripts_dir=scripts_dir,
            python_executable=python_executable,
            quiet=quiet,
            delay_seconds=delay_seconds,
        )
