# IRIS: Intelligent RNA Interaction Scorer

IRIS is a machine learning-based tool for improving RNA-ligand docking pose ranking. It generates structure-based features for a docked RNA-ligand complex and re-ranks docking poses using a trained machine learning model.

## Installation

Create and activate the Conda environment:

```bash
conda env create -f environment.yml
conda activate iris_env
```

To use the `iris` command-line interface, install IRIS from the repository root:

```bash
pip install -e ".[dev]"
```

## Additional Requirements

IRIS feature generation requires:

- A working local installation of PyMOL
- The PyMOL plugins:
  - `show_bumps.py`
  - `get_raw_distances.py`

Example PyMOL executable path on macOS:

```bash
/Applications/PyMOL.app/Contents/bin/pymol
```

Example PyMOL plugin startup path on macOS:

```bash
/Applications/PyMOL.app/Contents/lib/python3.9/site-packages/pmg_tk/startup
```

## Required Input Folder Structure

IRIS expects one folder per RNA-ligand complex. The folder name is typically a four-character identifier such as the receptor PDB ID.

Each complex folder should contain:

- `<folder_name>.pdb`  
  RNA receptor structure. If no PDB is available, copy the `.mol2` file and rename it as `<folder_name>.pdb`.

- `<folder_name>.mol2`  
  Receptor structure without ligand bound.

- `<folder_name>_lig.sd`  
  Ligand structure file.

- `<folder_name>_docking_out_sorted.sdf`  
  rDock docking output.

- `<folder_name>_cavity.log`  
  rDock cavity log file.

## Usage

IRIS currently runs in two stages.

### 1. Generate features

```bash
iris features <folder_path> <pymol_path> <pymol_plugin_path>
```

Example:

```bash
iris features \
IRIS_example/1uts \
/Applications/PyMOL.app/Contents/bin/pymol \
/Applications/PyMOL.app/Contents/lib/python3.9/site-packages/pmg_tk/startup
```

### 2. Re-rank docking poses

```bash
iris predict <folder_path>
```

Example:

```bash
iris predict IRIS_example/1uts
```

## Example

IRIS includes an example docking case in:

```text
IRIS_example/1uts
```

To test the installation:

```bash
cd IRIS
conda activate iris_env

iris features \
IRIS_example/1uts \
/<pymol_path> \
/<pymol_plugin_path>

iris predict IRIS_example/1uts
```

## Output

Successful execution produces a re-ranked SDF file:

```text
<folder_name>_IRIS_sorted.sdf
```

For the included example, the main output is:

```text
1uts_IRIS_sorted.sdf
```

This file contains the docking poses reordered by IRIS.

## Troubleshooting

Most errors come from one of the following:

- Missing required input files
- Incorrect PyMOL executable path
- Incorrect PyMOL plugin path
- Running prediction before feature generation
- Incorrect folder naming or file naming

If prediction fails, confirm that feature generation created:

```text
descriptors_all.txt
```

inside the target complex folder.

## Notes

- Run commands from the repository root unless using absolute paths.
- Feature generation and prediction are separate steps.
- PyMOL and its required plugins must be installed locally before running `iris features`.

## Citation

If you use IRIS in your work, please cite the associated manuscript at the following: https://doi.org/10.1021/acsomega.5c11891
