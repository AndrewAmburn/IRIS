# IRIS: Intelligent RNA Interaction Scorer  
A machine learning-based tool for improving RNA-ligand docking pose ranking.

## Installation Requirements  

### 1. Setting Up the Environment  
IRIS requires a Conda environment with Python 3.9 or higher. To install the necessary dependencies, follow these steps:

#### Install Conda Dependencies  
```bash
conda install --yes --file ~/IRIS/requirements.txt -c conda-forge
```

#### Install Additional Packages via Pip  
```bash
pip install spyrmsd==0.8.0 rdkit==2024.3.2 scikit-learn==1.3.2 joblib==1.3.2 pandas catboost==1.2.7
```

### 2. Additional Requirements  
In addition to the dependencies listed in `requirements.txt`, the following software and plugins are required:

- **PyMOL** (A working installation is needed for structure visualization and feature extraction).  
- **PyMOL Plugins**:
  - `show_bumps.py`
  - `get_raw_distances.py`  

Ensure that these plugins are accessible within your PyMOL installation before proceeding.  

---

## Data Organization  

IRIS requires a specific folder structure for proper execution. The dataset should be structured as follows:

- **Folder Name:** A four-character identifier, typically matching the **PDB ID** of the RNA receptor if applicable.
- The following files are required after **rDock** calculations:
  - `<folder_name>.pdb` (Include this if a PDB structure of the RNA receptor is available. If not, copy the `.mol2` file and rename the copied `.mol2` file as `<folder name>.pdb`.)
  - `<folder_name>.mol2` - receptor file without ligand bound
  - `<folder_name>_lig.sd` - ligand file
  - `<folder_name>_docking_out_sorted.sd` rDock docking output
  - `<folder_name>_cavity.log` - rDock cavity file 

Ensure that all required files are present before running IRIS.

---

## Running IRIS Scripts  

### 1. Generate Features for a Single RNA-Ligand Complex  
Navigate to the `IRIS/scripts` directory and run the following command:  
```
python get_features.py <folder_path> <pymol_path> <pymol_plugin_path>
```

### 2. Perform ML-Based Pose Re-Ranking via the following script:

- `IRIS.py`


#### Example Usage:  
```
python IRIS.py <folder_path>
```

Replace `<folder_path>` with the path to the directory containing the docking results and previously mentioned files.

---

## Notes  
- Ensure that all required dependencies and plugins are correctly installed before running the scripts.  
- Run all scripts from the `IRIS/scripts` directory.  
- If any errors occur, verify that the folder structure and input files match the expected format.  
