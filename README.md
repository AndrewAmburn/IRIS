README.md

# Step 1. Download requirements.txt in a conda-environment (python 3.9 and above).

# Install Conda dependencies
conda install --yes --file ~/IRIS/requirements.txt -c conda-forge

# Pip installs
pip install spyrmsd==0.8.0 rdkit==2024.3.2 scikit-learn==1.3.2 joblib==1.3.2 pandas 



# In addition to the requirements in the requirements.txt file, you will need to download a working version of PyMol, as well as the PyMol plugins: show_bumps.py and get_raw_distances.py

# Step 2. Structuring the data

# IRIS requires a specifc data structure for seamless script execution. The structure should follow this general flow:
- Folder name is a four character ID, usually corresponding to the same 4 character PDB ID of the RNA receptor if relevant
- The following files are required and should be available post-rDock calculations:
    - (folder name).pdb (Note: this file is intended for inclusion if there is a pdb of the corresponding RNA receptor. If no PDB exists of the RNA receptor, please copy the .mol2 and rename the extension as .pdb)
    - (folder name).mol2
    - (folder name)_lig.sd
    - (folder name)_docking_out_sorted.sd
    - (folder name)_cavity.log

# Step 3. Running the scripts

Note: run all scripts from IRIS/scripts directory


# Generate features for a single RNA-ligand complex:
1. get_features.py 

# Usage: python get_features.py <folder_path> <pymol_path> <pymol_plugin_path>"

2. IRIS_(docking_method).py 

# Valid IRIS re-ranking scripts: IRIS_RL_dock.py, IRIS_RL_dock_solv.py, IRIS_TS_dock.py, IRIS_TS_dock_solv.py 

# Usage example: python IRIS_RL_dock.py <folder_path>


#####



