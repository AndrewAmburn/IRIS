import os
import pickle
import pandas as pd
from sklearn.impute import SimpleImputer
import numpy as np
from rdkit import Chem
from rdkit.Chem import SDWriter
import sys

# Load the pickled model
def load_model(model_path):
    with open(model_path, 'rb') as model_file:
        model = pickle.load(model_file)
    return model

# Read the descriptors.txt file
def read_scores_file(file_path):
    with open(file_path, 'r') as f:
        lines = f.readlines()
    ligand_data_list = []
    ligand_data = {}
    for line in lines:
        if line.startswith('Ligand'):
            if ligand_data:
                ligand_data_list.append(ligand_data)
                ligand_data = {}
        else:
            parts = line.split(':')
            if len(parts) == 2:
                key, value = parts[0].strip(), parts[1].strip()
                try:
                    ligand_data[key] = float(value)
                except ValueError:
                    pass  # Handle non-numeric data
    if ligand_data:
        ligand_data_list.append(ligand_data)
    return pd.DataFrame(ligand_data_list)

# Preprocess the data (imputation)
def preprocess_data(df):
    imputer = SimpleImputer(strategy='most_frequent')
    df_imputed = imputer.fit_transform(df)
    return df_imputed

# Predict RMSD using the loaded model
def predict_rmsd(model, data):
    return model.predict(data)

# Read SDF, match scores, and assign predicted ranks
def read_sdf_and_match_scores(sdf_path, df_complex, sanitize=False):
    suppl = Chem.SDMolSupplier(sdf_path, sanitize=sanitize)
    mols = []
    for mol in suppl:
        if mol is not None:
            # Extract the SCORE property from the SDF
            score = float(mol.GetProp("SCORE"))
            # Find the row in the DataFrame that matches this SCORE
            matched_row = df_complex[df_complex['SCORE'].eq(score)]
            if not matched_row.empty:
                # If a match is found, store the Predicted_Rank in the molecule object
                mol.SetProp("Predicted_Rank", str(matched_row.iloc[0]['Predicted_Rank']))
                mols.append(mol)
    # Sort the molecules based on the Predicted_Rank
    sorted_mols = sorted(mols, key=lambda x: float(x.GetProp("Predicted_Rank")))
    return sorted_mols

# Write sorted molecules to a new SDF file
def write_sorted_sdf(sorted_mols, output_sdf_path):
    writer = SDWriter(output_sdf_path)
    for mol in sorted_mols:
        writer.write(mol)
    writer.close()

def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <path_to_folder>")
        sys.exit(1)

    directory = sys.argv[1]  # Get input path from command-line argument

    # Get the folder name (this will be used to construct the SDF file name)
    folder_name = os.path.basename(os.path.normpath(directory))


    # Set paths
    script_dir = os.path.dirname(os.path.realpath(__file__))  # Get the directory where the script is located
    model_path = os.path.join(script_dir, '../IRIS_models/TS_dock_solv.pkl')  # Model located in IRIS_models relative to this script
    descriptors_path = os.path.join(directory, 'descriptors_all.txt')
    sdf_file_path = os.path.join(directory, f"{folder_name}_docking_out_sorted.sdf")  # Construct SDF file name based on folder name
    output_sdf_path = os.path.join(directory, f"{folder_name}_IRIS_sorted.sdf")

    # Check if the provided directory exists
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        return

    # Check if the SDF file exists
    if not os.path.isfile(sdf_file_path):
        print(f"Error: {sdf_file_path} does not exist.")
        return

    # Step 1: Load the model
    model = load_model(model_path)

    # Step 2: Load the descriptors
    df = read_scores_file(descriptors_path)

    # Remove unnecessary columns (adjust this based on your model's input)
    required_features = ['MaxPartialCharge', 'SCORE.INTRA.POLAR.0', 'fr_lactone', 'fr_ArN', 'fr_Ar_NH', 
                        'SMR_VSA6', 'SlogP_VSA5', 'fr_Imine', 'PEOE_VSA7', 'BCUT2D_MWLOW', 'fr_bicyclic', 
                        'SCORE.INTRA.POLAR', 'RadiusOfGyration', 'fr_imide', 'com_C-C_mindist', 
                        'SCORE.INTRA.norm', 'fr_Ndealkylation1', 'SlogP_VSA8', 'SlogP_VSA3', 'VDW_Strain', 
                        'NumAliphaticCarbocycles', 'EState_VSA5', 'com_ring_to_ring_min_dist', 'VSA_EState6', 
                        'ClosestPolarContact', 'IntraPolarContacts', 'VSA_EState1', 'SpherocityIndex', 
                        'lig_C-C_min_dist', 'fr_guanido', 'MaxAbsEStateIndex', 'LigandHydrophobicContacts', 
                        'SCORE.SYSTEM.POLAR', 'NumAliphaticHeterocycles', 'fr_Al_COO', 'SlogP_VSA10', 
                        'fr_Nhpyrrole', 'EState_VSA2', 'RI', 'ligand_to_metal_min_dist', 'fr_sulfide', 
                        'LogP', 'SlogP_VSA12', 'lig_O-O_min_dist', 'fr_unbrch_alkane', 'SCORE.INTRA.REPUL.0', 
                        'NumHeteroatoms', 'receptor_to_metal_min_dist', 'NumAliphaticRings', 
                        'InertialShapeFactor', 'PBF', 'SCORE.INTER.VDW', 'VSA_EState7', 'SCORE.INTER.norm', 
                        'com_C-O_mindist', 'fr_allylic_oxid', 'fr_para_hydroxylation', 'fr_NH1', 
                        'SCORE.INTRA.VDW.0', 'PEOE_VSA9', 'SCORE.RESTR.norm', 'SCORE.INTRA.DIHEDRAL', 
                        'fr_piperzine', 'PEOE_VSA11', 'fr_NH0', 'SCORE.RESTR.CAVITY', 'lig_N-O_min_dist', 
                        'com_N-O_mindist', 'MolLogP', 'lig_N-N_min_dist', 'SCORE.RESTR', 
                        'SCORE.SYSTEM.DIHEDRAL', 'PEOE_VSA4', 'NumSaturatedCarbocycles', 'PEOE_VSA8', 
                        'com_N-N_mindist', 'RotatableBonds', 'PEOE_VSA10', 'MaxEStateIndex', 'fr_Ar_COO']

    # Ensure all expected features exist, filling missing ones with 0
    for feature in required_features:
        if feature not in df.columns:
            df[feature] = 0  # Fill missing columns with 0

    # Fill any NaN values with 0
    df.fillna(0, inplace=True)

    # Select the required features
    df_cleaned = df[required_features]
    # Step 3: Preprocess the data
    X_input = preprocess_data(df_cleaned)

    # Step 4: Make predictions
    predicted_rmsd = predict_rmsd(model, X_input)

    # Add the predictions to the DataFrame
    df['Predicted_Rank'] = predicted_rmsd

    # Step 5: Read, match, and sort the SDF file
    sorted_mols = read_sdf_and_match_scores(sdf_file_path, df, sanitize=False)

    # Step 6: Write the sorted molecules to a new SDF file
    write_sorted_sdf(sorted_mols, output_sdf_path)

    print(f"Processed and wrote sorted SDF file to {output_sdf_path}")

if __name__ == "__main__":
    main()
