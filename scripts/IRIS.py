#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import joblib
import pandas as pd
import numpy as np
from sklearn.impute import SimpleImputer
from rdkit import Chem
from rdkit.Chem import SDWriter


# ---------------------------
# Model loading
# ---------------------------
def load_model(model_path: str):
    model_path = os.path.abspath(model_path)

    if not os.path.exists(model_path):
        raise FileNotFoundError(f"Model file not found: {model_path}")

    if not os.path.isfile(model_path):
        raise IsADirectoryError(f"Model path is not a file: {model_path}")

    print(f"[IRIS] Loading model from: {model_path}")
    try:
        model = joblib.load(model_path)
    except Exception as e:
        raise RuntimeError(f"Failed to load model from {model_path}: {e}")

    return model


# ---------------------------
# Read descriptors_all.txt
# ---------------------------
def read_scores_file(file_path: str) -> pd.DataFrame:
    if not os.path.isfile(file_path):
        raise FileNotFoundError(f"Descriptor file not found: {file_path}")

    with open(file_path, "r") as f:
        lines = f.readlines()

    ligand_data_list = []
    ligand_data = {}

    for raw in lines:
        line = raw.strip()
        if not line:
            continue

        if line.startswith("Ligand ") and "Scores" in line:
            # Start of a new ligand block
            if ligand_data:
                ligand_data_list.append(ligand_data)
                ligand_data = {}
            continue

        # Parse "key: value"
        parts = line.split(":", 1)
        if len(parts) != 2:
            continue

        key, value = parts[0].strip(), parts[1].strip()
        try:
            ligand_data[key] = float(value)
        except ValueError:
            # ignore non-numeric fields
            pass

    if ligand_data:
        ligand_data_list.append(ligand_data)

    if not ligand_data_list:
        raise ValueError(f"No ligand descriptor blocks parsed from {file_path}")

    return pd.DataFrame(ligand_data_list)


# ---------------------------
# Preprocessing (KEEP feature names)
# ---------------------------
def preprocess_data(df_features: pd.DataFrame) -> pd.DataFrame:
    """
    Keep a pandas DataFrame all the way into model.predict() so sklearn
    transformers (e.g., StandardScaler) see valid feature names and don't warn.
    """
    imputer = SimpleImputer(strategy="most_frequent")

    # Fit/transform while preserving column names + index
    X = imputer.fit_transform(df_features)
    X_df = pd.DataFrame(X, columns=df_features.columns, index=df_features.index)

    return X_df


# ---------------------------
# Predict RMSD
# ---------------------------
def predict_rmsd(model, X) -> np.ndarray:
    return model.predict(X)


# ---------------------------
# SDF read / rank / write
# ---------------------------
def read_sdf_and_attach_ranks(
    sdf_path: str,
    df_complex: pd.DataFrame,
    sanitize: bool = False,
    score_field: str = "SCORE",
):
    if not os.path.isfile(sdf_path):
        raise FileNotFoundError(f"SDF file not found: {sdf_path}")

    if score_field not in df_complex.columns:
        raise ValueError(
            f"'{score_field}' column not found in descriptors dataframe; cannot map SDF poses."
        )

    # Build a lookup from SCORE -> Predicted_Rank.
    # If duplicates exist, keep the first occurrence.
    score_to_rank = {}
    for _, row in df_complex.iterrows():
        s = row.get(score_field, None)
        r = row.get("Predicted_Rank", None)
        if pd.isna(s) or pd.isna(r):
            continue
        if float(s) not in score_to_rank:
            score_to_rank[float(s)] = float(r)

    suppl = Chem.SDMolSupplier(sdf_path, sanitize=sanitize)
    mols = []

    for mol in suppl:
        if mol is None:
            continue

        if not mol.HasProp(score_field):
            continue

        try:
            score = float(mol.GetProp(score_field))
        except Exception:
            continue

        if score not in score_to_rank:
            continue

        mol.SetProp("Predicted_Rank", str(score_to_rank[score]))
        mols.append(mol)

    # Sort by Predicted_Rank (ascending)
    mols.sort(key=lambda m: float(m.GetProp("Predicted_Rank")))
    return mols


def write_sorted_sdf(sorted_mols, output_sdf_path: str):
    writer = SDWriter(output_sdf_path)
    try:
        for mol in sorted_mols:
            writer.write(mol)
    finally:
        writer.close()


# ---------------------------
# Main
# ---------------------------
def main():
    if len(sys.argv) not in (2, 3):
        print("Usage: python IRIS.py <path_to_folder> [optional_model_path]")
        sys.exit(1)

    directory = os.path.abspath(sys.argv[1])
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        sys.exit(1)

    folder_name = os.path.basename(os.path.normpath(directory))

    # Paths
    script_dir = os.path.dirname(os.path.realpath(__file__))
    default_model_path = os.path.abspath(os.path.join(script_dir, "..", "IRIS_model", "IRIS.pkl"))
    model_path = os.path.abspath(sys.argv[2]) if len(sys.argv) == 3 else default_model_path

    descriptors_path = os.path.join(directory, "descriptors_all.txt")
    sdf_file_path = os.path.join(directory, f"{folder_name}_docking_out_sorted.sdf")
    output_sdf_path = os.path.join(directory, f"{folder_name}_IRIS_sorted.sdf")

    # Input checks
    if not os.path.isfile(descriptors_path):
        print(f"Error: {descriptors_path} does not exist.")
        sys.exit(1)

    if not os.path.isfile(sdf_file_path):
        print(f"Error: {sdf_file_path} does not exist.")
        sys.exit(1)

    # 1) Load model
    model = load_model(model_path)

    # 2) Load descriptors
    df = read_scores_file(descriptors_path)

    # 3) Determine which features the model expects
    feature_names = getattr(model, "feature_names_in_", None)
    if feature_names is None:
        feature_names = list(df.columns)

    required_features = list(feature_names)

    # 4) Ensure all required feature columns exist
    for feat in required_features:
        if feat not in df.columns:
            df[feat] = 0.0

    # IMPORTANT: build feature matrix with correct column order + names
    df_features = (
        df[required_features]
        .replace([np.inf, -np.inf], 0.0)
        .fillna(0.0)
    )

    # 5) Preprocess and predict (keep DataFrame into predict to avoid sklearn warning)
    X_input = preprocess_data(df_features)   # DataFrame with column names
    predicted_rmsd = predict_rmsd(model, X_input)

    # 6) Add predictions and ranks
    df["Predicted_RMSD"] = predicted_rmsd
    df["Predicted_Rank"] = df["Predicted_RMSD"].rank(method="min", ascending=True)

    # 7) Reorder SDF and write output
    try:
        sorted_mols = read_sdf_and_attach_ranks(
            sdf_file_path, df, sanitize=False, score_field="SCORE"
        )
    except Exception as e:
        print(f"Error while mapping SDF poses to descriptor scores: {e}")
        sys.exit(1)

    if not sorted_mols:
        print(
            "Error: No poses could be mapped between SDF and descriptors.\n"
            "Most common causes:\n"
            "  - SCORE values in the SDF do not exactly match SCORE values in descriptors_all.txt\n"
            "  - SCORE is missing from one of the files\n"
        )
        sys.exit(1)

    write_sorted_sdf(sorted_mols, output_sdf_path)
    print(f"Processed and wrote sorted SDF file to {output_sdf_path}")


if __name__ == "__main__":
    main()
