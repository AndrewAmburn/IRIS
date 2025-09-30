#!/usr/bin/env python
# coding: utf-8
import os
import sys
import time
from rdkit import Chem
from rdkit.Chem import Descriptors

# Function to calculate descriptors for a single molecule
def getMolDescriptors(mol, missingVal=None):
    res = {}
    try:
        Chem.SanitizeMol(mol)
    except Exception as e:
        print(f"[WARN] Molecule sanitation failed: {e}")
    try:
        Chem.GetSymmSSSR(mol)  # Initialize ring info
    except Exception as e:
        print(f"[WARN] Could not initialize ring info: {e}")
    for nm, fn in Descriptors._descList:
        try:
            val = fn(mol)
        except Exception as e:
            print(f"[ERROR] Descriptor {nm} failed: {e}")
            val = missingVal
        res[nm] = val
    return res

# Load SDF file with retry logic
def load_sdf_with_retry(sdf_file, retries=3, delay=2):
    for attempt in range(retries):
        try:
            suppl = Chem.SDMolSupplier(sdf_file, sanitize=False, removeHs=False)
            if suppl is not None:
                return suppl
        except Exception as e:
            print(f"Attempt {attempt + 1} failed to load {sdf_file}: {e}")
            time.sleep(delay)
    print(f"Failed to load SDF file after {retries} attempts: {sdf_file}")
    return None

# Combine descriptors from descriptors.txt and rdkit_full.txt
def combine_descriptors(descriptors_file, rdkit_file, output_file):
    combined_data = []

    with open(descriptors_file, 'r') as f:
        descriptors_data = f.readlines()

    with open(rdkit_file, 'r') as f:
        rdkit_data = f.readlines()

    ligand_number = None
    rdkit_dict = {}

    for line in rdkit_data:
        line = line.strip()
        if line.startswith("Ligand ") and "Scores" in line:
            ligand_number = line
            rdkit_dict[ligand_number] = []
        elif ligand_number:
            rdkit_dict[ligand_number].append(f"  {line}")

    for line in descriptors_data:
        line = line.strip()
        if line.startswith("Ligand ") and "Scores" in line:
            ligand_number = line
            combined_data.append(line)
            if ligand_number in rdkit_dict:
                combined_data.extend(rdkit_dict[ligand_number])
        else:
            combined_data.append(f"  {line}" if line and not line.startswith("Ligand ") else line)

    with open(output_file, 'w') as f:
        f.write('\n'.join(combined_data))

# Main processing function for a given subfolder
def process_directory(directory):
    folder_name = os.path.basename(os.path.normpath(directory))

    # Identify a valid out_sorted.sdf file
    sdf_file = None
    for fname in os.listdir(directory):
        if fname.endswith("out_sorted.sdf"):
            sdf_file = os.path.join(directory, fname)
            break

    if not sdf_file or not os.path.exists(sdf_file):
        print(f"[SKIP] No SDF file found in: {directory}")
        return

    rdkit_file = os.path.join(directory, f"{folder_name}_rdkit_full.txt")
    descriptors_file = os.path.join(directory, "descriptors.txt")
    combined_output_file = os.path.join(directory, "descriptors_full.txt")

    suppl = load_sdf_with_retry(sdf_file)
    if suppl is None:
        print(f"[ERROR] Failed to load molecules from {sdf_file}")
        return

    mols = [mol for mol in suppl if mol is not None]
    if not mols:
        print(f"[WARN] No valid molecules found in {sdf_file}")
        return

    output_lines = []
    ligand_counter = 1
    for mol in mols:
        output_lines.append(f"Ligand {ligand_counter} Scores:")
        descriptors = getMolDescriptors(mol, missingVal="NaN")
        for descriptor, value in descriptors.items():
            output_lines.append(f"  {descriptor}: {value}")
        ligand_counter += 1

    with open(rdkit_file, 'w') as f:
        f.write('\n'.join(output_lines))
    print(f"[OK] RDKit descriptors saved to: {rdkit_file}")

    if os.path.exists(descriptors_file) and os.path.isfile(rdkit_file):
        combine_descriptors(descriptors_file, rdkit_file, combined_output_file)
        print(f"[OK] Combined descriptors saved to: {combined_output_file}")
    else:
        print(f"[SKIP] Could not combine descriptors in: {directory}")

# Entry point for walking through a parent directory
def main():
    if len(sys.argv) != 2:
        print("Usage: python rdkit_full.py <parent_directory>")
        sys.exit(1)

    parent_dir = sys.argv[1]
    if not os.path.isdir(parent_dir):
        print(f"Error: {parent_dir} is not a valid directory.")
        sys.exit(1)

    for subfolder in os.listdir(parent_dir):
        subfolder_path = os.path.join(parent_dir, subfolder)
        if os.path.isdir(subfolder_path):
            print(f"\n📂 Processing: {subfolder_path}")
            process_directory(subfolder_path)

if __name__ == "__main__":
    main()