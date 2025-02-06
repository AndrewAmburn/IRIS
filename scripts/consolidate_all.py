#!/usr/bin/env python
# coding: utf-8
import os
import sys 

# File names to read from the folder
file_names = [
    "hbonds.txt",
    "lig_pcd.txt",
    "msa_sasa.txt",
    "polar_contacts.txt",
    "scores_output.txt",
    "vdw_strain.txt"
]

# Function to read hbonds.txt data
def read_hbonds(file_path):
    data = {}
    with open(file_path) as f:
        for line in f:
            if "Ligand" in line:
                ligand, value = line.strip().split(": hydrogen bonds: ")
                data[ligand] = value
    return data

# Function to read lig_pcd.txt data
def read_lig_pcd(file_path):
    data = []
    with open(file_path) as f:
        header = next(f).strip().split("\t")  # Read and store the header
        for line in f:
            data.append(dict(zip(header, line.strip().split("\t"))))
    return data, header

# Function to read msa_sasa.txt data
def read_msa_sasa(file_path):
    data = {}
    with open(file_path) as f:
        ligand = ""
        for line in f:
            if "Ligand" in line:
                ligand = line.strip()[:-1]
                data[ligand] = {}
            elif "MSA" in line:
                data[ligand]["MSA"] = line.strip().split(": ")[1]
            elif "SASA" in line:
                data[ligand]["SASA"] = line.strip().split(": ")[1]
    return data

# Function to read polar_contacts.txt data
def read_polar_contacts(file_path):
    data = {}
    with open(file_path) as f:
        for line in f:
            if "Ligand" in line:
                ligand, value = line.strip().split(": polar contacts: ")
                data[ligand] = value
    return data

# Function to read scores_output.txt data
def read_scores_output(file_path):
    data = {}
    with open(file_path) as f:
        ligand = ""
        for line in f:
            if "Ligand" in line and "Scores" in line:
                ligand = line.strip()[:-8]
                data[ligand] = {}
            elif ": " in line:
                key, value = line.strip().split(": ")
                data[ligand][key] = value
    return data

# Function to read vdw_strain.txt data
def read_vdw_strain(file_path):
    data = {}
    with open(file_path) as f:
        for line in f:
            if "VDW Strain in state" in line:
                state, value = line.strip().split(": ")
                data[state] = value
    return data

def consolidate_descriptors(folder_path):
    # Check if all files are present before reading
    missing_files = [file for file in file_names if not os.path.exists(os.path.join(folder_path, file))]
    if missing_files:
        print(f"Error: Missing files in folder: {', '.join(missing_files)}")
        return
    
    # Read data from each file
    hbonds_data = read_hbonds(os.path.join(folder_path, "hbonds.txt"))
    lig_pcd_data, lig_pcd_headers = read_lig_pcd(os.path.join(folder_path, "lig_pcd.txt"))
    msa_sasa_data = read_msa_sasa(os.path.join(folder_path, "msa_sasa.txt"))
    polar_contacts_data = read_polar_contacts(os.path.join(folder_path, "polar_contacts.txt"))
    scores_output_data = read_scores_output(os.path.join(folder_path, "scores_output.txt"))
    vdw_strain_data = read_vdw_strain(os.path.join(folder_path, "vdw_strain.txt"))

    # Write consolidated data to descriptors.txt
    with open(os.path.join(folder_path, "descriptors.txt"), "w") as out_file:
        for i, ligand in enumerate(scores_output_data.keys(), 1):
            ligand_key = f"Ligand {i}"
            if scores_output_data.get(ligand_key):
                out_file.write(f"{ligand_key} Scores:\n")
                for key, value in scores_output_data[ligand_key].items():
                    out_file.write(f"  {key}: {value}\n")
                out_file.write(f"  PolarContactswithPocket: {hbonds_data.get(ligand_key, 'N/A')}\n")
                out_file.write(f"  IntraPolarContacts: {polar_contacts_data.get(ligand_key, 'N/A')}\n")
                out_file.write(f"  MSA: {msa_sasa_data.get(ligand_key, {}).get('MSA', 'N/A')}\n")
                out_file.write(f"  SASA: {msa_sasa_data.get(ligand_key, {}).get('SASA', 'N/A')}\n")
                out_file.write(f"  VDW_Strain: {vdw_strain_data.get(f'VDW Strain in state {i}', 'N/A')}\n")
                if i - 1 < len(lig_pcd_data):
                    for key, value in lig_pcd_data[i - 1].items():
                        out_file.write(f"  {key}: {value}\n")
                out_file.write("\n")

def main():
    # Ensure the correct number of arguments are passed
    if len(sys.argv) != 2:
        print("Usage: python <script_name> <directory_path>")
        sys.exit(1)

    # Get the folder path from command-line arguments
    folder_path = sys.argv[1]

    # Run the consolidation process
    consolidate_descriptors(folder_path)

if __name__ == "__main__":
    main()
