#!/usr/bin/env python
# coding: utf-8
import os
import sys 

file_names = [
    "hbonds.txt",
    "lig_pcd.txt",
    "msa_sasa.txt",
    "polar_contacts.txt",
    "scores_output.txt",
    "vdw_strain.txt"
]

def read_hbonds(file_path):
    data = {}
    with open(file_path) as f:
        for line in f:
            if "Ligand" in line:
                ligand, value = line.strip().split(": hydrogen bonds: ")
                data[ligand] = value
    return data

def read_lig_pcd(file_path):
    data = []
    with open(file_path) as f:
        header = next(f).strip().split("\t")
        for line in f:
            data.append(dict(zip(header, line.strip().split("\t"))))
    return data, header

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

def read_polar_contacts(file_path):
    data = {}
    with open(file_path) as f:
        for line in f:
            if "Ligand" in line:
                ligand, value = line.strip().split(": polar contacts: ")
                data[ligand] = value
    return data

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

def read_vdw_strain(file_path):
    data = {}
    with open(file_path) as f:
        for line in f:
            if "VDW Strain in state" in line:
                state, value = line.strip().split(": ")
                data[state] = value
    return data

def consolidate_descriptors(folder_path):
    missing_files = [file for file in file_names if not os.path.exists(os.path.join(folder_path, file))]
    if missing_files:
        print(f"Warning: Missing files in folder: {', '.join(missing_files)}. Continuing with available data.")

    hbonds_data = read_hbonds(os.path.join(folder_path, "hbonds.txt")) if "hbonds.txt" not in missing_files else {}
    lig_pcd_data, lig_pcd_headers = read_lig_pcd(os.path.join(folder_path, "lig_pcd.txt")) if "lig_pcd.txt" not in missing_files else ([], [])
    msa_sasa_data = read_msa_sasa(os.path.join(folder_path, "msa_sasa.txt")) if "msa_sasa.txt" not in missing_files else {}
    polar_contacts_data = read_polar_contacts(os.path.join(folder_path, "polar_contacts.txt")) if "polar_contacts.txt" not in missing_files else {}
    scores_output_data = read_scores_output(os.path.join(folder_path, "scores_output.txt")) if "scores_output.txt" not in missing_files else {}
    vdw_strain_data = read_vdw_strain(os.path.join(folder_path, "vdw_strain.txt")) if "vdw_strain.txt" not in missing_files else {}

    with open(os.path.join(folder_path, "descriptors.txt"), "w") as out_file:
        ligand_keys = list(scores_output_data.keys()) if scores_output_data else ["Ligand 1"]
        for i, ligand_key in enumerate(ligand_keys, 1):
            out_file.write(f"{ligand_key} Scores:\n")
            for key, value in scores_output_data.get(ligand_key, {}).items():
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
    if len(sys.argv) != 2:
        print("Usage: python <script_name> <parent_directory>")
        sys.exit(1)

    parent_dir = sys.argv[1]
    if not os.path.isdir(parent_dir):
        print(f"Error: {parent_dir} is not a valid directory.")
        sys.exit(1)

    for subfolder in os.listdir(parent_dir):
        subfolder_path = os.path.join(parent_dir, subfolder)
        if os.path.isdir(subfolder_path):
            print(f"📁 Processing: {subfolder_path}")
            consolidate_descriptors(subfolder_path)


if __name__ == "__main__":
    main()
