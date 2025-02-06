#!/usr/bin/env python
# coding: utf-8

import os
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors
import sys

# Function to calculate 3D descriptors for a ligand using rdMolDescriptors
def calculate_3d_descriptors(ligand):
    try:
        descriptors = {
            'PMI1': rdMolDescriptors.CalcPMI1(ligand),
            'PMI2': rdMolDescriptors.CalcPMI2(ligand),
            'PMI3': rdMolDescriptors.CalcPMI3(ligand),
            'NPR1': rdMolDescriptors.CalcNPR1(ligand),
            'NPR2': rdMolDescriptors.CalcNPR2(ligand),
            'RadiusOfGyration': rdMolDescriptors.CalcRadiusOfGyration(ligand),
            'InertialShapeFactor': rdMolDescriptors.CalcInertialShapeFactor(ligand),
            'Eccentricity': rdMolDescriptors.CalcEccentricity(ligand),
            'Asphericity': rdMolDescriptors.CalcAsphericity(ligand),
            'SpherocityIndex': rdMolDescriptors.CalcSpherocityIndex(ligand),
            'PBF': rdMolDescriptors.CalcPBF(ligand)
        }
        return descriptors
    except ValueError as e:
        print(f"Error calculating 3D descriptors: {e}")
        return {}

# Process the folder to calculate 3D descriptors
def process_directory(directory):
    output_file_path = os.path.join(directory, 'rdkit_3D_descriptors.txt')

    # Check if the file already exists
    if os.path.exists(output_file_path):
        print(f"{output_file_path} already exists. Skipping...")
        return

    sdf_files = [f for f in os.listdir(directory) if f.endswith('out_sorted.sdf')]
    descriptors = []
    for sdf_file in sdf_files:
        sdf_path = os.path.join(directory, sdf_file)
        suppl = Chem.SDMolSupplier(sdf_path, sanitize=False, removeHs=False)
        for ligand in suppl:
            if ligand:
                try:
                    ligand = Chem.AddHs(ligand)
                    descriptor_values = calculate_3d_descriptors(ligand)
                    descriptors.append(descriptor_values)
                except Exception as e:
                    print(f"Error processing ligand: {e}")
                    continue

    if descriptors:
        write_descriptors_to_file(descriptors, output_file_path)

# Write the calculated 3D descriptors to a file
def write_descriptors_to_file(descriptor_data, output_file):
    if not descriptor_data:
        return

    headers = descriptor_data[0].keys()
    with open(output_file, 'w') as f:
        f.write('\t'.join(headers) + '\n')
        for data in descriptor_data:
            line = '\t'.join([str(data[h]) for h in headers])
            f.write(line + '\n')

# Read the 3D descriptors from the file
def read_3d_descriptors(file_path):
    with open(file_path, 'r') as f:
        headers = f.readline().strip().split('\t')
        descriptors_list = []
        for line in f:
            values = line.strip().split('\t')
            descriptors_list.append(dict(zip(headers, values)))
    return descriptors_list

# Append the 3D descriptors to the descriptors.txt file
def append_3d_descriptors_to_file(descriptors_file, descriptors_data):
    with open(descriptors_file, 'r') as f:
        lines = f.readlines()

    keywords = ['PMI1:', 'PMI2:', 'PMI3:', 'NPR2:', 'RadiusOfGyration:', 'InertialShapeFactor:', 
                'Eccentricity:', 'SpherocityIndex:', 'Asphericity:', 'NPR1:', 'PBF:']
    
    existing_keywords = any(any(keyword in line for keyword in keywords) for line in lines)
    
    if existing_keywords:
        print(f"Skipping append: Descriptor lines already exist in {descriptors_file}")
        return

    ligand_index = 0
    with open(descriptors_file, 'w') as f:
        for line in lines:
            f.write(line)
            if line.startswith("Ligand") and "Scores:" in line:
                if ligand_index < len(descriptors_data):
                    descriptor_values = descriptors_data[ligand_index]
                    for key, value in descriptor_values.items():
                        f.write(f"  {key}: {value}\n")
                    ligand_index += 1

# Process the folder to append the 3D descriptors to descriptors.txt
def process_3d_append(directory):
    descriptors_file_path = os.path.join(directory, 'descriptors.txt')
    rdkit_3d_descriptors_file_path = os.path.join(directory, 'rdkit_3D_descriptors.txt')

    if os.path.exists(descriptors_file_path) and os.path.exists(rdkit_3d_descriptors_file_path):
        descriptors_data = read_3d_descriptors(rdkit_3d_descriptors_file_path)
        append_3d_descriptors_to_file(descriptors_file_path, descriptors_data)
    else:
        print(f"Files missing in directory {directory}: Skipping...")

def main():
    # Ensure the correct number of arguments are passed
    if len(sys.argv) != 2:
        print("Usage: python <script_name> <directory_path>")
        sys.exit(1)

    # Get the folder path from command-line arguments
    folder_path = sys.argv[1]

    # Step 1: Calculate 3D descriptors and write them to rdkit_3D_descriptors.txt
    process_directory(folder_path)
    print("3D descriptors extraction and file writing complete.")

    # Step 2: Append the 3D descriptors to descriptors.txt
    process_3d_append(folder_path)
    print("Appended 3D descriptors to descriptors.txt.")

if __name__ == "__main__":
    main()