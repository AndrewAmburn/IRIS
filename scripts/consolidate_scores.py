#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys

def extract_scores(ligand_block):
    score_matches = re.finditer(r'>\s+<(\w+(?:\.\w+)*)>\s+(-?\d+\.\d+|\d+)', ligand_block)
    unwanted_fields = {'CHROM.0', 'CHROM.1', 'Name', 'Rbt.Receptor'}
    scores_dict = {}
    for match in score_matches:
        field_name = match.group(1)
        if field_name not in unwanted_fields:
            scores_dict[field_name] = float(match.group(2))
    return scores_dict

def find_sdf_file(input_folder):
    # Look for any file that ends with "_docking_out_sorted.sdf" within the folder
    for file_name in os.listdir(input_folder):
        if file_name.endswith("_docking_out_sorted.sdf"):
            return os.path.join(input_folder, file_name)
    return None

def process_file(input_folder):
    # Find the correct SDF file
    sdf_file_path = find_sdf_file(input_folder)
    
    if not sdf_file_path:
        print(f"Error: No SDF file matching '_docking_out_sorted.sdf' found in {input_folder}.")
        return

    # Continue with the rest of the processing
    with open(sdf_file_path, 'r') as file:
        sdf_content = file.read()

    ligand_blocks = re.split(r'\$\$\$\$', sdf_content.strip())

    output_file_path = os.path.join(input_folder, "scores_output.txt")

    with open(output_file_path, 'w') as output_file:
        for i, ligand_block in enumerate(ligand_blocks):
            scores_dict = extract_scores(ligand_block)
            output_file.write(f'Ligand {i+1} Scores:\n')
            for key, value in scores_dict.items():
                output_file.write(f'  {key}: {value}\n')
            output_file.write('\n')



def main():
    if len(sys.argv) != 2:
        print("Usage: python consolidate_scores.py <directory_path>")
        sys.exit(1)

    input_folder = sys.argv[1]
    
    if os.path.isdir(input_folder):
        process_file(input_folder)
    else:
        print(f"Error: {input_folder} is not a valid directory.")

if __name__ == "__main__":
    main()

