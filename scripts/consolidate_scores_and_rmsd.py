#!/usr/bin/env python
# coding: utf-8

import os
import re
import sys

def extract_scores(ligand_block):
    """Extract scores from an SDF ligand block, excluding unwanted fields."""
    score_matches = re.finditer(r'>\s+<(\w+(?:\.\w+)*)>\s+(-?\d+\.\d+|\d+)', ligand_block)
    unwanted_fields = {'CHROM.0', 'CHROM.1', 'Name', 'Rbt.Receptor'}  # Define unwanted fields
    scores_dict = {}
    for match in score_matches:
        field_name = match.group(1)
        if field_name not in unwanted_fields:  # Only add if it's not an unwanted field
            scores_dict[field_name] = float(match.group(2))
    return scores_dict

def associate_rmsd_with_scores(scores, rmsd_file_path, ligand_index):
    """Associate RMSD values from the RMSD file with the ligand scores."""
    with open(rmsd_file_path, 'r') as rmsd_file:
        rmsd_lines = rmsd_file.readlines()

    # Ensure ligand_index is within the valid range
    if 0 <= ligand_index < len(rmsd_lines):
        rmsd_match = re.search(f'RMSD for molecule {ligand_index}: (\d+\.\d+)', rmsd_lines[ligand_index])
        if rmsd_match:
            scores['RMSD'] = float(rmsd_match.group(1))

def process_directory(directory_path):
    """Process the specified directory to extract scores and associate RMSD values."""
    folder_name = os.path.basename(directory_path)
    sdf_file_path = os.path.join(directory_path, f"{folder_name}_docking_out_sorted.sdf")
    rmsd_file_path = os.path.join(directory_path, f"{folder_name}_rmsd_output.txt")

    # Check if required files exist
    if not os.path.exists(sdf_file_path):
        print(f"Missing SDF file: {sdf_file_path}")
        return
    if not os.path.exists(rmsd_file_path):
        print(f"Missing RMSD file: {rmsd_file_path}")
        return

    # Read and split the SDF content
    with open(sdf_file_path, 'r') as file:
        sdf_content = file.read()
    ligand_blocks = re.split(r'\$\$\$\$', sdf_content.strip())

    output_file_path = os.path.join(directory_path, "scores_output.txt")

    # Process each ligand block and write results
    with open(output_file_path, 'w') as output_file:
        for i, ligand_block in enumerate(ligand_blocks):
            scores_dict = extract_scores(ligand_block)
            associate_rmsd_with_scores(scores_dict, rmsd_file_path, i)
            output_file.write(f'Ligand {i + 1} Scores:\n')
            for key, value in scores_dict.items():
                output_file.write(f'  {key}: {value}\n')
            output_file.write('\n')
    print(f"Scores and RMSD written to {output_file_path}")

def main():
    # Ensure the correct number of arguments are passed
    if len(sys.argv) != 2:
        print("Usage: python consolidate_scores_and_rmsd.py <directory_path>")
        sys.exit(1)

    # Get the directory from the command-line arguments
    directory_path = sys.argv[1]

    # Check if the input is a valid directory
    if os.path.isdir(directory_path):
        process_directory(directory_path)
    else:
        print(f"Error: {directory_path} is not a valid directory.")
        sys.exit(1)

if __name__ == "__main__":
    main()

