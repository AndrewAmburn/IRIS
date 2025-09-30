#!/usr/bin/env python
# coding: utf-8

import os
import subprocess
import tempfile
import re
from rdkit import Chem
import sys 

def get_number_of_states(sdf_path):
    """
    Use RDKit to count the number of molecules (states) in the SDF file.
    """
    suppl = Chem.SDMolSupplier(sdf_path, sanitize=False, removeHs=False)
    num_states = len([mol for mol in suppl if mol is not None])
    print(f"Number of states in {sdf_path}: {num_states}")
    return num_states

def calculate_hydrogen_bonds(rna_path, sdf_path, output_path, pymol_path, plugin_path, num_states):
    """
    Calculate hydrogen bonds between RNA and ligand for each state using PyMOL and write results to an output file.
    """
    print(f"Preparing to calculate hydrogen bonds for {rna_path} and {sdf_path} with {num_states} states.")

    plugin_full_path = os.path.join(plugin_path, "get_raw_distances.py")
    if not os.path.exists(plugin_full_path):
        print(f"Error: Plugin not found at {plugin_full_path}")
        return

    pymol_commands = f"""
run {plugin_full_path};
load {rna_path}, rna;
load {sdf_path}, lig;
cmd.h_add('rna');
cmd.h_add('lig');

select don_rna, (rna and (elem N or elem O) and (neighbor hydro));
select acc_rna, (rna and (elem O or (elem N and not (neighbor hydro))));
select don_lig, (lig and (elem N or elem O) and (neighbor hydro));
select acc_lig, (lig and (elem O or (elem N and not (neighbor hydro))));

distance HBD, (don_lig, acc_rna) or (don_rna, acc_lig), mode=2, cutoff=3.2;
"""

    for state in range(1, num_states + 1):
        pymol_commands += f"get_raw_distances HBD, state={state}\n"

    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.pml') as tmp_script:
        tmp_script.write(pymol_commands)
        tmp_script_path = tmp_script.name

    command = [pymol_path, '-cq', tmp_script_path]
    print(f"Running PyMOL command: {' '.join(command)}")
    process = subprocess.run(command, text=True, capture_output=True)

    if process.stdout:
        with open(output_path, 'w') as temp_file:
            temp_file.write(process.stdout)
    if process.stderr:
        print("STDERR:", process.stderr)

    os.unlink(tmp_script_path)

def process_output(output_path):
    with open(output_path, 'r') as file:
        lines = file.readlines()

    results = []
    current_state_bonds = 0
    last_state = 0

    for line in lines:
        if "get_raw_distances:" in line:
            current_state_bonds += 1

        state_line_match = re.search(r"state=(\d+)", line)
        if state_line_match:
            current_state = int(state_line_match.group(1))
            if current_state != last_state:
                if last_state != 0:
                    results.append(f"Ligand {last_state}: hydrogen bonds: {current_state_bonds}\n")
                current_state_bonds = 0
                last_state = current_state

    if last_state != 0:
        results.append(f"Ligand {last_state}: hydrogen bonds: {current_state_bonds}\n")

    with open(output_path, 'w') as f:
        f.writelines(results)

    if os.path.getsize(output_path) == 0:
        print(f"⚠️ No hydrogen bonds found in: {output_path}")

def process_directory_for_files(directory_path, pymol_path, plugin_path):
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith('out_sorted.sdf'):
                sdf_path = os.path.join(root, file)
                dir_name = os.path.basename(root)
                rna_path = os.path.join(root, f"{dir_name}.mol2")
                output_path = os.path.join(root, "hbonds.txt")
                print(f"Calculating hydrogen bonds for RNA {rna_path} and ligand {sdf_path}, writing results to {output_path}")

                num_states = get_number_of_states(sdf_path)
                if num_states == 0:
                    print(f"⚠️ Skipping {sdf_path} — no valid ligand states detected.")
                    continue

                calculate_hydrogen_bonds(rna_path, sdf_path, output_path, pymol_path, plugin_path, num_states)
                process_output(output_path)

def main():
    if len(sys.argv) != 4:
        print("Usage: python <script_name> <directory_path> <pymol_path> <plugin_path>")
        sys.exit(1)

    folder_path = sys.argv[1]
    pymol_path = sys.argv[2]
    plugin_path = sys.argv[3]

    process_directory_for_files(folder_path, pymol_path, plugin_path)

if __name__ == "__main__":
    main()