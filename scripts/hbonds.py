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
    suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
    num_states = len([mol for mol in suppl if mol is not None])
    print(f"Number of states in {sdf_path}: {num_states}")
    return num_states

def calculate_hydrogen_bonds(rna_path, sdf_path, output_path, pymol_path, plugin_path, num_states):
    """
    Calculate hydrogen bonds between RNA and ligand for each state using PyMOL and write results to an output file.
    """
    print(f"Preparing to calculate hydrogen bonds for {rna_path} and {sdf_path} with {num_states} states.")
    
    # Path to the get_raw_distances plugin
    plugin_full_path = os.path.join(plugin_path, "get_raw_distances.py")
    
    if not os.path.exists(plugin_full_path):
        print(f"Error: Plugin not found at {plugin_full_path}")
        return

    # Consolidate the PyMOL commands into a single block for better execution
    pymol_commands = f"""
# Load the get_raw_distances plugin
run {plugin_full_path};
load {rna_path}, rna;
load {sdf_path}, lig;
cmd.h_add('rna');
cmd.h_add('lig');

# Define donors and acceptors for each state
select don_rna, (rna and (elem N or elem O) and (neighbor hydro));
select acc_rna, (rna and (elem O or (elem N and not (neighbor hydro))));
select don_lig, (lig and (elem N or elem O) and (neighbor hydro));
select acc_lig, (lig and (elem O or (elem N and not (neighbor hydro))));

# Create a distance object to visualize the hydrogen bonds
distance HBD, don_lig, acc_rna, mode=2, cutoff=3.2;
"""

    # Add the hydrogen bond calculations for each state
    for state in range(1, num_states + 1):
        pymol_commands += f"get_raw_distances HBD, state={state}\n"


    # Write the PyMOL commands to a temporary script file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.pml') as tmp_script:
        tmp_script.write(pymol_commands)
        tmp_script_path = tmp_script.name

    # Construct the command to run PyMOL
    command = [pymol_path, '-cq', tmp_script_path]
    
    # Print the full PyMOL command for debugging
    print(f"Running PyMOL command: {' '.join(command)}")
    
    # Execute the PyMOL command and capture the output
    process = subprocess.run(command, text=True, capture_output=True)

    # Write PyMOL output to a file for further processing
    if process.stdout:
        with open(output_path, 'w') as temp_file:
            temp_file.write(process.stdout)
    if process.stderr:
        print("STDERR:", process.stderr)

    # Clean up the temporary script file
    os.unlink(tmp_script_path)

def process_output(output_path):
    """
    Process the PyMOL output file to count the number of hydrogen bonds for each ligand pose.
    """
    with open(output_path, 'r') as file:
        lines = file.readlines()

    results = []
    current_state_bonds = 0
    last_state = 0

    for line in lines:
        if "get_raw_distances:" in line:
            # Increment bond count for each distance line
            current_state_bonds += 1

        state_line_match = re.search(r"state=(\d+)", line)
        if state_line_match:
            current_state = int(state_line_match.group(1))
            if current_state != last_state:
                if last_state != 0:
                    results.append(f"Ligand {last_state}: hydrogen bonds: {current_state_bonds}\n")
                current_state_bonds = 0  # Reset bond count for the new state
                last_state = current_state

    # Handle the last state after the loop ends
    if last_state != 0:
        results.append(f"Ligand {last_state}: hydrogen bonds: {current_state_bonds}\n")

    #output_hbonds_path = output_path.replace('output', 'hbonds')
    with open(output_path, 'w') as f:
        f.writelines(results)

def process_directory_for_files(directory_path, pymol_path, plugin_path):
    """
    Process all SDF files in the specified directory for hydrogen bond calculations.
    """
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith('out_sorted.sdf'):
                sdf_path = os.path.join(root, file)
                dir_name = os.path.basename(root)
                rna_path = os.path.join(root, f"{dir_name}.mol2")
                output_path = os.path.join(root, "hbonds.txt")
                print(f"Calculating hydrogen bonds for RNA {rna_path} and ligand {sdf_path}, writing results to {output_path}")

                # Calculate the number of states using RDKit
                num_states = get_number_of_states(sdf_path)

                # Calculate hydrogen bonds for each state
                calculate_hydrogen_bonds(rna_path, sdf_path, output_path, pymol_path, plugin_path, num_states)
                process_output(output_path)

def main():
    # Ensure the correct number of arguments are passed
    if len(sys.argv) != 4:
        print("Usage: python <script_name> <directory_path> <pymol_path> <plugin_path>")
        sys.exit(1)

    # Get the folder path, PyMOL path, and plugin directory path from command-line arguments
    folder_path = sys.argv[1]
    pymol_path = sys.argv[2]
    plugin_path = sys.argv[3]

    # Run the process
    process_directory_for_files(folder_path, pymol_path, plugin_path)

if __name__ == "__main__":
    main()
