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

def calculate_polar_contacts(sdf_path, output_path, pymol_path, plugin_path, num_states):
    """
    Calculate polar contacts for each state using PyMOL and write results to an output file.
    """
    # Path to the PyMOL plugin for distance calculation
    plugin_full_path = os.path.join(plugin_path, "get_raw_distances.py")

    # PyMOL script to load the plugin, calculate distances, and iterate through states
    pymol_commands = f"""
# Load the get_raw_distances plugin
run {plugin_full_path};
load {sdf_path}, molecule;
# Add hydrogens to each state
cmd.h_add('molecule');

# Define donors and acceptors
select donors, (molecule and (elem N or elem O) and (neighbor hydro));
select acceptors, (molecule and (elem N or elem O) and !(neighbor hydro));

# Create a distance object to visualize the polar contacts
distance polar_contacts, donors, acceptors, mode=2, cutoff=3.6;
"""

    # Calculate polar contacts for each state
    for state in range(1, num_states + 1):
        pymol_commands += f"get_raw_distances polar_contacts, state={state}\n"

    # Write the PyMOL commands to a temporary script file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.pml') as tmp_script:
        tmp_script.write(pymol_commands)
        tmp_script_path = tmp_script.name

    # Construct the command to run PyMOL
    command = [pymol_path, '-cq', tmp_script_path]

    # Print command for debugging
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
    Process the PyMOL output file to count the number of polar contacts for each ligand pose.
    """
    with open(output_path, 'r') as file:
        lines = file.readlines()

    results = []
    current_state_contacts = 0
    last_state = 0

    for line in lines:
        if "get_raw_distances:" in line:
            # Increment contact count for each distance line
            current_state_contacts += 1

        state_line_match = re.search(r"state=(\d+)", line)
        if state_line_match:
            current_state = int(state_line_match.group(1))
            if current_state != last_state:
                if last_state != 0:
                    # Append the previous state's contact count to results
                    results.append(f"Ligand {last_state}: polar contacts: {current_state_contacts}\n")
                # Reset for new state
                current_state_contacts = 0
                last_state = current_state

    # Handle the last state after the loop ends
    if last_state != 0:
        results.append(f"Ligand {last_state}: polar contacts: {current_state_contacts}\n")

    # Write the final results to the output file
    with open(output_path, 'w') as f:
        f.writelines(results)

def process_directory_for_sdf_files(directory_path, pymol_path, plugin_path):
    """
    Process all SDF files in the specified directory for polar contact calculations.
    """
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith('out_sorted.sdf'):
                sdf_path = os.path.join(root, file)
                output_path = os.path.join(root, "polar_contacts.txt")
                print(f"Calculating polar contacts for {sdf_path} and writing results to {output_path}")

                # Calculate the number of states using RDKit
                num_states = get_number_of_states(sdf_path)

                # Calculate polar contacts for each state
                calculate_polar_contacts(sdf_path, output_path, pymol_path, plugin_path, num_states)
                process_output(output_path)

def main():
    if len(sys.argv) != 4:
        print("Usage: python polar.py <directory_path> <pymol_path> <plugin_path>")
        sys.exit(1)

    folder_path = sys.argv[1]
    pymol_path = sys.argv[2]
    plugin_path = sys.argv[3]

    # Call your process function with the provided arguments
    process_directory_for_sdf_files(folder_path, pymol_path, plugin_path)

if __name__ == "__main__":
    main()