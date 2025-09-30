import os
import sys
import subprocess
import tempfile
from rdkit import Chem

def get_number_of_states(sdf_path):
    suppl = Chem.SDMolSupplier(sdf_path, sanitize=False, removeHs=False)
    num_states = len([mol for mol in suppl if mol is not None])
    print(f"Number of states in {sdf_path}: {num_states}")
    return num_states

def calculate_surface_areas_with_command(sdf_path, output_path, pymol_path, max_states):
    """
    Calculate SASA and MSA using PyMOL and write results to an output file.
    """
    pymol_commands = f"""
load {sdf_path}, molecule
"""
    for state in range(1, max_states + 1):
        pymol_commands += f"""
set dot_solvent, off
cmd.set('dot_density', 3)
msa = cmd.get_area('all', state={state})
set dot_solvent, on
cmd.set('dot_density', 3)
sasa = cmd.get_area('all', state={state})
print('Ligand {state}:')
print('MSA:', msa)
print('SASA:', sasa)
"""

    # Write the pymol_commands to a temporary file
    with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.pml') as tmp_script:
        tmp_script.write(pymol_commands)
        tmp_script_path = tmp_script.name

    # Construct the final PyMOL command
    command = [pymol_path, '-cq', tmp_script_path]

    # Execute the PyMOL command and capture the output
    process = subprocess.run(command, text=True, capture_output=True)

    if process.returncode != 0:
        print(f"Error running PyMOL command: {process.stderr}")
        sys.exit(1)

    # Write the results to the output file
    if process.stdout:
        with open(output_path, 'w') as f:
            for line in process.stdout.split('\n'):
                if line.strip() and not line.strip().startswith('PyMOL>'):
                    if 'Ligand' in line or 'MSA:' in line or 'SASA:' in line:
                        f.write(line.strip() + '\n')

    # Clean up the temporary file
    os.unlink(tmp_script_path)

def process_directory_for_sdf_files(directory_path, pymol_path):
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith('sorted.sdf'):
                sdf_path = os.path.join(root, file)
                output_path = os.path.join(root, "msa_sasa.txt")
                print(f"Processing {sdf_path} and writing results to {output_path}")

                # Get the number of states
                num_states = get_number_of_states(sdf_path)

                # Calculate SASA and MSA
                calculate_surface_areas_with_command(sdf_path, output_path, pymol_path, num_states)

def main():
    # Ensure correct number of arguments are passed
    if len(sys.argv) != 3:
        print("Usage: python sasa_msa.py <directory_path> <pymol_path>")
        sys.exit(1)

    # Get the folder path and PyMOL path from command-line arguments
    folder_path = sys.argv[1]
    pymol_path = sys.argv[2]

    # Run the process
    process_directory_for_sdf_files(folder_path, pymol_path)

if __name__ == "__main__":
    main()
