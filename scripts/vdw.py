import os
import subprocess
import tempfile
import sys

def calculate_vdw_strain(sdf_path, output_path, pymol_path, plugin_path):
    """
    Calculate VDW strain using the show_bumps PyMOL plugin and write results to an output file.
    """
    # Path to the show_bumps plugin
    plugin_full_path = os.path.join(plugin_path, "show_bumps.py")

    # PyMOL script to load the plugin and calculate VDW strain
    pymol_commands = f"""
# Load the show_bumps plugin
run {plugin_full_path};
load {sdf_path}, molecule;
show_bumps;
delete all;
"""

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

    # Handle output and write results to the output file
    if process.stdout:
        with open(output_path, 'w') as f:
            for line in process.stdout.split('\n'):
                if "VDW Strain" in line:  # Filter to write only VDW strain values
                    f.write(line + '\n')
    if process.stderr:
        print("STDERR:", process.stderr)

    # Clean up the temporary script file
    os.unlink(tmp_script_path)

def process_directory_for_sdf_files(directory_path, pymol_path, plugin_path):
    """
    Process all SDF files in the specified directory for VDW strain calculations.
    """
    for root, dirs, files in os.walk(directory_path):
        for file in files:
            if file.endswith('out_sorted.sdf'):
                sdf_path = os.path.join(root, file)
                output_path = os.path.join(root, "vdw_strain.txt")
                print(f"Calculating VDW strain for {sdf_path} and writing results to {output_path}")
                calculate_vdw_strain(sdf_path, output_path, pymol_path, plugin_path)

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
    process_directory_for_sdf_files(folder_path, pymol_path, plugin_path)

if __name__ == "__main__":
    main()
