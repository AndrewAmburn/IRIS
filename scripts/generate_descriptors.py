import os
import subprocess
import sys

# Function to run a Python script
def run_script(script_path, args):
    """
    Run a Python script with the provided arguments.
    Args:
        script_path (str): Path to the Python script to execute.
        args (list): Arguments to pass to the script.
    """
    try:
        command = ["python", script_path] + args
        subprocess.run(command, check=True)
        print(f"Successfully ran: {script_path}")
    except subprocess.CalledProcessError as e:
        print(f"Error: {script_path} failed with return code {e.returncode}")
        sys.exit(1)

def main():
    """
    Main function to generate descriptors.
    Expects 3 command-line arguments:
        1. working_directory
        2. pymol_path
        3. pymol_plugin_path
    """
    if len(sys.argv) != 4:
        print("Usage: python generate_descriptors.py <working_directory> <pymol_path> <pymol_plugin_path>")
        sys.exit(1)

    # Get command-line arguments
    directory_path = sys.argv[1]
    pymol_path = sys.argv[2]
    pymol_plugin_path = sys.argv[3]

    # Validate the provided paths
    if not os.path.isdir(directory_path):
        print(f"Error: {directory_path} is not a valid directory.")
        sys.exit(1)

    if not os.path.isfile(pymol_path):
        print(f"Error: {pymol_path} is not a valid PyMOL application path.")
        sys.exit(1)

    if not os.path.isdir(pymol_plugin_path):
        print(f"Error: {pymol_plugin_path} is not a valid PyMOL plugin directory path.")
        sys.exit(1)

    # Define script paths (assuming all scripts are in the same directory as this script)
    scripts_dir = os.path.dirname(__file__)

    # List of scripts in the order they need to be run
    script_list = [
        "rdkit_descriptors.py",
        "sasa_msa.py",
        "polar.py",
        "vdw.py",
        "hbonds.py",
        "consolidate_all.py",
        "get_nonpolarcontacts.py",
        "rdkit_3d_descriptors.py"
    ]

    # Run each script in sequence with appropriate arguments
    for script_name in script_list:
        script_path = os.path.join(scripts_dir, script_name)

        if script_name in ["sasa_msa.py"]:
            # Pass 2 arguments: working_directory and pymol_path
            run_script(script_path, [directory_path, pymol_path])
        elif script_name in ["polar.py", "vdw.py", "hbonds.py"]:
            # Pass 3 arguments: working_directory, pymol_path, pymol_plugin_path
            run_script(script_path, [directory_path, pymol_path, pymol_plugin_path])
        else:
            # Pass only the working directory
            run_script(script_path, [directory_path])

    print("Descriptors generated successfully.")

if __name__ == "__main__":
    main()
