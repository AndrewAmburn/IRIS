import os
import sys
import subprocess
import time

def run_script(script_name, args):
    """
    Helper function to run a Python script with the given arguments silently,
    except for completion messages and errors.
    
    Args:
        script_name (str): The name of the Python script to execute.
        args (list): List of arguments to pass to the script.
    """
    try:
        subprocess.run(
            ["python", script_name] + args,
            check=True,
            stdout=subprocess.DEVNULL,  # Suppress standard output
            stderr=subprocess.DEVNULL   # Suppress standard error
        )
        time.sleep(3)  # 3-second delay to ensure file writing completion
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")
        sys.exit(1)

def main():
    """
    Automates the entire process for one complex.
    """
    # Check for required command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python get_features.py <working_directory> <pymol_path> <pymol_plugin_path>")
        sys.exit(1)

    # Get command-line arguments
    working_directory = sys.argv[1]
    pymol_path = sys.argv[2]
    pymol_plugin_path = sys.argv[3]

    # Validate directories
    if not os.path.isdir(working_directory):
        print(f"Error: {working_directory} is not a valid directory.")
        sys.exit(1)
    if not os.path.exists(pymol_path):
        print(f"Error: PyMOL executable not found at {pymol_path}.")
        sys.exit(1)
    if not os.path.exists(pymol_plugin_path):
        print(f"Error: PyMOL plugin directory not found at {pymol_plugin_path}.")
        sys.exit(1)

    # List of scripts to run in sequence with their arguments
    scripts = [
        ("sd2sdf.py", [working_directory]),
        ("consolidate_scores.py", [working_directory]),
        ("lig_sd2pdb.py", [working_directory]),
        ("consolidate_scores_and_rmsd.py", [working_directory]),
        ("generate_descriptors.py", [working_directory, pymol_path, pymol_plugin_path]),
        ("rdkit_full.py", [working_directory]),
        ("cav_volume.py", [working_directory]),
        ("distances.py", [working_directory]),
        ("receptor_descriptors.py", [working_directory])
    ]

    # Run each script in sequence
    print("Starting feature generation...\n")
    for script_name, args in scripts:
        run_script(script_name, args)

    print("Features generated successfully!")

if __name__ == "__main__":
    main()

