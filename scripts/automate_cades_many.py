import os
import sys
import subprocess
import time

def run_script(script_name, args):
    """
    Helper function to run a Python script with the given arguments.
    
    Args:
        script_name (str): The name of the Python script to execute.
        args (list): List of arguments to pass to the script.
    """
    try:
        print(f"Running: {script_name} with arguments: {' '.join(args)}")
        subprocess.run(["python", script_name] + args, check=True)
        print(f"Completed: {script_name}\n")
        time.sleep(3)  # 3-second delay to ensure file writing completion
    except subprocess.CalledProcessError as e:
        print(f"Error running {script_name}: {e}")
        sys.exit(1)

def main():
    """
    Automates the entire process iteratively for each subfolder in a parent directory.
    """
    # Check for required command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python automate_cades_many.py <parent_directory> <pymol_path> <pymol_plugin_path>")
        sys.exit(1)


    # Get command-line arguments
    parent_directory = sys.argv[1]
    pymol_path = sys.argv[2]
    pymol_plugin_path = sys.argv[3]

    # Validate directories
    if not os.path.isdir(parent_directory):
        print(f"Error: {parent_directory} is not a valid directory.")
        sys.exit(1)
    if not os.path.exists(pymol_path):
        print(f"Error: PyMOL executable not found at {pymol_path}.")
        sys.exit(1)
    if not os.path.exists(pymol_plugin_path):
        print(f"Error: PyMOL plugin directory not found at {pymol_plugin_path}.")
        sys.exit(1)

    # List of scripts to run in sequence with their arguments
    scripts = [
        ("sd2sdf.py", []),
        ("consolidate_scores.py", []),
        ("lig_sd2pdb.py", []),
        ("rmsd.py", []),
        ("consolidate_scores_and_rmsd.py", []),
        ("generate_descriptors.py", [pymol_path, pymol_plugin_path]),
        ("rdkit_full.py", []),
        ("cav_volume.py", []),
        ("distances.py", [])
    ]

    # Process each subfolder in the parent directory
    print("Starting the automated pipeline...\n")
    for subfolder in os.listdir(parent_directory):
        working_directory = os.path.join(parent_directory, subfolder)

        if not os.path.isdir(working_directory):
            continue  # Skip non-directory files

        # Check if the descriptors_all.txt file exists
        descriptors_file = os.path.join(working_directory, "descriptors_all.txt")
        if os.path.exists(descriptors_file):
            print(f"Skipping directory: {working_directory}. descriptors_all.txt already exists.\n")
            continue

        print(f"\nProcessing directory: {working_directory}\n")

        # Run each script for the current subfolder as the working directory
        for script_name, additional_args in scripts:
            script_args = [working_directory] + additional_args
            run_script(script_name, script_args)

    print("Pipeline completed successfully for all subfolders!")

if __name__ == "__main__":
    main()