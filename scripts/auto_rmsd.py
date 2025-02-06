import os
import subprocess

# Define the parent directory
PARENT_DIR = "/Users/jasonamburn/Desktop/RNA_extra/RNA_complexes_TS_done_500/"
SCRIPT_PATH = "rmsd.py"

def has_sdf_files(folder_path):
    """
    Check if a folder contains any rmsd files.
    """
    for file in os.listdir(folder_path):
        if file.endswith("_rmsd_output.txt"):
            return True
    return False

def process_folders():
    """
    Iterate through subfolders, check for .sdf files, and run the script if none exist.
    """
    if not os.path.exists(PARENT_DIR):
        print(f"Directory {PARENT_DIR} does not exist. Please check the path.")
        return

    for subfolder in os.listdir(PARENT_DIR):
        subfolder_path = os.path.join(PARENT_DIR, subfolder)
        if os.path.isdir(subfolder_path):
            print(f"Checking folder: {subfolder}")
            if has_sdf_files(subfolder_path):
                print(f"Skipping {subfolder}: rmsd file(s) already exist.")
            else:
                print(f"No rmsd files found in {subfolder}. Running the script...")
                try:
                    subprocess.run(
                        ["python", SCRIPT_PATH, f"{PARENT_DIR}/{subfolder}"],
                        check=True
                    )
                    print(f"Script executed successfully for {subfolder}.")
                except subprocess.CalledProcessError as e:
                    print(f"Error running script for {subfolder}: {e}")

if __name__ == "__main__":
    process_folders()
