import os
import re

def extract_cavity_volume(log_file_path):
    """
    Extract the total volume from the cavity log file.

    Parameters:
        log_file_path (str): Path to the cavity log file.

    Returns:
        float: The extracted total volume.
    """
    if not os.path.exists(log_file_path):
        print(f"Warning: {log_file_path} does not exist. Skipping.")
        return None

    try:
        with open(log_file_path, 'r') as file:
            for line in file:
                match = re.search(r'Total volume ([0-9]+(?:\.[0-9]+)?) A\^3', line)
                if match:
                    return float(match.group(1))
    except Exception as e:
        print(f"Error reading {log_file_path}: {e}")

    print(f"No total volume found in {log_file_path}. Skipping.")
    return None

def update_scores_with_cavity_volume(directory_path):
    """
    Update the descriptors_full.txt file with the cavity volume from the log file.

    Parameters:
        directory_path (str): Path to the directory containing the descriptors_full.txt and cavity log file.
    """
    folder_name = os.path.basename(os.path.normpath(directory_path))
    log_file_path = os.path.join(directory_path, f"{folder_name}_cavity.log")
    scores_file_path = os.path.join(directory_path, "descriptors_full.txt")

    cavity_volume = extract_cavity_volume(log_file_path)
    if cavity_volume is None:
        return  # Skip if cavity volume couldn't be extracted

    if not os.path.exists(scores_file_path):
        print(f"Warning: {scores_file_path} does not exist. Skipping.")
        return

    try:
        with open(scores_file_path, 'r') as scores_file:
            scores_lines = scores_file.readlines()

        updated_lines = []
        for line in scores_lines:
            updated_lines.append(line)
            if line.strip().startswith("Ligand") and "Scores:" in line:
                updated_lines.append(f"  cav_volume: {cavity_volume}\n")

        with open(scores_file_path, 'w') as scores_file:
            scores_file.writelines(updated_lines)

        print(f"Updated {scores_file_path} with cavity volume for each ligand.")
    except Exception as e:
        print(f"Error updating {scores_file_path}: {e}")

if __name__ == "__main__":
    import sys

    if len(sys.argv) != 2:
        print("Usage: python cav_volume.py <directory_path>")
        sys.exit(1)

    directory_path = sys.argv[1]

    if os.path.isdir(directory_path):
        update_scores_with_cavity_volume(directory_path)
    else:
        print(f"Error: {directory_path} is not a valid directory.")
        sys.exit(1)
