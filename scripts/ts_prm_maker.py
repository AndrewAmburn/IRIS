import os
import sys
import subprocess

def replace_z_in_content(content, replacement):
    """
    Replace 'z' in the content with the specified replacement string.
    """
    return content.replace('z', replacement)

def get_first_coords(ligand_file, pymol_path):
    """
    Use PyMOL to extract the first set of coordinates from a ligand file.

    Args:
        ligand_file (str): Path to the ligand file.
        pymol_path (str): Path to the PyMOL executable.

    Returns:
        list: A list of x, y, z coordinates as strings.
    """
    temp_file_name = 'temp_coords.txt'
    command = f"{pymol_path} -cq -d 'load {ligand_file}; save {temp_file_name}, format=xyz'"
    subprocess.run(command, shell=True)

    # Read the coordinates from the temp file
    with open(temp_file_name, 'r') as temp_file:
        lines = temp_file.readlines()
        if len(lines) > 2:
            coords_line = lines[2].split()
            coords = coords_line[1:4]  # Get x, y, z coordinates
        else:
            raise ValueError("Unable to extract coordinates from ligand file.")
        
    os.remove(temp_file_name)  # Clean up the temp file
    return coords

def replace_center_in_content(content, coordinates):
    """
    Replace the 'CENTER' line in the content with new coordinates.

    Args:
        content (str): The file content to update.
        coordinates (list): List of x, y, z coordinates as strings.

    Returns:
        str: Updated content with the new CENTER line.
    """
    lines = content.split('\n')
    for i, line in enumerate(lines):
        if line.strip().startswith('CENTER'):
            leading_spaces = len(line) - len(line.lstrip())  # Preserve indentation
            lines[i] = ' ' * leading_spaces + f'CENTER ({coordinates[0]},{coordinates[1]},{coordinates[2]})'
            break
    return '\n'.join(lines)

def process_folders(parent_directory, pymol_path, parameter_file_path):
    """
    Process each subfolder in the parent directory to generate updated .prm files.

    Args:
        parent_directory (str): Path to the parent directory containing subfolders.
        pymol_path (str): Path to the PyMOL executable.
        parameter_file_path (str): Path to the template z_rdock.prm file.
    """
    # Read the template parameter file
    with open(parameter_file_path, 'r') as original_file:
        original_content = original_file.read()

    # Process each folder in the parent directory
    for folder_name in os.listdir(parent_directory):
        folder_path = os.path.join(parent_directory, folder_name)

        if os.path.isdir(folder_path):
            # Replace "z" with the folder name in the template
            content_with_folder_name = replace_z_in_content(original_content, folder_name)

            # Get the ligand file and extract coordinates
            ligand_file = os.path.join(folder_path, f"{folder_name}_lig.sd")
            if not os.path.exists(ligand_file):
                print(f"Warning: Ligand file {ligand_file} not found. Skipping folder {folder_name}.")
                continue

            try:
                coords = get_first_coords(ligand_file, pymol_path)
            except Exception as e:
                print(f"Error extracting coordinates for {folder_name}: {e}")
                continue

            # Replace CENTER line with the new coordinates
            new_content = replace_center_in_content(content_with_folder_name, coords)

            # Write the updated content to a new .prm file
            new_file_name = f'{folder_name}.prm'
            new_file_path = os.path.join(folder_path, new_file_name)
            with open(new_file_path, 'w') as new_file:
                new_file.write(new_content)

            print(f'Created file {new_file_path}')

if __name__ == "__main__":
    # Validate command-line arguments
    if len(sys.argv) != 4:
        print("Usage: python script.py <parent_directory> <pymol_path> <z_rdock.prm_path>")
        sys.exit(1)

    parent_directory = sys.argv[1]
    pymol_path = sys.argv[2]
    parameter_file_path = sys.argv[3]

    # Verify input paths
    if not os.path.isdir(parent_directory):
        print(f"Error: {parent_directory} is not a valid directory.")
        sys.exit(1)

    if not os.path.exists(pymol_path):
        print(f"Error: PyMOL executable not found at {pymol_path}.")
        sys.exit(1)

    if not os.path.exists(parameter_file_path):
        print(f"Error: Parameter file not found at {parameter_file_path}.")
        sys.exit(1)

    # Process the folders
    process_folders(parent_directory, pymol_path, parameter_file_path)
