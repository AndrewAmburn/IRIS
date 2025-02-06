import os
import sys
from rdkit import Chem
from rdkit.Chem import Descriptors
import time

# Function to calculate descriptors for a single molecule
def getMolDescriptors(mol, missingVal=None):
    """Calculate the full list of descriptors for a molecule.

    Args:
        mol: RDKit molecule object
        missingVal: Value to assign if descriptor calculation fails

    Returns:
        dict: Dictionary of descriptor names and values
    """
    res = {}
    for nm, fn in Descriptors._descList:
        try:
            val = fn(mol)
        except:
            val = missingVal
        res[nm] = val
    return res

# Function to load the SDF file with retry mechanism
def load_sdf_with_retry(sdf_file, retries=3, delay=2):
    """Attempt to load an SDF file with retries on failure.

    Args:
        sdf_file: Path to the SDF file
        retries: Number of retry attempts
        delay: Delay between retries in seconds

    Returns:
        Chem.SDMolSupplier or None if loading fails
    """
    for attempt in range(retries):
        try:
            suppl = Chem.SDMolSupplier(sdf_file)
            if suppl is not None:
                return suppl
        except Exception as e:
            print(f"Attempt {attempt + 1} failed to load {sdf_file}: {e}")
            time.sleep(delay)
    print(f"Failed to load SDF file after {retries} attempts: {sdf_file}")
    return None

# Function to combine descriptors from two files
def combine_descriptors(descriptors_file, rdkit_file, output_file):
    """Combine descriptors from descriptors.txt and rdkit_full.txt.

    Args:
        descriptors_file (str): Path to descriptors.txt file.
        rdkit_file (str): Path to rdkit_full.txt file.
        output_file (str): Path to save the combined descriptors.

    Returns:
        None
    """
    combined_data = []

    # Read descriptors.txt
    with open(descriptors_file, 'r') as f:
        descriptors_data = f.readlines()

    # Read rdkit_full.txt
    with open(rdkit_file, 'r') as f:
        rdkit_data = f.readlines()

    # Combine data by ligand number
    ligand_number = None
    rdkit_dict = {}

    # Parse rdkit_full.txt into a dictionary
    for line in rdkit_data:
        line = line.strip()
        if line.startswith("Ligand ") and "Scores" in line:
            ligand_number = line
            rdkit_dict[ligand_number] = []
        elif ligand_number:
            rdkit_dict[ligand_number].append(f"  {line}")  # Indent lines

    # Combine descriptors.txt with rdkit_full.txt
    for line in descriptors_data:
        line = line.strip()
        if line.startswith("Ligand ") and "Scores" in line:
            ligand_number = line
            combined_data.append(line)
            if ligand_number in rdkit_dict:
                combined_data.extend(rdkit_dict[ligand_number])
        else:
            combined_data.append(f"  {line}" if line and not line.startswith("Ligand ") else line)

    # Save combined data
    with open(output_file, 'w') as f:
        f.write('\n'.join(combined_data))

# Main function to process a single input directory
def process_directory(directory):
    folder_name = os.path.basename(directory)
    sdf_file = os.path.join(directory, f"{folder_name}_docking_out_sorted.sdf")

    # Check if the SDF file exists
    if os.path.exists(sdf_file):
        # Read the SDF file with retry
        suppl = load_sdf_with_retry(sdf_file)
        if suppl is None:
            return  # Skip if the file could not be loaded

        # Prepare to store descriptors
        output_lines = []

        ligand_counter = 1
        for mol in suppl:
            if mol is None:
                continue  # Skip invalid molecules

            # Add ligand header
            output_lines.append(f"Ligand {ligand_counter} Scores:")

            # Calculate descriptors for the molecule
            descriptors = getMolDescriptors(mol, missingVal=None)
            for descriptor, value in descriptors.items():
                output_lines.append(f"  {descriptor}: {value}")

            ligand_counter += 1

        # Save descriptors to a TXT file in the same directory
        rdkit_file = os.path.join(directory, f"{folder_name}_rdkit_full.txt")
        with open(rdkit_file, 'w') as f:
            f.write('\n'.join(output_lines))

        print(f"Descriptors saved to {rdkit_file}")

    # Combine descriptors if descriptors.txt exists
    descriptors_file = os.path.join(directory, "descriptors.txt")
    combined_output_file = os.path.join(directory, "descriptors_full.txt")

    if os.path.exists(descriptors_file):
        combine_descriptors(descriptors_file, rdkit_file, combined_output_file)
        print(f"Combined descriptors saved to {combined_output_file}")

# Command-line entry point
def main():
    if len(sys.argv) != 2:
        print("Usage: python script.py <directory>")
        sys.exit(1)

    input_directory = sys.argv[1]

    if not os.path.isdir(input_directory):
        print(f"Error: {input_directory} is not a valid directory.")
        sys.exit(1)

    process_directory(input_directory)

if __name__ == "__main__":
    main()
