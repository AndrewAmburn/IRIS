import os
import sys
import time
from rdkit import Chem
from rdkit.Chem import AllChem
import numpy as np
import time

def calculate_min_distances(mol1, mol2=None, atom_pairs=[('N', 'C')]):
    """
    Calculates the minimum distances between specified atom pairs.

    Parameters:
        mol1 (Mol): The first molecule (e.g., ligand).
        mol2 (Mol): The second molecule (e.g., receptor), or None if only within mol1.
        atom_pairs (list of tuples): List of atom pairs to calculate distances for, e.g., [('N', 'C')].

    Returns:
        dict: Minimum distances for each atom pair.
    """
    if mol2 is None:
        mol2 = mol1

    # Get 3D coordinates of atoms
    conf1 = mol1.GetConformer()
    coords1 = np.array([list(conf1.GetAtomPosition(i)) for i in range(mol1.GetNumAtoms())])

    conf2 = mol2.GetConformer()
    coords2 = np.array([list(conf2.GetAtomPosition(i)) for i in range(mol2.GetNumAtoms())])

    # Get atom types
    atom_types1 = [atom.GetSymbol() for atom in mol1.GetAtoms()]
    atom_types2 = [atom.GetSymbol() for atom in mol2.GetAtoms()]

    # Initialize distances dictionary
    distances = {}

    for atom_type1, atom_type2 in atom_pairs:
        # Get indices of the specified atom types
        indices1 = [i for i, atype in enumerate(atom_types1) if atype == atom_type1]
        indices2 = [i for i, atype in enumerate(atom_types2) if atype == atom_type2]

        # Skip calculation if either list is empty
        if not indices1 or not indices2:
            distances[f"{atom_type1}-{atom_type2}"] = 0
            continue

        # Expand indices to select relevant atom coordinates
        selected_coords1 = coords1[indices1]
        selected_coords2 = coords2[indices2]

        # Compute pairwise distances
        dist_matrix = np.linalg.norm(
            selected_coords1[:, None, :] - selected_coords2[None, :, :], axis=2
        )

        # Exclude self-pairing distances and covalently bonded atoms if mol1 == mol2
        if mol1 == mol2:
            for i, idx1 in enumerate(indices1):
                for j, idx2 in enumerate(indices2):
                    if idx1 == idx2 or mol1.GetBondBetweenAtoms(idx1, idx2):
                        dist_matrix[i, j] = np.inf

        # Store the minimum distance
        min_dist = np.min(dist_matrix)
        distances[f"{atom_type1}-{atom_type2}"] = min_dist if np.isfinite(min_dist) else 0

    return distances

def calculate_ring_atom_distances(mol1, mol2=None):
    """
    Calculates the minimum distance between ring atoms and other atoms, as well as ring-to-ring atom distances.

    Parameters:
        mol1 (Mol): The first molecule (e.g., ligand).
        mol2 (Mol): The second molecule (e.g., receptor), or None if only within mol1.

    Returns:
        dict: Minimum distances for ring-other atom and ring-ring atom distances.
    """
    if mol2 is None:
        mol2 = mol1

    # Get 3D coordinates of atoms
    conf1 = mol1.GetConformer()
    coords1 = np.array([list(conf1.GetAtomPosition(i)) for i in range(mol1.GetNumAtoms())])

    conf2 = mol2.GetConformer()
    coords2 = np.array([list(conf2.GetAtomPosition(i)) for i in range(mol2.GetNumAtoms())])

    # Get ring atom indices
    ring_atoms1 = [atom.GetIdx() for atom in mol1.GetAtoms() if atom.IsInRing()]
    other_atoms1 = [atom.GetIdx() for atom in mol1.GetAtoms() if not atom.IsInRing()]

    ring_atoms2 = [atom.GetIdx() for atom in mol2.GetAtoms() if atom.IsInRing()]
    other_atoms2 = [atom.GetIdx() for atom in mol2.GetAtoms() if not atom.IsInRing()]

    distances = {}

    # Calculate ring-to-other distances for mol1
    if ring_atoms1 and other_atoms1:
        ring_coords1 = coords1[ring_atoms1]
        other_coords1 = coords1[other_atoms1]
        dist_matrix1 = np.linalg.norm(
            ring_coords1[:, None, :] - other_coords1[None, :, :], axis=2
        )
        distances["lig_ring_min_dist"] = np.min(dist_matrix1)
    else:
        distances["lig_ring_min_dist"] = 0.0

    # Calculate ring-to-other distances between mol1 and mol2
    if ring_atoms1 and other_atoms2:
        ring_coords1 = coords1[ring_atoms1]
        other_coords2 = coords2[other_atoms2]
        dist_matrix2 = np.linalg.norm(
            ring_coords1[:, None, :] - other_coords2[None, :, :], axis=2
        )
        distances["com_ring_min_dist"] = np.min(dist_matrix2)
    else:
        distances["com_ring_min_dist"] = 0.0

    # Calculate ring-to-ring distances for mol1
    if ring_atoms1:
        ring_coords1 = coords1[ring_atoms1]
        dist_matrix_ring1 = np.linalg.norm(
            ring_coords1[:, None, :] - ring_coords1[None, :, :], axis=2
        )
        np.fill_diagonal(dist_matrix_ring1, np.inf)  # Exclude self-distances
        distances["lig_ring_to_ring_min_dist"] = np.min(dist_matrix_ring1) if np.isfinite(np.min(dist_matrix_ring1)) else 0.0
    else:
        distances["lig_ring_to_ring_min_dist"] = 0.0

    # Calculate ring-to-ring distances between mol1 and mol2
    if ring_atoms1 and ring_atoms2:
        ring_coords1 = coords1[ring_atoms1]
        ring_coords2 = coords2[ring_atoms2]
        dist_matrix_ring2 = np.linalg.norm(
            ring_coords1[:, None, :] - ring_coords2[None, :, :], axis=2
        )
        distances["com_ring_to_ring_min_dist"] = np.min(dist_matrix_ring2) if np.isfinite(np.min(dist_matrix_ring2)) else 0.0
    else:
        distances["com_ring_to_ring_min_dist"] = 0.0

    return distances

def extract_metal_coordinates(pdb_file):
    """
    Extracts coordinates of non-bonded metal atoms from a PDB file.

    Parameters:
        pdb_file (str): Path to the PDB file.

    Returns:
        np.ndarray: Array of metal atom coordinates.
    """
    metal_coords = []
    valid_metals = ["MO", "NI", "HG", "CU", "CD", "ZN", "AG", "CO", "FE", "OS", "NA", "K", "MG", "MN", "SR", "CS", "BA", "CA", "PB"]
    with open(pdb_file, "r") as file:
        for line in file:
            if line.startswith("HETATM") and line[76:78].strip() in valid_metals:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                metal_coords.append([x, y, z])
    return np.array(metal_coords)

def calculate_metal_distances(ligand, receptor, metal_coords):
    """
    Calculates minimum distances from ligand and receptor atoms to metal atoms.

    Parameters:
        ligand (Mol): RDKit molecule for the ligand.
        receptor (Mol): RDKit molecule for the receptor.
        metal_coords (np.ndarray): Coordinates of metal atoms.

    Returns:
        dict: Minimum distances from ligand and receptor atoms to metal atoms.
    """
    ligand_conf = ligand.GetConformer()
    ligand_coords = np.array([list(ligand_conf.GetAtomPosition(i)) for i in range(ligand.GetNumAtoms())])

    receptor_conf = receptor.GetConformer()
    receptor_coords = np.array([list(receptor_conf.GetAtomPosition(i)) for i in range(receptor.GetNumAtoms())])

    distances = {}

    if metal_coords.size > 0:
        # Compute distances from ligand to metals
        ligand_metal_distances = np.linalg.norm(ligand_coords[:, None, :] - metal_coords[None, :, :], axis=2)
        distances["ligand_to_metal_min_dist"] = np.min(ligand_metal_distances)

        # Compute distances from receptor to metals
        receptor_metal_distances = np.linalg.norm(receptor_coords[:, None, :] - metal_coords[None, :, :], axis=2)
        distances["receptor_to_metal_min_dist"] = np.min(receptor_metal_distances)
    else:
        distances["ligand_to_metal_min_dist"] = 0.0
        distances["receptor_to_metal_min_dist"] = 0.0

    return distances

def process_directory(base_dir):
    """
    Process a single directory to calculate distances and write to distances.txt.

    Args:
        base_dir (str): Path to the input directory.
    """
    atom_pairs = [
        ('N', 'C'), ('N', 'O'), ('N', 'N'),
        ('C', 'C'), ('C', 'O'), ('O', 'O')
    ]

    # Define file paths in the given directory
    folder_name = os.path.basename(os.path.normpath(base_dir))

    output_file = os.path.join(base_dir, "distances.txt")
    # Dynamically find the SDF file ending with out_sorted.sdf
    ligand_file = None
    for fname in os.listdir(base_dir):
        if fname.endswith("out_sorted.sdf"):
            ligand_file = os.path.join(base_dir, fname)
            break

    if not ligand_file or not os.path.exists(ligand_file):
        print(f"No ligand SDF file ending with 'out_sorted.sdf' found in {base_dir}.")
        return

    receptor_file = os.path.join(base_dir, f"{folder_name}.mol2")
    pdb_file = os.path.join(base_dir, f"{folder_name}_h.pdb")

    # Check for existing output
    if os.path.exists(output_file):
        with open(output_file, "r") as f:
            content = f.read()
            if "Ligand 100 Scores:" in content:
                print(f"Skipping {folder_name}: distances.txt is already fully processed.")
                return

    start_time = time.time()
    print(f"Processing {folder_name}...")

    # Validate that required input files exist
    if not os.path.exists(ligand_file) or not os.path.exists(receptor_file):
        print(f"Missing ligand or receptor files in {base_dir}.")
        return

    # Load receptor
# Attempt receptor loading with and without sanitization
    receptor = Chem.MolFromMol2File(receptor_file)
    if receptor is None:
        print(f"[!] Default sanitization failed for receptor {folder_name}, trying without sanitization...")
        receptor = Chem.MolFromMol2File(receptor_file, sanitize=False, removeHs=False)

    if receptor is None:
        print(f"[ERROR] Receptor could not be loaded even after disabling sanitization. Skipping {folder_name}.")
        return

    # Extract metal coordinates
    metal_coords = extract_metal_coordinates(pdb_file) if os.path.exists(pdb_file) else np.array([])

    # Prepare output file
    with open(output_file, "w") as f_out:
        # Load ligands
        supplier = Chem.SDMolSupplier(ligand_file, sanitize=False, removeHs=False)


        for i, ligand in enumerate(supplier, start=1):  # Start count from 1
            if ligand is None:
                print(f"Failed to load ligand {i} in {folder_name}")
                continue

            # Compute intra-ligand distances
            ligand_distances = calculate_min_distances(ligand, atom_pairs=atom_pairs)

            # Compute ligand-receptor distances
            ligand_receptor_distances = calculate_min_distances(ligand, receptor, atom_pairs=atom_pairs)

            # Compute ring atom distances
            ring_distances = calculate_ring_atom_distances(ligand, receptor)

            # Compute metal distances
            metal_distances = calculate_metal_distances(ligand, receptor, metal_coords)

            # Write results
            f_out.write(f"Ligand {i} Scores:\n")
            for pair, dist in ligand_distances.items():
                f_out.write(f"    lig_{pair}_min_dist: {dist:.3f}\n")
            for pair, dist in ligand_receptor_distances.items():
                f_out.write(f"    com_{pair}_mindist: {dist:.3f}\n")
            for key, dist in ring_distances.items():
                f_out.write(f"    {key}: {dist:.3f}\n")
            for key, dist in metal_distances.items():
                f_out.write(f"    {key}: {dist:.3f}\n")
            f_out.write("\n")

    elapsed_time = time.time() - start_time
    print(f"{folder_name} complete. Time taken: {elapsed_time:.2f} seconds.")


# Function to combine descriptors after calculating distances
def combine_descriptors(descriptors_file, rdkit_file, output_file):
    """
    Combine descriptors from descriptors_full.txt and distances.txt.

    Args:
        descriptors_file (str): Path to descriptors_full.txt file.
        rdkit_file (str): Path to distances.txt file.
        output_file (str): Path to save the combined descriptors.

    Returns:
        None
    """
    combined_data = []

    # Read descriptors_full.txt
    with open(descriptors_file, 'r') as f:
        descriptors_data = f.readlines()

    # Read distances.txt
    with open(rdkit_file, 'r') as f:
        rdkit_data = f.readlines()

    # Combine data by ligand number
    ligand_number = None
    rdkit_dict = {}

    # Parse distances.txt into a dictionary
    for line in rdkit_data:
        line = line.strip()
        if line.startswith("Ligand ") and "Scores" in line:
            ligand_number = line
            rdkit_dict[ligand_number] = []
        elif ligand_number:
            rdkit_dict[ligand_number].append(f"  {line}")  # Indent lines

    # Combine descriptors_full.txt with distances.txt
    for line in descriptors_data:
        line = line.strip()
        if line.startswith("Ligand ") and "Scores" in line:
            ligand_number = line
            combined_data.append(line)
            if ligand_number in rdkit_dict:
                combined_data.extend(rdkit_dict[ligand_number])
        else:
            combined_data.append(f"  {line}" if line and not line.startswith("Ligand ") else line)

    # Ensure no extra blank lines
    cleaned_data = [line for line in combined_data if line.strip() != ""]

    # Save combined data
    with open(output_file, 'w') as f:
        f.write('\n'.join(cleaned_data))


def combine_after_distances(input_dir):
    """
    Combine distances.txt with descriptors_full.txt to produce descriptors_all.txt.

    Args:
        input_dir (str): Directory containing distances.txt and descriptors_full.txt.

    Returns:
        None
    """
    folder_name = os.path.basename(input_dir)
    descriptors_file = os.path.join(input_dir, "descriptors_full.txt")
    rdkit_file = os.path.join(input_dir, "distances.txt")
    output_file = os.path.join(input_dir, "descriptors_all.txt")

    # Check if both files exist
    if os.path.exists(descriptors_file) and os.path.exists(rdkit_file):
        combine_descriptors(descriptors_file, rdkit_file, output_file)
        print(f"Combined descriptors saved to {output_file}")
    else:
        print(f"Missing files for combining in {input_dir}. Ensure both descriptors_full.txt and distances.txt exist.")


# Entry point for walking through a parent directory
def wait_for_file(path):
    """Wait indefinitely until a file appears."""
    print(f"Waiting for {path} to appear...")
    while not os.path.exists(path):
        time.sleep(0.5)
    return True

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python rdkit_full.py <parent_directory>")
        sys.exit(1)

    parent_dir = sys.argv[1]
    if not os.path.isdir(parent_dir):
        print(f"Error: {parent_dir} is not a valid directory.")
        sys.exit(1)

    for subfolder in os.listdir(parent_dir):
        subfolder_path = os.path.join(parent_dir, subfolder)
        if os.path.isdir(subfolder_path):
            print(f"\n📂 Processing: {subfolder_path}")
            process_directory(subfolder_path)

            # Wait for distances.txt before combining
            distances_file = os.path.join(subfolder_path, "distances.txt")
            if wait_for_file(distances_file):
                combine_after_distances(subfolder_path)
