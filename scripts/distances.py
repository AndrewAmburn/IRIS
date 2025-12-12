#!/usr/bin/env python
# coding: utf-8

import os
import sys
import time
from rdkit import Chem
from rdkit.Chem import AllChem  # kept in case you extend later
import numpy as np


def calculate_min_distances(mol1, mol2=None, atom_pairs=[("N", "C")]):
    """
    Calculates the minimum distances between specified atom pairs.

    Parameters:
        mol1 (Mol): The first molecule (e.g., ligand).
        mol2 (Mol): The second molecule (e.g., receptor), or None if only within mol1.
        atom_pairs (list of tuples): List of atom pairs, e.g., [("N", "C")].

    Returns:
        dict: Minimum distances for each atom pair.
    """
    if mol2 is None:
        mol2 = mol1

    conf1 = mol1.GetConformer()
    coords1 = np.array(
        [list(conf1.GetAtomPosition(i)) for i in range(mol1.GetNumAtoms())]
    )

    conf2 = mol2.GetConformer()
    coords2 = np.array(
        [list(conf2.GetAtomPosition(i)) for i in range(mol2.GetNumAtoms())]
    )

    atom_types1 = [atom.GetSymbol() for atom in mol1.GetAtoms()]
    atom_types2 = [atom.GetSymbol() for atom in mol2.GetAtoms()]

    distances = {}

    for atom_type1, atom_type2 in atom_pairs:
        indices1 = [i for i, atype in enumerate(atom_types1) if atype == atom_type1]
        indices2 = [i for i, atype in enumerate(atom_types2) if atype == atom_type2]

        if not indices1 or not indices2:
            distances[f"{atom_type1}-{atom_type2}"] = 0.0
            continue

        selected_coords1 = coords1[indices1]
        selected_coords2 = coords2[indices2]

        dist_matrix = np.linalg.norm(
            selected_coords1[:, None, :] - selected_coords2[None, :, :],
            axis=2,
        )

        # Exclude self-pairs and covalent bonds if mol1 == mol2
        if mol1 == mol2:
            for i, idx1 in enumerate(indices1):
                for j, idx2 in enumerate(indices2):
                    if idx1 == idx2 or mol1.GetBondBetweenAtoms(idx1, idx2):
                        dist_matrix[i, j] = np.inf

        min_dist = np.min(dist_matrix)
        distances[f"{atom_type1}-{atom_type2}"] = (
            float(min_dist) if np.isfinite(min_dist) else 0.0
        )

    return distances


def calculate_ring_atom_distances(mol1, mol2=None):
    """
    Calculates ring–other and ring–ring minimum distances within and between molecules.
    """
    if mol2 is None:
        mol2 = mol1

    conf1 = mol1.GetConformer()
    coords1 = np.array(
        [list(conf1.GetAtomPosition(i)) for i in range(mol1.GetNumAtoms())]
    )

    conf2 = mol2.GetConformer()
    coords2 = np.array(
        [list(conf2.GetAtomPosition(i)) for i in range(mol2.GetNumAtoms())]
    )

    ring_atoms1 = [atom.GetIdx() for atom in mol1.GetAtoms() if atom.IsInRing()]
    other_atoms1 = [atom.GetIdx() for atom in mol1.GetAtoms() if not atom.IsInRing()]

    ring_atoms2 = [atom.GetIdx() for atom in mol2.GetAtoms() if atom.IsInRing()]
    other_atoms2 = [atom.GetIdx() for atom in mol2.GetAtoms() if not atom.IsInRing()]

    distances = {}

    # ligand ring to ligand non-ring
    if ring_atoms1 and other_atoms1:
        ring_coords1 = coords1[ring_atoms1]
        other_coords1 = coords1[other_atoms1]
        dist_matrix1 = np.linalg.norm(
            ring_coords1[:, None, :] - other_coords1[None, :, :],
            axis=2,
        )
        distances["lig_ring_min_dist"] = float(np.min(dist_matrix1))
    else:
        distances["lig_ring_min_dist"] = 0.0

    # ligand ring to receptor non-ring
    if ring_atoms1 and other_atoms2:
        ring_coords1 = coords1[ring_atoms1]
        other_coords2 = coords2[other_atoms2]
        dist_matrix2 = np.linalg.norm(
            ring_coords1[:, None, :] - other_coords2[None, :, :],
            axis=2,
        )
        distances["com_ring_min_dist"] = float(np.min(dist_matrix2))
    else:
        distances["com_ring_min_dist"] = 0.0

    # ligand ring to ligand ring
    if ring_atoms1:
        ring_coords1 = coords1[ring_atoms1]
        dist_matrix_ring1 = np.linalg.norm(
            ring_coords1[:, None, :] - ring_coords1[None, :, :],
            axis=2,
        )
        np.fill_diagonal(dist_matrix_ring1, np.inf)
        min_val = np.min(dist_matrix_ring1)
        distances["lig_ring_to_ring_min_dist"] = (
            float(min_val) if np.isfinite(min_val) else 0.0
        )
    else:
        distances["lig_ring_to_ring_min_dist"] = 0.0

    # ligand ring to receptor ring
    if ring_atoms1 and ring_atoms2:
        ring_coords1 = coords1[ring_atoms1]
        ring_coords2 = coords2[ring_atoms2]
        dist_matrix_ring2 = np.linalg.norm(
            ring_coords1[:, None, :] - ring_coords2[None, :, :],
            axis=2,
        )
        min_val = np.min(dist_matrix_ring2)
        distances["com_ring_to_ring_min_dist"] = (
            float(min_val) if np.isfinite(min_val) else 0.0
        )
    else:
        distances["com_ring_to_ring_min_dist"] = 0.0

    return distances


def extract_metal_coordinates(pdb_file):
    """
    Extract coordinates of metal atoms from a PDB file.
    """
    metal_coords = []
    valid_metals = [
        "MO",
        "NI",
        "HG",
        "CU",
        "CD",
        "ZN",
        "AG",
        "CO",
        "FE",
        "OS",
        "NA",
        "K",
        "MG",
        "MN",
        "SR",
        "CS",
        "BA",
        "CA",
        "PB",
    ]
    if not os.path.exists(pdb_file):
        return np.array([])

    with open(pdb_file, "r") as fh:
        for line in fh:
            if line.startswith("HETATM") and line[76:78].strip() in valid_metals:
                x = float(line[30:38].strip())
                y = float(line[38:46].strip())
                z = float(line[46:54].strip())
                metal_coords.append([x, y, z])
    return np.array(metal_coords)


def calculate_metal_distances(ligand, receptor, metal_coords):
    """
    Calculate minimum ligand–metal and receptor–metal distances.
    """
    ligand_conf = ligand.GetConformer()
    ligand_coords = np.array(
        [list(ligand_conf.GetAtomPosition(i)) for i in range(ligand.GetNumAtoms())]
    )

    receptor_conf = receptor.GetConformer()
    receptor_coords = np.array(
        [list(receptor_conf.GetAtomPosition(i)) for i in range(receptor.GetNumAtoms())]
    )

    distances = {}

    if metal_coords.size > 0:
        ligand_metal_distances = np.linalg.norm(
            ligand_coords[:, None, :] - metal_coords[None, :, :],
            axis=2,
        )
        receptor_metal_distances = np.linalg.norm(
            receptor_coords[:, None, :] - metal_coords[None, :, :],
            axis=2,
        )
        distances["ligand_to_metal_min_dist"] = float(
            np.min(ligand_metal_distances)
        )
        distances["receptor_to_metal_min_dist"] = float(
            np.min(receptor_metal_distances)
        )
    else:
        distances["ligand_to_metal_min_dist"] = 0.0
        distances["receptor_to_metal_min_dist"] = 0.0

    return distances


def process_directory(base_dir):
    """
    Process a single complex directory to calculate distances and write distances.txt.
    """
    atom_pairs = [
        ("N", "C"),
        ("N", "O"),
        ("N", "N"),
        ("C", "C"),
        ("C", "O"),
        ("O", "O"),
    ]

    folder_name = os.path.basename(os.path.normpath(base_dir))
    output_file = os.path.join(base_dir, "distances.txt")

    # Find ligand SDF (out_sorted.sdf)
    ligand_file = None
    for fname in os.listdir(base_dir):
        if fname.endswith("out_sorted.sdf"):
            ligand_file = os.path.join(base_dir, fname)
            break

    if not ligand_file or not os.path.exists(ligand_file):
        print(f"[SKIP] No ligand SDF ending with 'out_sorted.sdf' in {base_dir}.")
        return

    receptor_file = os.path.join(base_dir, f"{folder_name}.mol2")
    pdb_file = os.path.join(base_dir, f"{folder_name}_h.pdb")

    if not os.path.exists(receptor_file):
        print(f"[SKIP] Receptor file {receptor_file} not found. Skipping {folder_name}.")
        return

    # If distances.txt already looks complete, skip
    if os.path.exists(output_file):
        with open(output_file, "r") as fh:
            content = fh.read()
        if "Ligand 100 Scores:" in content:
            print(
                f"[SKIP] {folder_name}: distances.txt appears fully processed. "
                "Delete it to recompute."
            )
            return

    print(f"Processing {folder_name} for distance features...")
    start_time = time.time()

    # Load receptor
    receptor = Chem.MolFromMol2File(receptor_file)
    if receptor is None:
        print(
            f"[WARN] Default sanitization failed for receptor {folder_name}, "
            "retrying without sanitization..."
        )
        receptor = Chem.MolFromMol2File(receptor_file, sanitize=False, removeHs=False)

    if receptor is None:
        print(
            f"[ERROR] Receptor could not be loaded even after disabling sanitization. "
            f"Skipping {folder_name}."
        )
        return

    metal_coords = extract_metal_coordinates(pdb_file)

    supplier = Chem.SDMolSupplier(ligand_file, sanitize=False, removeHs=False)

    with open(output_file, "w") as f_out:
        for i, ligand in enumerate(supplier, start=1):
            if ligand is None:
                print(f"[WARN] Failed to load ligand {i} in {folder_name}")
                continue

            ligand_distances = calculate_min_distances(
                ligand, atom_pairs=atom_pairs
            )
            ligand_receptor_distances = calculate_min_distances(
                ligand, receptor, atom_pairs=atom_pairs
            )
            ring_distances = calculate_ring_atom_distances(ligand, receptor)
            metal_distances = calculate_metal_distances(
                ligand, receptor, metal_coords
            )

            f_out.write(f"Ligand {i} Scores:\n")
            for pair, dist in ligand_distances.items():
                f_out.write(f"  lig_{pair}_min_dist: {dist:.3f}\n")
            for pair, dist in ligand_receptor_distances.items():
                f_out.write(f"  com_{pair}_mindist: {dist:.3f}\n")
            for key, dist in ring_distances.items():
                f_out.write(f"  {key}: {dist:.3f}\n")
            for key, dist in metal_distances.items():
                f_out.write(f"  {key}: {dist:.3f}\n")
            f_out.write("\n")

    elapsed_time = time.time() - start_time
    print(f"{folder_name} distance calculation complete in {elapsed_time:.2f} s.")
    print(f"[OK] Wrote distances to {output_file}")


def combine_descriptors(descriptors_file, distances_file, output_file):
    """
    Combine descriptors_full.txt with distances.txt into descriptors_all.txt.
    """
    combined_data = []

    with open(descriptors_file, "r") as f:
        descriptors_data = f.readlines()

    with open(distances_file, "r") as f:
        dist_data = f.readlines()

    ligand_number = None
    dist_dict = {}

    # Build dict from distances.txt
    for line in dist_data:
        line = line.strip()
        if line.startswith("Ligand ") and "Scores" in line:
            ligand_number = line
            dist_dict[ligand_number] = []
        elif ligand_number:
            dist_dict[ligand_number].append(f"  {line}")

    # Merge into descriptors_full.txt
    for line in descriptors_data:
        stripped = line.strip()
        if stripped.startswith("Ligand ") and "Scores" in stripped:
            ligand_number = stripped
            combined_data.append(stripped)
            if ligand_number in dist_dict:
                combined_data.extend(dist_dict[ligand_number])
        else:
            if stripped and not stripped.startswith("Ligand "):
                combined_data.append(f"  {stripped}")
            else:
                combined_data.append(stripped)

    cleaned_data = [ln for ln in combined_data if ln.strip() != ""]

    with open(output_file, "w") as f:
        f.write("\n".join(cleaned_data))

    print(f"[OK] Combined descriptors saved to {output_file}")


def combine_after_distances(input_dir):
    """
    Wrapper to combine distances.txt with descriptors_full.txt.
    """
    descriptors_file = os.path.join(input_dir, "descriptors_full.txt")
    distances_file = os.path.join(input_dir, "distances.txt")
    output_file = os.path.join(input_dir, "descriptors_all.txt")

    if os.path.exists(descriptors_file) and os.path.exists(distances_file):
        combine_descriptors(descriptors_file, distances_file, output_file)
    else:
        print(
            f"[WARN] Skipping combination in {input_dir}: "
            "descriptors_full.txt or distances.txt missing."
        )


def wait_for_file(path):
    """Wait until a file appears."""
    while not os.path.exists(path):
        time.sleep(0.5)
    return True


def main():
    if len(sys.argv) != 2:
        print("Usage: python distances.py <directory_path>")
        sys.exit(1)

    directory = sys.argv[1]
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        sys.exit(1)

    print(f"\nProcessing directory: {directory}")
    process_directory(directory)

    distances_file = os.path.join(directory, "distances.txt")
    if wait_for_file(distances_file):
        combine_after_distances(directory)


if __name__ == "__main__":
    main()
