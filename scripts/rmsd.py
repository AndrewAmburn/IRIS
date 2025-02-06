#!/usr/bin/env python
# coding: utf-8

import os
from spyrmsd import io, rmsd
import sys
import numpy as np

def calculate_rmsd(folder_path):
    folder_name = os.path.basename(folder_path)

    pdb_path = os.path.join(folder_path, f"{folder_name}_lig.pdb")
    sdf_path = os.path.join(folder_path, f"{folder_name}_docking_out_sorted.sdf")
    output_file_path = os.path.join(folder_path, f"{folder_name}_rmsd_output.txt")

    if os.path.exists(output_file_path):
        with open(output_file_path, 'r') as output_file:
            if "RMSD for molecule 99:" in output_file.read():
                print(f"RMSD already calculated for all molecules in {folder_name}. Skipping.")
                return

    try:
        ref = io.loadmol(pdb_path)
        mols = io.loadallmols(sdf_path)

        ref.strip()
        coords_ref = ref.coordinates
        anum_ref = ref.atomicnums

        # Debugging Reference Molecule
        #print(f"[DEBUG] Reference Molecule: {len(coords_ref)} atoms, Atomic numbers: {len(anum_ref)}")

        successful_rmsd_count = 0
        rmsd_results = []

        for i, mol in enumerate(mols[:100]):
            try:
                mol.strip()
                coords = mol.coordinates
                anum = mol.atomicnums

                # Debugging Molecule Details
                #print(f"[DEBUG] Molecule {i}: {len(coords)} atoms, Atomic numbers: {len(anum)}")

                # Validate Coordinate Shapes and Atomic Numbers
                if coords_ref.shape != coords.shape:
                    print(f"[SKIP] Molecule {i}: Coordinate shape mismatch ({coords_ref.shape} != {coords.shape})")
                    continue

                if len(anum_ref) != len(anum) or not np.array_equal(anum_ref, anum):
                    print(f"[SKIP] Molecule {i}: Atomic number mismatch")
                    continue

                # Calculate RMSD (Ignore Adjacency Matrix for Now)
                RMSD = rmsd.rmsd(coords_ref, coords, anum_ref, anum)
                rmsd_results.append(f"RMSD for molecule {i}: {RMSD:.3f}\n")
                successful_rmsd_count += 1

            except Exception as e:
                print(f"[ERROR] Molecule {i}: {e}")

        if successful_rmsd_count > 0:
            with open(output_file_path, 'w') as output_file:
                output_file.writelines(rmsd_results)
            print(f"RMSD values for {folder_name} have been written to {output_file_path}")
        else:
            print(f"No successful RMSD calculations for {folder_name}.")
    except Exception as e:
        print(f"[ERROR] Folder {folder_name}: {e}")

def main():
    if len(sys.argv) != 2:
        print("Usage: python rmsd.py <directory_path>")
        sys.exit(1)

    folder_path = sys.argv[1]

    if os.path.isdir(folder_path):
        calculate_rmsd(folder_path)
    else:
        print(f"Error: {folder_path} is not a valid directory.")
        sys.exit(1)

if __name__ == "__main__":
    main()
