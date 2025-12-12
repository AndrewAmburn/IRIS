#!/usr/bin/env python3
# coding: utf-8

import os
import sys
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def calculate_3d_descriptors(ligand):
    try:
        return {
            "PMI1": rdMolDescriptors.CalcPMI1(ligand),
            "PMI2": rdMolDescriptors.CalcPMI2(ligand),
            "PMI3": rdMolDescriptors.CalcPMI3(ligand),
            "NPR1": rdMolDescriptors.CalcNPR1(ligand),
            "NPR2": rdMolDescriptors.CalcNPR2(ligand),
            "RadiusOfGyration": rdMolDescriptors.CalcRadiusOfGyration(ligand),
            "InertialShapeFactor": rdMolDescriptors.CalcInertialShapeFactor(ligand),
            "Eccentricity": rdMolDescriptors.CalcEccentricity(ligand),
            "Asphericity": rdMolDescriptors.CalcAsphericity(ligand),
            "SpherocityIndex": rdMolDescriptors.CalcSpherocityIndex(ligand),
            "PBF": rdMolDescriptors.CalcPBF(ligand),
        }
    except Exception as e:
        print(f"[WARN] 3D descriptor calc failed: {e}")
        return {}

def find_out_sorted_sdf(directory):
    for fname in os.listdir(directory):
        if fname.endswith("out_sorted.sdf"):
            return os.path.join(directory, fname)
    return None

def write_descriptors_to_file(descriptor_data, output_file):
    if not descriptor_data:
        return
    headers = list(descriptor_data[0].keys())
    with open(output_file, "w") as f:
        f.write("\t".join(headers) + "\n")
        for row in descriptor_data:
            f.write("\t".join(str(row.get(h, "")) for h in headers) + "\n")

def read_3d_descriptors(file_path):
    with open(file_path, "r") as f:
        headers = f.readline().strip().split("\t")
        out = []
        for line in f:
            vals = line.strip().split("\t")
            out.append(dict(zip(headers, vals)))
    return out

def append_3d_descriptors_to_descriptors_txt(descriptors_txt, descriptors_data):
    with open(descriptors_txt, "r") as f:
        lines = f.readlines()

    # If already appended once, do nothing
    keywords = [
        "PMI1:", "PMI2:", "PMI3:", "NPR1:", "NPR2:",
        "RadiusOfGyration:", "InertialShapeFactor:", "Eccentricity:",
        "Asphericity:", "SpherocityIndex:", "PBF:"
    ]
    if any(any(k in line for k in keywords) for line in lines):
        print(f"[SKIP] 3D descriptors already appended in: {descriptors_txt}")
        return

    ligand_idx = 0
    with open(descriptors_txt, "w") as f:
        for line in lines:
            f.write(line)
            if line.startswith("Ligand") and "Scores:" in line:
                if ligand_idx < len(descriptors_data):
                    d = descriptors_data[ligand_idx]
                    # preserve consistent indentation with the rest of your pipeline
                    for key, val in d.items():
                        f.write(f"  {key}: {val}\n")
                    ligand_idx += 1

    print(f"[OK] Appended 3D RDKit descriptors into: {descriptors_txt}")

def process_complex_dir(directory):
    sdf_path = find_out_sorted_sdf(directory)
    if not sdf_path:
        print(f"[SKIP] No *out_sorted.sdf found in: {directory}")
        return

    out_3d = os.path.join(directory, "rdkit_3D_descriptors.txt")
    descriptors_txt = os.path.join(directory, "descriptors.txt")

    # Generate rdkit_3D_descriptors.txt
    if os.path.exists(out_3d):
        print(f"[SKIP] {out_3d} exists")
    else:
        suppl = Chem.SDMolSupplier(sdf_path, sanitize=False, removeHs=False)
        rows = []
        for i, lig in enumerate(suppl, start=1):
            if lig is None:
                print(f"[WARN] Ligand {i} is None, skipping")
                rows.append({})
                continue
            try:
                lig_h = Chem.AddHs(lig, addCoords=True)
            except Exception:
                lig_h = lig
            rows.append(calculate_3d_descriptors(lig_h))

        # choose headers from first non-empty dict
        first_nonempty = next((r for r in rows if r), None)
        if not first_nonempty:
            print(f"[ERROR] No valid 3D descriptors computed in: {directory}")
            return

        # ensure every row has same keys
        keys = list(first_nonempty.keys())
        normalized = [{k: r.get(k, "") for k in keys} for r in rows]
        write_descriptors_to_file(normalized, out_3d)
        print(f"[OK] Wrote: {out_3d}")

    # Append into descriptors.txt
    if not os.path.exists(descriptors_txt):
        print(f"[SKIP] Missing {descriptors_txt}; cannot append 3D descriptors.")
        return

    if not os.path.exists(out_3d):
        print(f"[SKIP] Missing {out_3d}; cannot append 3D descriptors.")
        return

    data = read_3d_descriptors(out_3d)
    append_3d_descriptors_to_descriptors_txt(descriptors_txt, data)

def main():
    if len(sys.argv) != 2:
        print("Usage: python rdkit_3d_descriptors.py <complex_directory>")
        sys.exit(1)

    directory = sys.argv[1]
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a valid directory.")
        sys.exit(1)

    process_complex_dir(directory)

if __name__ == "__main__":
    main()
