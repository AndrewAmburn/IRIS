#!/usr/bin/env python
# coding: utf-8

import os
import sys
import time
import shutil
import numpy as np
from math import sqrt
from collections import defaultdict
from rdkit import Chem
from rdkit.Chem import AllChem
from scipy.spatial import cKDTree

# =====================
#  CONSTANTS & MAPPINGS
# =====================

bases_fixed = ["A", "U", "C", "G"]
base_aliases = {
    "A": "A", "ADE": "A", "DA": "A", "RA": "A",
    "U": "U", "URA": "U", "DU": "U", "RU": "U",
    "C": "C", "CYT": "C", "DC": "C", "RC": "C",
    "G": "G", "GUA": "G", "DG": "G", "RG": "G"
}

monovalent_metals = ["NA", "K"]
divalent_metals = ["MG", "CA", "ZN", "MN"]

element_charge_map = {
    "O": -0.4,
    "N": -0.3,
    "P": 0.2,
    "S": 0.1,
    "C": 0.0
}

OUTPUT_DESCRIPTOR_NAME = "receptor_naffinity_descriptors.txt"

# ============================
#   PARSE MOL2 RECEPTOR FILE
# ============================

def parse_mol2_receptor(mol2_path):
    coords = []
    meta = []
    section = None

    with open(mol2_path, "r") as f:
        for line in f:
            s = line.rstrip()
            if s.startswith("@<TRIPOS>"):
                section = s
                continue

            if section == "@<TRIPOS>ATOM":
                p = s.split()
                if len(p) < 6:
                    continue
                try:
                    x, y, z = float(p[2]), float(p[3]), float(p[4])
                except Exception:
                    continue
                atom_name = p[1]
                atom_type = p[5]

                elem = atom_type.split(".")[0]
                elem = "".join([c for c in elem if not c.isdigit()]).upper()
                if elem == "":
                    elem = atom_name[0].upper()

                subst = p[7].upper() if len(p) >= 8 else ""

                coords.append((x, y, z))
                meta.append({
                    "name": atom_name,
                    "element": elem,
                    "resname": subst if subst else "UNK",
                    "resid": p[6] if len(p) >= 7 else None,
                    "is_hetatm": False
                })

    return np.array(coords), meta

# ============================
#  COMPUTE RECEPTOR FEATURES
# ============================

def compute_naffinity_descriptors(lig_mol, receptor_coords, receptor_meta):
    try:
        Chem.SanitizeMol(lig_mol)
    except Exception:
        try:
            for atom in lig_mol.GetAtoms():
                atom.SetIsAromatic(False)
            for bond in lig_mol.GetBonds():
                bond.SetIsAromatic(False)
            Chem.SanitizeMol(lig_mol)
        except Exception:
            Chem.Kekulize(lig_mol, clearAromaticFlags=True)

    AllChem.ComputeGasteigerCharges(lig_mol)
    conf = lig_mol.GetConformer()
    lig_coords = conf.GetPositions()

    ligand_charges = []
    for atom in lig_mol.GetAtoms():
        try:
            ligand_charges.append(float(atom.GetProp("_GasteigerCharge")))
        except Exception:
            ligand_charges.append(0.0)
    ligand_total_charge = sum(ligand_charges)

    lig_tree = cKDTree(lig_coords)
    rec_tree = cKDTree(receptor_coords) if receptor_coords.size else None

    close_atoms = set()
    if rec_tree is not None:
        hits = lig_tree.query_ball_tree(rec_tree, r=6.0)
        close_atoms = set(idx for group in hits for idx in group)

    features = defaultdict(int)

    distances = []
    mono_dists = []
    div_dists = []
    phos_dists = []
    base_res = set()

    dummy_charges = []
    weighted_q = []
    inv_d = []

    hacc = 0
    hdon = 0
    clashes = 0

    lig_centroid = np.mean(lig_coords, axis=0)

    for idx in close_atoms:
        meta = receptor_meta[idx]
        coord = receptor_coords[idx]

        d = float(np.min(np.linalg.norm(lig_coords - coord, axis=1)))
        distances.append(d)

        # Base residues
        if meta["resname"] in base_aliases:
            letter = base_aliases[meta["resname"]]
            features[f"Residue_{letter}"] += 1
            base_res.add(meta["resid"])

        if meta["resname"] in ("HOH", "WAT"):
            features["Residue_HOH"] += 1

        if meta["resname"] in monovalent_metals:
            features["MonovalentMetalCount"] += 1
            mono_dists.append(d)

        if meta["resname"] in divalent_metals:
            features["DivalentMetalCount"] += 1
            div_dists.append(d)

        if meta["name"].upper() in ("P", "OP1", "OP2", "OP3"):
            phos_dists.append(d)

        q = element_charge_map.get(meta["element"], 0.0)
        dummy_charges.append(q)
        weighted_q.append(q / (d + 1e-3))
        inv_d.append(1.0 / (d + 1e-3))

        if meta["element"] in ("O", "N"):
            hacc += 1
        if meta["name"].upper().startswith("H"):
            hdon += 1

        if d < 2.0:
            clashes += 1

    # Summary stats
    features["NumNearbyAtoms"] = len(close_atoms)
    features["MeanDistanceToNearbyAtoms"] = round(float(np.mean(distances)), 3) if distances else 0.0
    features["MinDistanceToReceptor"] = round(float(np.min(distances)), 3) if distances else 0.001
    features["MaxDistanceToReceptor"] = round(float(np.max(distances)), 3) if distances else 0.001
    features["DistanceToClosestPhosphate"] = round(min(phos_dists), 3) if phos_dists else 10.0
    features["DistanceToClosestMonovalent"] = round(min(mono_dists), 3) if mono_dists else 10.0
    features["DistanceToClosestDivalent"] = round(min(div_dists), 3) if div_dists else 10.0

    features["NumBaseResidues"] = len(base_res)

    if receptor_coords.size:
        rec_centroid = np.mean(receptor_coords, axis=0)
        features["LigandToReceptorCOM"] = round(float(np.linalg.norm(lig_centroid - rec_centroid)), 3)
    else:
        features["LigandToReceptorCOM"] = 0.0

    features["LigandCoordRange"] = round(float(np.ptp(lig_coords)), 3)

    total_rec_q = sum(dummy_charges)
    features["TotalReceptorPartialCharge"] = round(total_rec_q, 3)
    features["NetChargeBalance"] = round(ligand_total_charge - total_rec_q, 3)
    features["AbsChargeSum"] = round(sum(abs(x) for x in dummy_charges), 3)
    features["ChargeVariance"] = round(float(np.var(dummy_charges)), 5) if dummy_charges else 0.0
    features["ChargeDensity"] = round(float(np.mean(dummy_charges)), 5) if dummy_charges else 0.0
    features["DistWeightedChargeSum"] = round(sum(weighted_q), 3)
    features["InverseDistanceSum"] = round(sum(inv_d), 3)

    features["NearbyHBondAcceptors"] = hacc
    features["NearbyHBondDonors"] = hdon
    features["StericClashCount"] = clashes

    return features

# ========================
#  WRITE OUTPUT BLOCK
# ========================

def write_output_block(fout, pose_idx, features):
    fout.write(f"Ligand {pose_idx} Scores:\n")

    keys = [
        "NumNearbyAtoms","MeanDistanceToNearbyAtoms","MinDistanceToReceptor",
        "MaxDistanceToReceptor","DistanceToClosestPhosphate",
        "MonovalentMetalCount","DivalentMetalCount",
        "DistanceToClosestMonovalent","DistanceToClosestDivalent",
        "NumBaseResidues","LigandToReceptorCOM","LigandCoordRange",
        "TotalReceptorPartialCharge","NetChargeBalance","AbsChargeSum",
        "ChargeVariance","ChargeDensity","DistWeightedChargeSum",
        "InverseDistanceSum","NearbyHBondAcceptors","NearbyHBondDonors",
        "StericClashCount"
    ] + [f"Residue_{b}" for b in bases_fixed] + ["Residue_HOH"]

    for key in keys:
        fout.write(f"    {key}: {features.get(key, 0)}\n")

    fout.write("\n")

# ========================
#       MERGING LOGIC
# ========================

def load_blocks(path):
    blocks = {}
    order = []
    cur = None
    if not os.path.exists(path):
        return blocks, order

    with open(path, "r") as f:
        for line in f:
            s = line.rstrip("\n")
            if s.startswith("Ligand ") and "Scores" in s:
                cur = s
                if cur not in blocks:
                    blocks[cur] = []
                    order.append(cur)
            elif cur is not None:
                t = s.strip()
                if t != "":
                    blocks[cur].append(t)

    return blocks, order

def merge_receptor_into_all(descriptors_all_path, receptor_path):
    all_blocks, order = load_blocks(descriptors_all_path)
    rec_blocks, _ = load_blocks(receptor_path)

    if not all_blocks:
        print("No ligand blocks in descriptors_all.txt — skipping merge.")
        return

    combined = []
    for hdr in order:
        combined.append(hdr)

        # Original
        for line in all_blocks.get(hdr, []):
            combined.append(f"    {line}")

        # Receptor descriptors (no duplicates)
        if hdr in rec_blocks:
            for line in rec_blocks[hdr]:
                if line not in all_blocks.get(hdr, []):
                    combined.append(f"    {line}")

        combined.append("")

    while combined and combined[-1].strip() == "":
        combined.pop()

    with open(descriptors_all_path, "w") as f:
        f.write("\n".join(combined) + "\n")

# ========================
#          MAIN
# ========================

def main():
    if len(sys.argv) != 2:
        print("Usage: python receptor_naffinity_descriptors.py <directory_path>")
        sys.exit(1)

    directory = sys.argv[1]
    if not os.path.isdir(directory):
        print(f"Error: {directory} is not a directory.")
        sys.exit(1)

    folder_name = os.path.basename(directory)

    # Find ligand SDF (out_sorted)
    ligand_sdf = None
    for f in os.listdir(directory):
        if f.endswith("out_sorted.sd") or f.endswith("out_sorted.sdf"):
            ligand_sdf = os.path.join(directory, f)
            break

    if not ligand_sdf:
        print("No out_sorted.sd(f) file found.")
        sys.exit(1)

    mol2_file = os.path.join(directory, f"{folder_name}.mol2")
    if not os.path.exists(mol2_file):
        print(f"No receptor .mol2 found: {mol2_file}")
        sys.exit(1)

    # ===== Step 1: Create receptor_naffinity_descriptors.txt =====
    print(f"Generating receptor descriptors in {directory}")

    coords, meta = parse_mol2_receptor(mol2_file)
    supplier = Chem.SDMolSupplier(ligand_sdf, removeHs=False, sanitize=False)

    out_file = os.path.join(directory, OUTPUT_DESCRIPTOR_NAME)
    with open(out_file, "w") as fout:
        for i, lig in enumerate(supplier, start=1):
            if lig is None:
                print(f"Pose {i} is None — skipping.")
                continue
            feats = compute_naffinity_descriptors(lig, coords, meta)
            write_output_block(fout, i, feats)

    print("Done generating receptor descriptors.")

    # ===== Step 2: Merge with descriptors_all.txt =====
    descriptors_all = os.path.join(directory, "descriptors_all.txt")
    if not os.path.exists(descriptors_all):
        print("descriptors_all.txt not found — skipping merge.")
        return

    backup = os.path.join(directory, "descriptors_all_old.txt")
    try:
        shutil.copy2(descriptors_all, backup)
    except Exception as e:
        print(f"Backup failed: {e}")

    merge_receptor_into_all(descriptors_all, out_file)
    print("Merged receptor descriptors into descriptors_all.txt")

if __name__ == "__main__":
    main()
