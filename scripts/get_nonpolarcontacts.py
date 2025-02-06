import os
from rdkit import Chem
import sys 

# Calculate hydrophobic contacts between ligand and protein
def calculate_ligandprotein_hydrophobic_contacts(ligand, protein):
    hydrophobic_contacts = 0
    hydrophobic_atoms = {6, 1, 15, 16}  # Carbon, Hydrogen, Phosphorus, Sulfur are typically considered hydrophobic
    for l_atom in ligand.GetAtoms():
        if l_atom.GetAtomicNum() in hydrophobic_atoms:
            l_atom_pos = ligand.GetConformer().GetAtomPosition(l_atom.GetIdx())
            for p_atom in protein.GetAtoms():
                if p_atom.GetAtomicNum() in hydrophobic_atoms:
                    p_atom_pos = protein.GetConformer().GetAtomPosition(p_atom.GetIdx())
                    if (l_atom_pos - p_atom_pos).Length() <= 3.5:
                        hydrophobic_contacts += 1
    return hydrophobic_contacts

# Calculate hydrophobic contacts within the ligand
def calculate_ligand_hydrophobic_contacts(ligand):
    hydrophobic_contacts = 0
    hydrophobic_atoms = {6, 1, 15, 16}  # Carbon, Hydrogen, Phosphorus, Sulfur are typically considered hydrophobic
    ligand_conformer = ligand.GetConformer()

    for i, l_atom in enumerate(ligand.GetAtoms()):
        if l_atom.GetAtomicNum() in hydrophobic_atoms:
            l_atom_pos = ligand_conformer.GetAtomPosition(l_atom.GetIdx())
            for j in range(i + 1, len(ligand.GetAtoms())):
                l2_atom = ligand.GetAtomWithIdx(j)
                if l2_atom.GetAtomicNum() in hydrophobic_atoms:
                    l2_atom_pos = ligand_conformer.GetAtomPosition(l2_atom.GetIdx())
                    if (l_atom_pos - l2_atom_pos).Length() <= 3.5:
                        hydrophobic_contacts += 1
    return hydrophobic_contacts

# Load and fix protein from .mol2 file
def load_and_fix_protein(mol2_path):
    protein = Chem.MolFromMol2File(mol2_path, sanitize=False, removeHs=False)
    if not protein:
        print(f"Could not load molecule from {mol2_path}")
        return None
    protein = Chem.AddHs(protein)  # Add hydrogens
    return protein

# Process the folder to calculate hydrophobic contacts
def process_directory(directory):
    """
    Process all SDF and MOL2 files in the directory for hydrophobic contact calculations.
    """
    for root, dirs, files in os.walk(directory):
        output_file_path = os.path.join(root, 'nonpolars.txt')

        # Check if the file already exists
        if os.path.exists(output_file_path):
            print(f"{output_file_path} already exists. Skipping...")
            continue

        mol2_files = [f for f in files if f.endswith('.mol2')]
        if mol2_files:
            mol2_path = os.path.join(root, mol2_files[0])
            protein = load_and_fix_protein(mol2_path)
            if not protein:
                continue

            sdf_files = [f for f in files if f.endswith('out_sorted.sdf')]
            for sdf_file in sdf_files:
                sdf_path = os.path.join(root, sdf_file)
                suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
                descriptors = []
                for ligand in suppl:
                    if ligand:
                        ligand = Chem.AddHs(ligand)
                        Chem.SanitizeMol(ligand)
                        ligandprotein_hydrophobic_contacts = calculate_ligandprotein_hydrophobic_contacts(ligand, protein)
                        ligand_hydrophobic_contacts = calculate_ligand_hydrophobic_contacts(ligand)
                        descriptors.append({
                            'LigandProteinHydrophobicContacts': ligandprotein_hydrophobic_contacts,
                            'LigandHydrophobicContacts': ligand_hydrophobic_contacts
                        })

                if descriptors:
                    write_descriptors_to_file(descriptors, output_file_path)

# Write the hydrophobic contact descriptors to a file
def write_descriptors_to_file(descriptor_data, output_file):
    with open(output_file, 'w') as f:
        headers = ['LigandProteinHydrophobicContacts', 'LigandHydrophobicContacts']
        f.write('\t'.join(headers) + '\n')
        for data in descriptor_data:
            line = '\t'.join([str(data[h]) for h in headers])
            f.write(line + '\n')
    print(f"Descriptors successfully written to {output_file}")

# Extract nonpolar contact data from nonpolars.txt
def extract_nonpolar_contacts(nonpolars_file):
    nonpolar_data = []
    with open(nonpolars_file, 'r') as f:
        headers = f.readline().strip().split('\t')
        for line in f:
            values = line.strip().split('\t')
            nonpolar_data.append(dict(zip(headers, values)))
    return nonpolar_data

# Append the nonpolar contact data to the descriptors.txt file
def append_to_descriptors(descriptors_file, nonpolar_data):
    with open(descriptors_file, 'r') as f:
        lines = f.readlines()

    ligand_index = 0
    with open(descriptors_file, 'w') as f:
        for line in lines:
            f.write(line)
            if line.startswith("Ligand") and "Scores:" in line:
                ligand_index += 1
                if ligand_index <= len(nonpolar_data):
                    np_data = nonpolar_data[ligand_index - 1]
                    f.write(f"  LigandProteinHydrophobicContacts: {np_data['LigandProteinHydrophobicContacts']}\n")
                    f.write(f"  LigandHydrophobicContacts: {np_data['LigandHydrophobicContacts']}\n")

# Process the directory to append nonpolar contacts to descriptors
def process_nonpolar_append(directory):
    """
    Appends nonpolar contact data from nonpolars.txt to descriptors.txt.
    """
    for root, dirs, files in os.walk(directory):
        nonpolars_file_path = os.path.join(root, 'nonpolars.txt')
        descriptors_file_path = os.path.join(root, 'descriptors.txt')

        if os.path.exists(nonpolars_file_path) and os.path.exists(descriptors_file_path):
            nonpolar_data = extract_nonpolar_contacts(nonpolars_file_path)
            append_to_descriptors(descriptors_file_path, nonpolar_data)
        else:
            print(f"Files missing in directory {root}: Skipping...")

def main():
    # Ensure the correct number of arguments are passed
    if len(sys.argv) != 2:
        print("Usage: python get_nonpolarcontacts.py <directory_path>")
        sys.exit(1)

    # Get the folder path from command-line arguments
    folder_path = sys.argv[1]

    # Step 1: Calculate hydrophobic contacts and write them to nonpolars.txt
    process_directory(folder_path)
    print("Hydrophobic contact descriptors extraction and file writing complete.")

    # Step 2: Append the nonpolar contact descriptors to descriptors.txt
    process_nonpolar_append(folder_path)
    print("Appended nonpolar contact descriptors to descriptors.txt.")

if __name__ == '__main__':
    main()
