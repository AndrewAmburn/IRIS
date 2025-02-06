import os
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors, Descriptors, rdmolops
import time
import sys 

def load_and_fix_protein(mol2_path):
    protein = Chem.MolFromMol2File(mol2_path, sanitize=False, removeHs=False)
    if not protein:
        print(f"Could not load molecule from {mol2_path}")
        return None

    # Add hydrogens explicitly
    protein = rdmolops.AddHs(protein)
    return protein

def calculate_hydrophobic_contacts(ligand, protein):
    hydrophobic_contacts = 0
    hydrophobic_atoms = {6, 1}  # Carbon, Hydrogen are typically considered hydrophobic
    for l_atom in ligand.GetAtoms():
        if l_atom.GetAtomicNum() in hydrophobic_atoms:
            l_atom_pos = ligand.GetConformer().GetAtomPosition(l_atom.GetIdx())
            for p_atom in protein.GetAtoms():
                if p_atom.GetAtomicNum() in hydrophobic_atoms:
                    p_atom_pos = protein.GetConformer().GetAtomPosition(p_atom.GetIdx())
                    if (l_atom_pos - p_atom_pos).Length() <= 5.0:
                        hydrophobic_contacts += 1
    return hydrophobic_contacts

def calculate_closest_polar_contact(ligand, protein):
    polar_contacts = float('inf')  # Start with a large number
    polar_atoms = {7, 8, 16}  # Nitrogen, Oxygen, and Sulfur atom types
    for l_atom in ligand.GetAtoms():
        if l_atom.GetAtomicNum() in polar_atoms:
            l_atom_pos = ligand.GetConformer().GetAtomPosition(l_atom.GetIdx())
            for p_atom in protein.GetAtoms():
                p_atom_pos = protein.GetConformer().GetAtomPosition(p_atom.GetIdx())
                distance = (l_atom_pos - p_atom_pos).Length()
                if distance < polar_contacts:
                    polar_contacts = distance
    return polar_contacts if polar_contacts != float('inf') else 0  # Return 0 if no contact is made

def get_ligand_descriptors(ligand, protein):
    """ Calculate unique chemical descriptors for each ligand pose including interaction-based descriptors. """
    descriptors = {}
    descriptors['HBondDonors'] = rdMolDescriptors.CalcNumHBD(ligand)
    descriptors['HBondAcceptors'] = rdMolDescriptors.CalcNumHBA(ligand)
    descriptors['MolecularWeight'] = Descriptors.MolWt(ligand)
    descriptors['LogP'] = Descriptors.MolLogP(ligand)
    descriptors['TPSA'] = rdMolDescriptors.CalcTPSA(ligand)
    descriptors['HydrophobicContacts'] = calculate_hydrophobic_contacts(ligand, protein)
    descriptors['ClosestPolarContact'] = calculate_closest_polar_contact(ligand, protein)
    descriptors['RotatableBonds'] = rdMolDescriptors.CalcNumRotatableBonds(ligand)
    return descriptors

def write_descriptors_to_file(descriptor_data, output_file):
    with open(output_file, 'w') as f:
        headers = ['HBondDonors', 'HBondAcceptors', 'MolecularWeight', 'LogP', 'TPSA', 'HydrophobicContacts', 'ClosestPolarContact', 
                   'RotatableBonds']
        f.write('\t'.join(headers) + '\n')
        for data in descriptor_data:
            line = '\t'.join([str(data[h]) for h in headers])
            f.write(line + '\n')
    print(f"Descriptors written to {output_file}")

def process_directory(directory, max_retries=3, retry_delay=5):
    """
    Process all SDF files in the specified directory for descriptor calculations.
    """
    for root, dirs, files in os.walk(directory):
        for file in files:
            if file.endswith('out_sorted.sdf'):
                sdf_path = os.path.join(root, file)
                output_file = os.path.join(root, "lig_pcd.txt")
                mol2_files = [f for f in os.listdir(root) if f.endswith('.mol2')]

                if not mol2_files:
                    print(f"No Mol2 files found in {root}. Skipping...")
                    continue

                mol2_path = os.path.join(root, mol2_files[0])
                print(f"Mol2 file found: {mol2_path}")
                protein = load_and_fix_protein(mol2_path)
                if not protein:
                    print(f"Could not load protein from {mol2_path}. Skipping...")
                    continue

                print(f"Processing SDF file: {sdf_path}")

                # Retry logic for loading SDF file
                for attempt in range(max_retries):
                    try:
                        suppl = Chem.SDMolSupplier(sdf_path, removeHs=False)
                        if not suppl:
                            raise OSError(f"Could not load SDF file from {sdf_path}")

                        descriptors = []
                        valid_ligands = False
                        for ligand in suppl:
                            if ligand:
                                valid_ligands = True
                                try:
                                    ligand = rdmolops.AddHs(ligand)
                                    Chem.SanitizeMol(ligand)
                                    desc = get_ligand_descriptors(ligand, protein)
                                    if desc:
                                        descriptors.append(desc)
                                except Exception as e:
                                    print(f"Error processing ligand in {sdf_path}: {e}")
                            else:
                                print(f"No valid ligand found in {sdf_path}")

                        if descriptors:
                            print(f"Writing {len(descriptors)} descriptors to {output_file}")
                            write_descriptors_to_file(descriptors, output_file)
                        elif valid_ligands:
                            print(f"Ligands were found but no descriptors were generated. Skipping writing...")
                        else:
                            print(f"No valid ligands found in {sdf_path}. Skipping...")

                        break  # Break out of the retry loop if successful
                    except OSError as e:
                        print(f"Attempt {attempt + 1} failed for {sdf_path}: {e}")
                        if attempt < max_retries - 1:
                            time.sleep(retry_delay)
                        else:
                            print(f"Max retries reached for {sdf_path}. Skipping this file.")

# Directory containing SDF and Mol2 files
def main():
    if len(sys.argv) != 2:
        print("Usage: python rdkit_descriptors.py <directory_path>")
        sys.exit(1)

    # Get the directory path from the command-line arguments
    directory_path = sys.argv[1]
    
    process_directory(directory_path)

if __name__ == "__main__":
    main()
