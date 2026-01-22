import pandas as pd
from rdkit import Chem 
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Draw 
from mordred import Calculator, descriptors


def load_csv(file_path):
    df_compounds = pd.read_csv(file_path)
    print(df_compounds.head())
    return df_compounds

def basic_inspection(df_compounds):
    print(df_compounds.describe())
    print(df_compounds.head())
    return df_compounds


def clean_molecules(df_compounds):
    mols = []
    unique_smiles = set()
    invalid_smiles = []

    for smile in df_compounds["smiles"]:
        mol = Chem.MolFromSmiles(smile)
        if mol is not None:
            # mol = rdMolStandardize.Cleanup(mol)
            canonical_smile = Chem.MolToSmiles(mol)
            if canonical_smile not in unique_smiles:
                unique_smiles.add(canonical_smile)
                mols.append(mol)
            else:
                invalid_smiles.append(smile)
    print("Number of initial molecules:", len(df_compounds["smiles"]))
    print("Number of unique valid molecules:", len(mols))
    print("Number of unique canonical smiles:", len(unique_smiles))
    if invalid_smiles:
        print("Invalid smiles:")
        for i, smile in enumerate(invalid_smiles, 1):
            print(f"{i}. {smile}")
    return mols
    

def draw_molecule(mol):
    Draw.MolToFile(mol, "molecule.png")
    return mol

