import pandas as pd
from rdkit import Chem 
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Draw 


def load_data(file_path):
    df_compounds = pd.read_csv(file_path)
    print(df_compounds.head())
    return df_compounds

def basic_inspection(df_compounds):
    print(df_compounds.describe())
    print(df_compounds.head())
    return df_compounds


def cleaning_molecules(df_compounds):
    mols = []
    unique_smiles = set()

    for smile in df_compounds["smiles"]:
        mol = Chem.MolFromSmiles(smile)
        if mol is not None:
            mol = rdMolStandardize.Cleanup(mol)
            canonical_smile = Chem.MolToSmiles(mol)
            if canonical_smile not in unique_smiles:
                unique_smiles.add(canonical_smile)
                mols.append(mol)
    print("Number of initial molecules:", len(df_compounds["smiles"]))
    print("Number of unique valid molecules:", len(mols))
    print("Number of unique canonical smiles:", len(unique_smiles))
    return mols

def draw_molecule(mol):
    Draw.MolToFile(mol, "molecule.png")
    return mol

