import pandas as pd
from rdkit import Chem 
from rdkit.Chem.MolStandardize import rdMolStandardize
from rdkit.Chem import Draw 


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
    invalid_smiles = []

    for smile in df_compounds["smiles"]:
        mol = Chem.MolFromSmiles(smile)
        if mol is None:
            invalid_smiles.append(smile)
            continue
        mols.append(mol)
    print("Number of initial molecules:", len(df_compounds["smiles"]))
    print("Number of valid molecules:", len(mols))
    print("Number of invalid molecules:", len(invalid_smiles))
    return mols
    

def draw_molecule(mol):
    Draw.MolToFile(mol, "molecule.png")
    return mol

