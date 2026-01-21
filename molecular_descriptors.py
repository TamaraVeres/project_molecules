import pandas as pd
from rdkit.Chem import Descriptors 
import seaborn as sns 
import matplotlib.pyplot as plt

def adding_descriptors(mols):
    molecular_descriptors = []
    for mol in mols:
        molecular_descriptors.append({
            'molecular_weight': Descriptors.MolWt(mol),
            'exact_molecular_weight': Descriptors.ExactMolWt(mol),
            'h_bond_donors': Descriptors.NumHDonors(mol),
            'h_bond_acceptors': Descriptors.NumHAcceptors(mol),
            'rotatable_bonds': Descriptors.NumRotatableBonds(mol),
            'aromatic_rings': Descriptors.NumAromaticRings(mol),
            'saturated_rings': Descriptors.NumSaturatedRings(mol),
            'aromatic_carbocycles': Descriptors.NumAromaticCarbocycles(mol),
            'aromatic_heterocycles': Descriptors.NumAromaticHeterocycles(mol),
            'tpsa': Descriptors.TPSA(mol),
            'logp': Descriptors.MolLogP(mol),
        })
    
    df_molecular_descriptors = pd.DataFrame(molecular_descriptors)
    print("Descriptors summary:")
    print(df_molecular_descriptors.describe())
    print("First molecule descriptors:")
    print(df_molecular_descriptors.iloc[0])
    return df_molecular_descriptors

def distribution(df_molecular_descriptors):
    sns.histplot(df_molecular_descriptors['molecular_weight'])
    plt.title("Distribution of molecular weight")
    plt.show()
    sns.histplot(df_molecular_descriptors['exact_molecular_weight'])
    plt.title("Distribution of exact molecular weight")
    plt.show()
    sns.histplot(df_molecular_descriptors['tpsa'])
    plt.title("Distribution of tpsa")
    plt.show()
    sns.histplot(df_molecular_descriptors['logp'])
    plt.title("Distribution of logp")
    plt.show()
    sns.scatterplot(x='tpsa', y='logp', data=df_molecular_descriptors)
    plt.title("TPSA vs logP")
    plt.show()

def lipinski_rule(mol):
    return {
        "MW <= 500": Descriptors.MolWt(mol) <= 500,
        "logP <= 5": Descriptors.MolLogP(mol) <= 5,
        "HBD <= 5": Descriptors.NumHDonors(mol) <= 5,
        "HBA <= 10": Descriptors.NumHAcceptors(mol) <= 10,
        "Rotatable bonds <= 10": Descriptors.NumRotatableBonds(mol) <= 10,
        "PSA <= 140": Descriptors.TPSA(mol) <= 140   
    }
 