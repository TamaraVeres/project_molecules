import pandas as pd
from rdkit.Chem import Descriptors, MolToSmiles
import seaborn as sns 
import matplotlib.pyplot as plt
from mordred import Calculator
from mordred.Weight import Weight
from mordred.SLogP import SLogP
from mordred.TopoPSA import TopoPSA
from mordred.HydrogenBond import HBondDonor, HBondAcceptor
from mordred.RotatableBond import RotatableBondsCount

def calculate_rdkit_descriptors(mols):
    molecular_descriptors = []
    for mol in mols:
        molecular_descriptors.append({
            'Molecular_Weight': Descriptors.MolWt(mol),
            'H-Bond_Donor_Count': Descriptors.NumHDonors(mol),
            'H-Bond_Acceptor_Count': Descriptors.NumHAcceptors(mol),
            'Rotatable_Bond_Count': Descriptors.NumRotatableBonds(mol),
            'Polar_Area': Descriptors.TPSA(mol),
            'XLogP': Descriptors.MolLogP(mol),
        })

    df_molecular_descriptors = pd.DataFrame(molecular_descriptors)
    print("Descriptors summary:")
    print(df_molecular_descriptors.describe())
    print("First molecule descriptors:")
    print(df_molecular_descriptors.iloc[0])
    return df_molecular_descriptors

def distribution(df_molecular_descriptors):
    sns.histplot(df_molecular_descriptors['Molecular_Weight'])
    plt.title("Distribution of molecular weight")
    plt.show()
    sns.histplot(df_molecular_descriptors['Polar_Area'])
    plt.title("Distribution of tpsa")
    plt.show()
    sns.histplot(df_molecular_descriptors['XLogP'])
    plt.title("Distribution of XLogP")
    plt.show()
    sns.scatterplot(x='Polar_Area', y='XLogP', data=df_molecular_descriptors)
    plt.title("Polar Area vs XLogP")
    plt.show()

def lipinski_rule(mol):
    if mol is None:
        return False

    conditions = [
        Descriptors.MolWt(mol) <= 500,
        Descriptors.MolLogP(mol) <= 5,
        Descriptors.NumHDonors(mol) <= 5,
        Descriptors.NumHAcceptors(mol) <= 10,
        Descriptors.NumRotatableBonds(mol) <= 10,
        Descriptors.TPSA(mol) <= 140
    ]
    for condition in conditions:
        if not condition:
            return False
    print(f'The following compound: "{MolToSmiles(mol)}" passes Lipinski\'s rule')
    return True

def calculate_mordred_descriptors(mols):
    calculator = Calculator([
        Weight,
        HBondDonor,
        HBondAcceptor,
        RotatableBondsCount,
        TopoPSA,
        SLogP
    ], ignore_3D=True)
    mordred_descriptors = []
    for mol in mols:
        result = calculator(mol)
        mordred_descriptors.append(result.asdict())
        
    df_mordred_descriptors = pd.DataFrame(mordred_descriptors)
    print("Mordred column names:", df_mordred_descriptors.columns.tolist())
    print("Mordred descriptors summary:")
    print(df_mordred_descriptors.describe())
    print("First molecule descriptors:")
    print(df_mordred_descriptors.iloc[0])
    return df_mordred_descriptors