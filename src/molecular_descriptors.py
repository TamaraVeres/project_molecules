import pandas as pd
import numpy as np
from rdkit.Chem import Descriptors, MolToSmiles
import seaborn as sns 
import matplotlib.pyplot as plt
from mordred import Calculator
from mordred.SLogP import SLogP
from rdkit.Chem import AllChem


def calculate_rdkit_descriptors(mols):
    molecular_descriptors = []
    for mol in mols:
        molecular_descriptors.append({
            'XLogP': Descriptors.MolLogP(mol),
        })

    df_molecular_descriptors = pd.DataFrame(molecular_descriptors)
    print("Descriptors summary:")
    print(df_molecular_descriptors.describe())
    print("First molecule descriptor:")
    print(df_molecular_descriptors.iloc[0])
    return df_molecular_descriptors

def distribution(df_molecular_descriptors):
    sns.histplot(df_molecular_descriptors['XLogP'])
    plt.title("Distribution of XLogP")
    plt.show()
  

def calculate_mordred_descriptors(mols):
    calculator = Calculator([
        SLogP
    ], ignore_3D=True)
    mordred_descriptors = []
    for mol in mols:
        result = calculator(mol)
        mordred_descriptors.append(result.asdict())
        
    df_mordred_descriptors = pd.DataFrame(mordred_descriptors)
    print("Mordred column names:", df_mordred_descriptors.columns.tolist())
    print("Mordred descriptor summary:")
    print(df_mordred_descriptors.describe())
    print("First molecule descriptor:")
    print(df_mordred_descriptors.iloc[0])
    return df_mordred_descriptors

def calculate_morgan_fingerprint(mols, radius=3, nbits=2048):
    fingerprints = []
    for mol in mols:
        fp = AllChem.GetMorganFingerprintAsBitVect(mol, radius, nbits)
        fingerprints.append(list(fp))
    return np.array(fingerprints)

