import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.metrics import r2_score

def compare_descriptors(df_molecular_descriptors, df_mordred_descriptors, df_compounds):
    comparison_df = pd.DataFrame()
    comparison_df["smiles"] = df_compounds["smiles"]
    comparison_df["logp_experimental"] = df_compounds["XLogP"]
    comparison_df["logp_rdkit"] = df_molecular_descriptors["XLogP"]
    comparison_df["logp_mordred"] = df_mordred_descriptors["SLogP"]
    comparison_df["logp_error"] = comparison_df["logp_experimental"] - comparison_df["logp_rdkit"]
    return comparison_df
   
    
def scatter_plot(df_compounds, df_molecular_descriptors):
    plt.scatter(df_compounds["XLogP"], df_molecular_descriptors["XLogP"])
    m, b = np.polyfit(df_compounds["XLogP"], df_molecular_descriptors["XLogP"], 1)
    plt.plot(df_compounds["XLogP"], m*df_compounds["XLogP"] + b, color='red')
    r2 = r2_score(df_compounds["XLogP"], df_molecular_descriptors["XLogP"])
    plt.text(1, -4, f"R² = {r2:.3f}", fontsize=12)
    plt.xlabel("Experimental XLogP")
    plt.ylabel("RDKit XLogP")
    plt.title("Experimental vs RDKit XLogP")
    plt.savefig("./results/logp_results.png") 
    plt.show()