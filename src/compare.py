import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

def compare_descriptors(df_molecular_descriptors, df_mordred_descriptors, df_compounds):
    comparison_df = pd.DataFrame()
    comparison_df["smiles"] = df_compounds["smiles"]
    comparison_df["molecular_weight_experimental"] = df_compounds["Molecular_Weight"]
    comparison_df["molecular_weight_rdkit"] = df_molecular_descriptors["Molecular_Weight"]
    comparison_df["molecular_weight_mordred"] = df_mordred_descriptors["MW"]
    comparison_df["molecular_weight_error"] = comparison_df["molecular_weight_experimental"] - comparison_df["molecular_weight_rdkit"]
    comparison_df["tpsa_experimental"] = df_compounds["Polar_Area"]
    comparison_df["tpsa_rdkit"] = df_molecular_descriptors["Polar_Area"]
    comparison_df["tpsa_mordred"] = df_mordred_descriptors["TopoPSA"]
    comparison_df["tpsa_error"] = comparison_df["tpsa_experimental"] - comparison_df["tpsa_rdkit"]
    comparison_df["logp_experimental"] = df_compounds["XLogP"]
    comparison_df["logp_rdkit"] = df_molecular_descriptors["XLogP"]
    comparison_df["logp_mordred"] = df_mordred_descriptors["SLogP"]
    comparison_df["logp_error"] = comparison_df["logp_experimental"] - comparison_df["logp_rdkit"]
    return comparison_df
   
    
def scatter_plot(df_compounds, df_molecular_descriptors):
    plt.scatter(df_compounds["XLogP"], df_molecular_descriptors["XLogP"])
    plt.xlabel("Experimental XLogP")
    plt.ylabel("RDKit XLogP")
    plt.title("Experimental vs RDKit XLogP")
    plt.savefig("./results/logp_results.png") 
    plt.show()