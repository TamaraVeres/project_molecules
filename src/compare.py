import pandas as pd

def compare_descriptors(df_molecular_descriptors, df_compounds):
    comparison_df = pd.DataFrame()
    comparison_df["smiles"] = df_compounds["smiles"]
    comparison_df["molecular_weight_experimental"] = df_compounds["Molecular_Weight"]
    comparison_df["molecular_weight_rdkit"] = df_molecular_descriptors["Molecular_Weight"]
    comparison_df["molecular_weight_error"] = comparison_df["molecular_weight_experimental"] - comparison_df["molecular_weight_rdkit"]
    comparison_df["h_bond_donors_experimental"] = df_compounds["H-Bond_Donor_Count"]
    comparison_df["h_bond_donors_rdkit"] = df_molecular_descriptors["H-Bond_Donor_Count"]
    comparison_df["h_bond_donors_error"] = comparison_df["h_bond_donors_experimental"] - comparison_df["h_bond_donors_rdkit"]
    comparison_df["h_bond_acceptors_experimental"] = df_compounds["H-Bond_Acceptor_Count"]
    comparison_df["h_bond_acceptors_rdkit"] = df_molecular_descriptors["H-Bond_Acceptor_Count"]
    comparison_df["h_bond_acceptors_error"] = comparison_df["h_bond_acceptors_experimental"] - comparison_df["h_bond_acceptors_rdkit"]
    comparison_df["rotatable_bonds_experimental"] = df_compounds["Rotatable_Bond_Count"]
    comparison_df["rotatable_bonds_rdkit"] = df_molecular_descriptors["Rotatable_Bond_Count"]
    comparison_df["rotatable_bonds_error"] = comparison_df["rotatable_bonds_experimental"] - comparison_df["rotatable_bonds_rdkit"]
    comparison_df["tpsa_experimental"] = df_compounds["Polar_Area"]
    comparison_df["tpsa_rdkit"] = df_molecular_descriptors["Polar_Area"]
    comparison_df["tpsa_error"] = comparison_df["tpsa_experimental"] - comparison_df["tpsa_rdkit"]
    comparison_df["logp_experimental"] = df_compounds["XLogP"]
    comparison_df["logp_rdkit"] = df_molecular_descriptors["XLogP"]
    comparison_df["logp_error"] = comparison_df["logp_experimental"] - comparison_df["logp_rdkit"]
    return comparison_df
   
       
       