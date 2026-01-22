from src.utils import load_csv, basic_inspection, clean_molecules, draw_molecule
from src.molecular_descriptors import calculate_rdkit_descriptors, distribution, lipinski_rule
from src.compare import compare_descriptors


def molecules_main():
    df_compounds = load_csv(file_path="./pubchemdb.csv")
    basic_inspection(df_compounds)
    mols = clean_molecules(df_compounds)
    draw_molecule(mols[5])
    df_molecular_descriptors = calculate_rdkit_descriptors(mols)
    distribution(df_molecular_descriptors)
    for mol in mols:
        lipinski_rule(mol)
    comparison_table = compare_descriptors(df_molecular_descriptors, df_compounds)
    print(comparison_table.head(10))
    comparison_table.to_csv("./results/comparison_results.csv", index=False)
    comparison_table.to_html("./results/comparison_results.html", index=False)


if __name__ == "__main__":
    molecules_main()