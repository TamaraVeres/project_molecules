from basic_inspection_compounds import load_data, basic_inspection, cleaning_molecules, draw_molecule
from molecular_descriptors import adding_descriptors, distribution, lipinski_rule

def molecules_main():
    df_compounds = load_data(file_path="./compounds_with_names.csv")
    basic_inspection(df_compounds)
    mols = cleaning_molecules(df_compounds)
    df_molecular_descriptors = adding_descriptors(mols)
    distribution(df_molecular_descriptors)
    lipinski_rule(mols[0])
    print(lipinski_rule(mols[0]))
    print(df_molecular_descriptors['logp'].head(10))
    draw_molecule(mols[5])

if __name__ == "__main__":
    molecules_main()