# Molecular Descriptor Analysis with RDKit

The aim of this project is to become familiar with RDKit by computing basic molecular descriptors and comparing them with reference descriptor values from PubChem.

---

## Installation

```bash
pip install -r requirements.txt
```

## Usage

Run the main analysis:

```bash
python molecules_main.py
```

---

## Project Structure

```
exercises_deepchem/
│
├── molecules_main.py          # Main script to run the analysis
├── requirements.txt            # Project dependencies
├── pubchemdb.csv             
├── README.md                  
│
├── src/                       
│   ├── utils.py              
│   ├── molecular_descriptors.py  
│   └── compare.py           
│
└── results/                 
    ├── comparison_results.csv
    └── comparison_results.html
```

---

## Basic Inspection and Preparation

First, the CSV file is loaded and inspected to understand the dataset. An example database containing molecular descriptors is included in this repository (pubchem.csv). The SMILES strings are then converted into RDKit molecule objects. Invalid molecules are removed, and the remaining molecules are cleaned and standardized. Duplicate molecules are identified using canonical SMILES and removed. Finally, individual molecules can be visualized by saving their structure as an image.

## Computing Molecular Descriptors with RDKit

Molecular descriptors are calculated for each molecule using RDKit. These include molecular weight, hydrogen bond donors and acceptors, rotatable bonds, polar surface area (TPSA), and XLogP. Summary statistics are used to give an overview of the descriptor values, and simple plots are created to visualize their distributions. At the end, Lipinski's rule of five is applied to check whether each molecule meets basic drug-likeness criteria.

## Comparison of the RDKit Data with the Experimental Data

Finally, the molecular descriptors calculated with RDKit are compared with the reference values from the dataset. For each descriptor, the experimental value, the RDKit value, and the error between them are calculated and stored in a table. The comparison results are saved as output files in the results folder.

---

## Conclusion

Overall, RDKit performs well when calculating basic molecular descriptors.