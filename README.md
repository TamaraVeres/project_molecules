# LogP Prediction and Analysis Using RDKit

The aim of this project is to use RDKit to calculate logP values and compare them with experimental data. Linear regression and deep neural network models are then used to test whether better predictions can be achieved than with the RDKit calculation alone.

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
├── molecules_main.py              # Main script to run the full analysis pipeline
├── requirements.txt               # Project dependencies   
├── logP_dataset.csv              
├── README.md                     
│
├── src/                        
│   ├── utils.py                  
│   ├── molecular_descriptors.py 
│   ├── compare.py               
│   ├── linear_regression.py   
│   └── dnn.py                    
│
└── results/                      
    ├── comparison_results.csv    
    ├── comparison_results.html   
    ├── logp_results.png          
    ├── logp_linear_regression.png 
    └── logp_dnn.png              
```

---

## Data Inspection and Preparation

First, the CSV file is loaded and inspected to understand the dataset. The SMILES strings are then converted into RDKit molecule objects, and invalid molecules are removed.

## LogP Calculation with RDKit and Mordred

 LogP is calculated for each molecule using RDKit and Mordred. A comparison between RDKit and Mordred shows that both libraries produce similar descriptor values.

## Comparison with Experimental LogP Values 

Finally, the molecular descriptors calculated with RDKit are compared with the reference values from the dataset. For each descriptor, the experimental value, the RDKit value, and the error between them are calculated and stored in a table. The comparison results are saved as output files in the results folder. The comparison between RDKit-calculated and experimental logP values resulted in an R² of 0.55, showing that the similarity is not strong.

## Predictive Modeling: Linear Regression and DNN

To investigate whether predictive performance could be improved, a linear regression model and a deep neural network were trained using the available data. The linear regression model achieved an R² of 0.80, while the DNN further improved the performance to an R² of 0.84.

### Model Performance Results

| Metric | Linear Regression | Deep Neural Network |
|--------|------------------|---------------------|
| **Train MAE** | 0.345 | 0.117 |
| **Test MAE** | 0.435 | 0.348 |
| **Train MSE** | 0.215 | - |
| **Test MSE** | 0.340 | 0.261 |
| **RMSE** | 0.583 | 0.511 |
| **R²** | **0.802** | **0.848** |

The DNN model outperformed Linear Regression across all metrics, achieving the highest R² score of 0.848 and the lowest test error (RMSE: 0.511). Both models significantly improved upon the baseline RDKit calculation (R² = 0.553).

---

## Conclusion

In this project, RDKit was used to calculate logP values and compare them with experimental data. The comparison showed that RDKit logP does not closely match the experimental values. To improve the predictions, a linear regression model and a deep neural network were developed. Both models showed improved performance compared to the RDKit calculation alone, with the deep neural network giving the best results. However, experimental data is still needed for reliable logP values.
