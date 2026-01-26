import pandas as pd
from src.utils import load_csv, basic_inspection, clean_molecules, prepare_data
from src.molecular_descriptors import calculate_rdkit_descriptors, distribution, calculate_mordred_descriptors, calculate_morgan_fingerprint
from src.compare import compare_descriptors, scatter_plot
from sklearn.model_selection import train_test_split
from src.linear_regression import linear_regression_model, scatter_plot_linear_regression
from src.dnn import prepare_torch_data, deep_neural_network_model, train_model, evaluate_model, scatter_plot_dnn
from rdkit import RDLogger
RDLogger.DisableLog('rdApp.*') # disable RDKit logging

def molecules_main():
    df_compounds = load_csv(file_path="./logP_dataset.csv")
    basic_inspection(df_compounds)
    mols = clean_molecules(df_compounds)
    df_molecular_descriptors = calculate_rdkit_descriptors(mols)
    distribution(df_molecular_descriptors)
    df_mordred_descriptors = calculate_mordred_descriptors(mols)
    comparison_table = compare_descriptors(df_molecular_descriptors, df_mordred_descriptors, df_compounds)
    print(comparison_table.head(10))
    scatter_plot(df_compounds, df_molecular_descriptors)
    comparison_table.to_csv("./results/comparison_results.csv", index=False)
    comparison_table.to_html("./results/comparison_results.html", index=False)
    
    fingerprints = calculate_morgan_fingerprint(mols)
    X, y = prepare_data(fingerprints, df_compounds)
    X_train, X_test, y_train, y_test = train_test_split(X, y, test_size=0.2, random_state=42)
    linear_test_predictions = linear_regression_model(X_train, y_train, X_test, y_test)
    X_train_tensor, X_test_tensor, y_train_tensor, y_test_tensor = prepare_torch_data(
        X_train, y_train, X_test, y_test)
    dnn_model = deep_neural_network_model(X_train_tensor.shape[1])
    dnn_model = train_model(dnn_model, X_train_tensor, y_train_tensor)
    dnn_preds = evaluate_model(dnn_model, X_train_tensor, X_test_tensor, y_train_tensor, y_test_tensor)
    scatter_plot_linear_regression(y_test, linear_test_predictions)
    scatter_plot_dnn(y_test, dnn_preds)
  
   # compare the predicted logP by linear regression and deep neural network
    predicted_logp = pd.DataFrame({
        "Actual logP": y_test,
        "Linear Regression Prediction": linear_test_predictions,
        "Deep Neural Network Prediction": dnn_preds
    })
    print("Predicted logP:\n", predicted_logp.head(10))

  
if __name__ == "__main__":
    molecules_main()