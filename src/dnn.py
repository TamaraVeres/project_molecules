
import numpy as np
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error
import torch
import torch.nn as nn
import torch.optim as optim
import matplotlib.pyplot as plt


def prepare_torch_data(X_train, y_train, X_test, y_test):

    X_train_tensor = torch.tensor(X_train, dtype=torch.float32)
    y_train_tensor = torch.tensor(
        y_train, dtype=torch.float32).view(-1, 1)

    X_test_tensor = torch.tensor(X_test, dtype=torch.float32)
    y_test_tensor = torch.tensor(
        y_test, dtype=torch.float32).view(-1, 1)

    return X_train_tensor, X_test_tensor, y_train_tensor, y_test_tensor


def deep_neural_network_model(input_size):
    model_sequential = nn.Sequential(
        nn.Linear(input_size, 512), # 1024 is the number of features
        nn.ReLU(), 
        nn.Dropout(0.3), # prevent overfitting
        nn.Linear(512, 256),
        nn.ReLU(),
        nn.Dropout(0.3),
        nn.Linear(256, 128),
        nn.ReLU(),
        nn.Dropout(0.2),
        nn.Linear(128, 64),
        nn.ReLU(),
        nn.Linear(64, 1),
    
    )
    return model_sequential


def train_model(model_sequential, X_train_tensor, y_train_tensor, num_epochs=1000, lr=0.001):
    loss_fn = nn.MSELoss()
    optimizer = optim.Adam(model_sequential.parameters(), lr=lr)

    for epoch in range(num_epochs):
        predictions = model_sequential(X_train_tensor) 
        MSE = loss_fn(predictions, y_train_tensor) 
        optimizer.zero_grad() 
        MSE.backward() 
        optimizer.step() 
        
        if (epoch + 1) % 100 == 0:
            print(f"Epoch {epoch+1}/{num_epochs}, MSE: {MSE.item():.4f}")
    return model_sequential



def evaluate_model(model_sequential, X_train_tensor, X_test_tensor, y_train_tensor, y_test_tensor):
    model_sequential.eval()
    loss_fn = nn.MSELoss()
    with torch.no_grad():  
        train_predictions = model_sequential(X_train_tensor) 
        test_predictions = model_sequential(X_test_tensor) 
        test_loss = loss_fn(test_predictions, y_test_tensor)

        print("\nNeural network model evaluation:")
        print("Train MAE:", mean_absolute_error(y_train_tensor.numpy(), train_predictions.numpy()))
        print("Test MAE:", mean_absolute_error(y_test_tensor.numpy(), test_predictions.numpy()))
        print("Train MSE:", mean_squared_error(y_train_tensor.numpy(), train_predictions.numpy()))
        print("Test MSE:", test_loss.item())
        print("RMSE:", np.sqrt(test_loss.item()))
        print("R²:", r2_score(y_test_tensor.numpy(), test_predictions.numpy()))

        
        nn_preds = test_predictions.numpy().flatten()
        return nn_preds

def scatter_plot_dnn(y_test, dnn_preds):
    plt.scatter(y_test, dnn_preds)
    m, b = np.polyfit(y_test, dnn_preds, 1)
    plt.plot(y_test, m*y_test + b, color='red')
    r2 = r2_score(y_test, dnn_preds)
    plt.text(0.05, 0.95, f"R² = {r2:.3f}", transform=plt.gca().transAxes, fontsize=12)
    plt.xlabel("Experimental XLogP")
    plt.ylabel("Predicted XLogP")
    plt.title("Experimental vs DNN Predicted XLogP")
    plt.savefig("./results/logp_dnn.png") 
    plt.show()