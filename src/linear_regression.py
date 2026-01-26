import numpy as np
from sklearn.linear_model import LinearRegression
from sklearn.metrics import mean_squared_error, r2_score, mean_absolute_error
import matplotlib.pyplot as plt

def linear_regression_model(X_train, y_train, X_test, y_test):
    linear_model = LinearRegression()
    linear_model.fit(X_train, y_train)

    linear_train_predictions = linear_model.predict(X_train)
    linear_train_mse = mean_squared_error(y_train, linear_train_predictions)
    linear_test_predictions = linear_model.predict(X_test)
    test_mse = mean_squared_error(y_test, linear_test_predictions)
    
    print("Linear regression model evaluation:")
    print("Train MAE:", mean_absolute_error(y_train, linear_train_predictions))
    print("Test MAE:", mean_absolute_error(y_test, linear_test_predictions))
    print("Train MSE:", linear_train_mse)
    print("Test MSE:", test_mse)
    print("RMSE:", np.sqrt(test_mse))
    print("R^2:", r2_score(y_test, linear_test_predictions))
    
    return linear_test_predictions
    
    
def scatter_plot_linear_regression(y_test, linear_test_predictions):
    plt.scatter(y_test, linear_test_predictions)
    m, b = np.polyfit(y_test, linear_test_predictions, 1)
    plt.plot(y_test, m*y_test + b, color='red')
    r2 = r2_score(y_test, linear_test_predictions)
    plt.text(0.05, 0.95, f"R² = {r2:.3f}", transform=plt.gca().transAxes, fontsize=12)
    plt.xlabel("Experimental XLogP")
    plt.ylabel("Predicted XLogP")
    plt.title("Experimental vs Linear Regression Predicted XLogP")
    plt.savefig("./results/logp_linear_regression.png") 
    plt.show()