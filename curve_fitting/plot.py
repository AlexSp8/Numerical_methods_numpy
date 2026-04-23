
from typing import List
import matplotlib.pyplot as plt

def plot_experiment_model(xi: List[float] = None,
    yi: List[float] = None, y_model: List[float] = None):
    """Plots the predictions of the model along with
    experimental data."""

    plt.scatter(xi, yi, color='red', label='Experiment', marker='o')

    plt.plot(xi, y_model, color='blue', linestyle='-', label='Model')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Experiment vs Model Predictions')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)

    plt.show()
