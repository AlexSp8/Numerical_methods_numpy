
from typing import List
import matplotlib.pyplot as plt

def plot_ode_evolution(x: List[float] = None, y: List[float] = None):
    """Plots the evolution of an ODE."""

    plt.plot(x, y, color='red', label='y(x)', marker='o')

    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('ODE Evolution')
    plt.legend()
    plt.grid(True, linestyle=':', alpha=0.6)

    plt.show()

def plot_ode_system_evolution(x: List[float] = None, y: List[List[float]] = None):
    """Plots the evolution of a system of ODEs."""

    neq = len(y)

    for j in range(neq):
        plt.plot(x, y[j], marker='o', label=f'y_{j+1}')

    plt.title("Solution of Coupled BVP System")
    plt.xlabel("x")
    plt.ylabel("y variables")
    plt.grid(True, linestyle='--', alpha=0.7)
    plt.legend()
    plt.show()
