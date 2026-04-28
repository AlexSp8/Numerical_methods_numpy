# Numerical Methods Library

A comprehensive Python library implementing classical numerical methods for solving mathematical problems, designed for educational and research purposes.

## Features

This library provides implementations of fundamental numerical algorithms:

### Root Finding (`roots`)
- **Bracketing Methods**: Bisection, False Position, Multi-bracketing
- **Open Methods**: Fixed Point, Secant, IQI, Newton-Raphson, Ralston-Rabinowitz, Multi-Open
- **Hybrid methods**: Brent, Ridders, Chandrupatla, Multi-Hybrid

### Linear Systems (`linear_systems`)
- Linear solver abstract class
- **Direct Solver classes**: Cramer's Rule, Gaussian Elimination, LU Decomposition, Thomas Algorithm (tri-diagonal routine)
- **Iterative Solver classes**: Jacobi, Successive Over-Relaxation (SOR), Steepest Descent, Conjugate Gradient

### Non-Linear Systems (`non_linear_systems`)
- Newton-Raphson Solver class
- Levenberg-Marquardt Solver class

### Optimization (`optimization`)
#### 1D optimization (`one_dimensional`)
- **Bracketing Methods**: Golden-section search, Multi-Bracketing
- **Open Methods**: Parabolic Interpolation, Secant, Newton-Raphson, Multi-open
- **Hybrid Methods**: Brent, Multi-hybrid

#### Unconstrained
- **Direct Methods**: Powell
- **Gradient Methods**: Steepest Descent, Conjugate Gradient, Newton, Marquardt, BFGS
- Line search helper class for 1D optimization.

#### Constrained
- **Equality Constraints Methods**: Lagrange multipliers, Augmented Lagrangian

### Curve Fitting (`curve_fitting`)
- **Regression**: Linear, Polynomial, Multi-linear, Non-linear (Levenberg-Marquardt)

- **Interpolation**: Vandermonde, Newton-Gregory polynomials, Lagrange polynomials, Splines (linear, quadratic, cubic)

- **Fourier Transform**: Discrete Fourier Transform (DFT), Fast Fourier Transform (FFT) Sande-Tukey

### Integration
- **Function Integration**: Trapezoidal, Simpson 1/3, Simpson 3/8, Mixed Simpson, Romberg, Gauss-Legendre (1D, 2D, multi-dimensional, points-weights coefficients evaluation)
- **Data Integration**: Trapezoidal, Simpson 1/3, Simpson Mixed (1D, 2D, multi-dimensional)

### Differentiation
- **Finite Differences**: Forward ($O(h), O(h^2)$), Backward ($O(h), O(h^2)$), Central ($O(h^2), O(h^4)$) up to fourth order.
- **Partial Derivatives**: Grad, Hessian.
- **Function Differentiation**: Richardson Extrapolation
- **Data Differentiation**: Divided Differences

### ODEs
#### Explicit Methods
**RK1**: Explicit Euler \
**RK2**: General method, Heun's method, Midpoint \
**RK3**: General method (+coefficients evaluation) \
**RK4**: Classic RK4 method \
**RK5**: Butcher's RK5 method \
**Adaptive step methods**: RK45 Fehlberg, RK45 Cash-Karp, RK45 Dormand-Prince \
**Adams-Bashforth**: Explicit multistep \
**Implicit**: Euler, Midpoint, Crank-Nicolson, Gauss-Legendre (2- and 3-stage), Radau-IIA (2- and 3-stage), Lobatto IIIC-2, Lobatto IIIA-3, Adams-Moulton \
**BVP 1D**: $2^{nd}$ order non-linear BVP 1D
**EVP 1D**: $2^{nd}$ order linear EVP 1D (Sturm-Liouville)

### Utilities (`utilities`)
- Matrix operations (addition, multiplication, transpose, determinant, inverse, etc.)
- Vector operations (norm, dot product)
- Statistics (mean, STD, variance, etc.)
- Indexing: nearest index, starting index, Bubble Sort
- I/O: read, write files

## Installation

### Requirements
- Python 3.8 or higher
- Standard library
- numpy (2.4.4)
- scipy (1.17.1) for built-in functions
- matplotlib (3.10.8) for graphical solutions

### Install from source
```bash
git clone <repository-url>
cd numerical-methods
pip install .
```

### Direct usage
Simply clone the repository and import the modules directly:
```python
from roots.bracketing import bisection
from linear_systems.direct_solvers import gauss_elimination
```

## Quick Start

### Root Finding Example
```python
from roots.bracketing import bisection
import math

# Find root of f(x) = x^2 - 2 in [1, 2]
def f(x):
    return x**2 - 2

root = bisection(f, 1.0, 2.0)
print(f"Root: {root}")  # Should be approximately 1.4142
```

### Linear System Example
```python
from linear_systems.direct_solvers import gauss_elimination

# Solve Ax = b where A is 3x3 matrix, b is vector
A = [[2, 1, -1], [1, 3, 2], [-1, 2, 4]]
b = [8, 13, 10]

x = gauss_elimination(A, b)
print(f"Solution: {x}")
```

### Non-Linear System Example
```python
from non_linear_systems.non_linear_system_iterative_solvers import newton_raphson
import math

# Solve system: x^2 + y^2 = 1, x + y = 1
def F(xy):
    x, y = xy
    return [x**2 + y**2 - 1, x + y - 1]

x0 = [0.5, 0.5]  # Initial guess
solution = newton_raphson(F, x0)
print(f"Solution: {solution}")
```

## Project Structure

```
numerical_methods/
├── main_*.py          # Demos
├── roots/
│   ├── bracketing.py
│   ├── hybrid.py
│   └── open_methods.py
├── linear_systems/
│   ├── direct_solvers.py
│   └── iterative_solvers.py
├── non_linear_systems/
│   └── non_linear_system_iterative_solvers.py
├── optimization/
│   ├── constrained/
│       └── equality_constraints.py
│   ├── unconstrained_1D/
│       ├── bracketing.py
│       └── open_methods.py
│   └── unconstrained_multi/
│       ├── direct_methods.py
│       └── gradient_methods.py
├── curve_fitting/
│       ├── fourier_transform.py
│       ├── interpolation.py
│       ├── plot.py
│       └── regression.py
├── integration/
│       ├── data_integration.py
│       ├── function_integration.py
│       └── gauss_legendre_data.py
├── differentiation/
│       ├── backward_fd.py
│       ├── central_fd.py
│       ├── forward_fd.py
│       ├── data_differentiation.py
│       ├── function_differentiation.py
│       └── partial_derivatives.py
├── ode/
│       ├── explicit_ode.py
│       ├── implicit_ode.py
│       ├── butcher_tableaus.py
│       ├── bvp_1d.py
│       ├── evp_1d.py
│       └── plot.py

└── utilities/
    ├── indexing.py             # Index finding, sort
    ├── io_utils.py             # Read/write files
    ├── matrix_operations.py    # Matrix algebra
    ├── statistics.py           # Statistics quantities
    └── vector_operations.py    # Vector operations
```

## Educational Purpose

This library is designed for:
- Learning numerical methods algorithms
- Comparing different solution approaches
- Understanding convergence behavior
- Research and prototyping

## Future Development

- [ ] Add PDEs module
- [ ] Add more advanced algorithms
