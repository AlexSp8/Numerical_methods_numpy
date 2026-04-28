
from pathlib import Path

import numpy as np
import numpy.typing as npt

from differentiation import forward_fd as ffd, backward_fd as bfd
from differentiation import central_fd as cfd, function_differentiation
from differentiation import data_differentiation, partial_derivatives

# Path to the Datasets directory
DATASETS_DIR = Path(__file__).parent / 'integration/data'

def f(x: float) -> float:
    return 1.2 - 0.25*x - 0.5*(x**2) - 0.15*(x**3) - 0.1*(x**4)

def f_2D(xi: npt.NDArray) -> float:
    # x, y = xi[0], xi[1]
    # return 2*x*y + 2*x - x**2 -2*(y**2) + 72
    return np.sum(xi**2)

def main():

    print('Function Differentiation')
    print('\nForward')
    print(ffd.df_h2(f, x=0.5, h=0.25))
    print('\nBackward')
    print(bfd.df_h2(f, x=0.5, h=0.25))
    print('\nCentral')
    print(cfd.df_h4(f, x=0.5, h=0.25))
    print('\nRichardson')
    print(function_differentiation.richardson_extrapolation(cfd.df_h2, f, x=0.5, h=0.5))

    print('\nPartial Function Differentiation')
    x = np.array([1.0, -2.0, 3.0])
    df = cfd.df_h2
    print(partial_derivatives.grad_f_fd(df, f_2D, x))
    print(partial_derivatives.hessian_f_fd(df, f_2D, x))

    print('\nData Differentiation')
    x = np.array([0.0, 1.25, 3.75])
    y = np.array([10.0, 12.0, 13.5])
    print(data_differentiation.divided_difference(x, y))

if __name__ == '__main__':
    main()
