
import sys
from pathlib import Path

import numpy as np
import numpy.typing as npt

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from optimization.constrained import equality_constraints
from differentiation import forward_fd as ffd

def f(xp: npt.NDArray[np.float64]) -> float:
    x, y = xp[0], xp[1]
    # return x**2 + y**2
    # return x + y
    # return (x+y)**2
    return (x**2)*y

def f_max(xp: npt.NDArray[np.float64]) -> float:
    return -f(xp)

def g_con(xp: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    x, y = xp[0], xp[1]
    # return np.array([x + y - 1])
    # return np.array([x**2 + y**2 - 1])
    # return np.array([x**2 + y**2 - 1])
    return np.array([x**2 + y**2 - 3])

def main():

    print('\nConstrained optimization')

    print('Lagrange Multipliers')
    x0, l0_m = np.array([1.0, 1.0]), np.array([1.0])
    df = ffd.df_h
    x, l = equality_constraints.lagrange_multipliers(f, g_con, x0, df=df, l0_m=l0_m)
    print(f"x = {x}, l = {l}, f = {f(x)}")

    print('Augmented Lagrangian')
    x0 = np.array([1.0, 1.0])
    df = ffd.df_h
    x, l = equality_constraints.augmented_lagrangian(f_max, g_con, x0, df=df)
    print(f"x = {x}, l = {l}, f = {f(x)}")

if __name__ == '__main__':
    main()
