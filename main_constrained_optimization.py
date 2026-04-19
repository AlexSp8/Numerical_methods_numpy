
import numpy as np
import numpy.typing as npt

from optimization.constrained import equality_constraints

def f(xp: npt.NDArray[np.float64]) -> float:
    x, y = xp[0], xp[1]
    # return x**2 + y**2
    # return x + y
    # return (x+y)**2
    return (x**2)*y

def g_con(xp: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
    x, y = xp[0], xp[1]
    # return np.array([x + y - 1])
    # return np.array([x**2 + y**2 - 1])
    # return np.array([x**2 + y**2 - 1])
    return np.array([x**2 + y**2 - 3])

def main():

    print('\nConstrained optimization')

    print('Lagrange Multipliers')
    x0, l0 = np.array([1.0, 1.0]), np.array([1.0])
    x, l = equality_constraints.lagrange_multipliers(f, g_con, x0, l0, fd_type='ffd')
    print(f"x = {x}, l = {l}, f = {f(x)}")

    print('Augmented Lagrangian')
    x0 = np.array([1.0, 1.0])
    x, l = equality_constraints.augmented_lagrangian(f, g_con, x0, mode='max', output=True)
    print(f"x = {x}, l = {l}, f = {f(x)}")

if __name__ == '__main__':
    main()
