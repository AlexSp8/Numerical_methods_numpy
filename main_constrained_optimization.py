
from typing import List

from optimization.constrained import equality_constraints

def f_2D(xp: List[float]) -> float:
    x, y = xp[0], xp[1]
    # return x**2 + y**2
    # return x + y
    # return (x+y)**2
    return (x**2)*y

def g_constraints(xp: List[float]) -> List[float]:
    x, y = xp[0], xp[1]
    # return [x + y - 1]
    # return [x**2 + y**2 - 1]
    # return [x**2 + y**2 - 1]
    return [x**2 + y**2 - 3]

def main():

    print('\nConstrained optimization')
    print('Lagrange Multipliers')
    x0, l0 = [1.0, 1.0], [1.0]
    print(equality_constraints.lagrange_multipliers(f_2D, g_constraints, x0, l0, mode='max'))

if __name__ == '__main__':
    main()
