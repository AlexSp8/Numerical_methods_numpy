
from typing import Callable

from linear_systems import direct_solver, iterative_solver
from utilities import matrix_operations
from ode import butcher_tableaus, explicit_ode

def implicit_solver(f: Callable[[float, List[float]], List[float]],
    x0: float, xf: float, h: float, y0: List[float],
    method: str = 'euler') -> tuple[List[float], List[List[float]]]:
    """Returns the solution of a system of ODEs using a single-stage implicit method."""

    n = int(round((xf-x0)/h))
    neq = len(y0)
    x = [x0+i*h for i in range(n+1)]
    y = [ [0.0 for _ in range(neq)] for _ in range(n+1) ]
    y[0] = y0[:]

    method_dict = butcher_tableaus.IMPLICIT_METHODS[method]
    ns = method_dict['ns']
    b = method_dict['b']
    p = method_dict['p']

    for i in range(n):

        xn, yi = x[i+1], y[i]

        u = newton_raphson(f, xn, yi, h, method_dict)

        y_next = yi[:]
        for j in range(ns):
            Y_j = u[j*neq:(j+1)*neq]
            t_j = (xn - h) + p[j]*h
            k_j = f(t_j, Y_j)
            for ieq in range(neq):
                y_next[ieq] += h*b[j]*k_j[ieq]

        y[i+1] = y_next

    return x, y

def newton_raphson(f: Callable[[float, List[float]], List[float]],
    xn: float, yo: List[float], h: float, method_dict: dict,
    eps: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> List[float]:
    """Returns the solution of a non-linear system of algebraic equations, F,
    around an initial guess, u0, using the Newton-Raphson method"""

    ns = method_dict['ns']
    n = ns*len(yo)

    u = yo*ns

    for k in range(1, k_max+1):

        res = residual(f, u, xn, yo, h, method_dict)
        res_norm = vector_operations.norm2(res)

        jac = jacobian(u, res, f, xn, yo, h, method_dict)
        # print(matrix_operations.condition_number(jac))

        b = [-val for val in res]

        du = direct_solvers.lu_decomposition_solve(jac, b)
        # du = iterative_solvers.sor_solver(jac, b, x0=None, w=0.5, output=False)

        cor_norm = vector_operations.norm2(du)

        # print(f'k = {k}, Res Norm: {res_norm:.4e}, Cor Norm: {cor_norm:.4e}')
        if cor_norm < eps and res_norm < eps:
            break

        u = [ u[i]+r*du[i] for i in range(n) ]

    return u

def jacobian(u: List[float], res: List[float],
    f: Callable[[float, List[float]], List[float]], xn: float,
    yo: List[float], h: float, method_dict: dict, e: float = 1e-8) -> List[List[float]]:
    """Returns the Jacobian of a system of non-linear algebraic equations
    around the values, u, of the unknowns."""

    n, m = len(u), len(res)

    jac = [ [0.0 for _ in range(n)] for _ in range(m) ]
    for j in range(n):

        u_f = u[:]
        u_f[j] += e
        res_f = residual(f, u_f, xn, yo, h, method_dict)

        # u_b = u[:]
        # u_b[j] -= e
        # #res_b = residual(f, u_b, xn, yo, h, method_dict)

        for i in range(m):
            jac[i][j] = (res_f[i]-res[i])/e
            # jac[i][j] = (res_f[i]-res_b[i])/(2*e)

    return jac

def residual(f: Callable[[float, List[float]], List[float]],
    yn: List[float], xn: float, yo: List[float], h: float,
    method_dict: dict) -> List[float]:
    """Returns the residual of an implicit method given the
    current and previous states and the step, h."""

    neq = len(yo)
    ns = method_dict['ns']
    q = method_dict['q']
    p = method_dict['p']

    Ys = [yn[i*neq:(i+1)*neq] for i in range(ns)]

    t_i = xn - h
    k_j = [0.0]*ns
    for j in range(ns):
        k_j[j] = f(t_i + p[j]*h, Ys[j])

    res = [0.0]*(ns*neq)
    for j in range(ns):
        for ieq in range(neq):
            row = j*neq + ieq
            s_sum = 0.0
            for m in range(ns):
                s_sum += q[j][m]*k_j[m][ieq]
            res[row] = Ys[j][ieq] - yo[ieq] - h*s_sum

    return res

def adams_moulton(f: Callable[[float, List[float]], List[float]],
    x0: float, xf: float, h: float, y0: List[float],
    order: int = 2) -> tuple[List[float], List[List[float]]]:
    """Returns the solution of a system of ODEs using an
    Adams-Moulton implicit method of specified order."""

    AM_COEFFS = {
        2: [0.5, 0.5],
        3: [5.0/12.0, 8.0/12.0, -1.0/12.0],
        4: [9.0/24.0, 19.0/24.0, -5.0/24.0, 1.0/24.0],
        5: [251.0/720.0, 646.0/720.0, -264.0/720.0,
            106.0/720.0, -19.0/720.0]
    }

    if order not in AM_COEFFS:
        raise ValueError(f"Order {order} is not supported.")

    b = AM_COEFFS[order]

    n = int(round((xf-x0)/h))
    if n < order:
        raise ValueError(f"Order {order} larger than n.")

    neq = len(y0)
    x = [x0+i*h for i in range(n+1)]
    y = [ [0.0 for _ in range(neq)] for _ in range(n+1) ]
    y[0] = y0[:]

    for i in range(order-1):
        xi = x[i]
        _, y_rk = explicit_ode.rk4(f, xi, xi+h, h, y[i])
        y[i+1] = y_rk[-1]

    f_past = [f(x[i], y[i]) for i in range(order)]

    for i in range(order-1, n):
        y_next = y[i][:]
        # for k in range(order-1):
        for k in range(1,order):
            # f_j = f_past[-(k+1)]
            f_j = f_past[-k]
            for ieq in range(neq):
                y_next[ieq] += h*b[k]*f_j[ieq]

        y[i+1] = newton_raphson_am(f, x[i+1], y_next, h, b[0])
        if i+1 < n:
            f_past = f_past[1:] + [f(x[i+1], y[i+1])]

    return x, y

def newton_raphson_am(f: Callable[[float, List[float]], List[float]],
    xn: float, yo: List[float], h: float, b0: float,
    eps: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> List[float]:
    """Returns the solution of a non-linear system of algebraic equations, F,
    around an initial guess, u0, using the Newton-Raphson method"""

    n = len(yo)

    u = yo[:]

    for k in range(1, k_max+1):

        res = residual_am(f, u, xn, yo, h, b0)
        res_norm = vector_operations.norm2(res)

        jac = jacobian_am(u, res, f, xn, yo, h, b0)
        # print(matrix_operations.condition_number(jac))

        b = [-val for val in res]

        du = direct_solvers.lu_decomposition_solve(jac, b)
        # du = iterative_solvers.sor_solver(jac, b, x0=None, w=0.5, output=False)

        cor_norm = vector_operations.norm2(du)

        # print(f'k = {k}, Res Norm: {res_norm:.4e}, Cor Norm: {cor_norm:.4e}')
        if cor_norm < eps and res_norm < eps:
            break

        u = [ u[i]+r*du[i] for i in range(n) ]

    return u

def jacobian_am(u: List[float], res: List[float],
    f: Callable[[float, List[float]], List[float]], xn: float,
    yo: List[float], h: float, b0: float, e: float = 1e-8) -> List[List[float]]:
    """Returns the Jacobian of a system of non-linear algebraic equations
    around the values, u, of the unknowns."""

    n, m = len(u), len(res)

    jac = [ [0.0 for _ in range(n)] for _ in range(m) ]
    for j in range(n):

        u_f = u[:]
        u_f[j] += e
        res_f = residual_am(f, u_f, xn, yo, h, b0)

        # u_b = u[:]
        # u_b[j] -= e
        # #res_b = residual_am(f, u_b, xn, yo, h, method_dict)

        for i in range(m):
            jac[i][j] = (res_f[i]-res[i])/e
            # jac[i][j] = (res_f[i]-res_b[i])/(2*e)

    return jac

def residual_am(f: Callable[[float, List[float]], List[float]],
    yn: List[float], xn: float, yo: List[float], h: float,
    b0: float) -> List[float]:
    """Returns the residual of an implicit method given the
    current and previous states and the step, h."""

    neq = len(yo)

    f_val = f(xn,yn)

    res = [0.0]*neq
    for ieq in range(neq):
        res[ieq] = yn[ieq] - yo[ieq] - h*b0*f_val[ieq]

    return res
