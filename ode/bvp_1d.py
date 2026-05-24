
from typing import Callable, List, Tuple, Dict, Any
import math

from linear_systems import direct_solvers, iterative_solvers
from utilities import vector_operations, matrix_operations

def create_mesh(nnodes: int, xd: List[float],
    p: int = 1) -> List[float]:
    """Returns the mesh coordinates given the number of nodes,
    the domain and the stretching factor, p, for power law
    refinement. p > 1: refinement towards x0,
    0 < p < 1: refinement towards xf."""

    x0, xf = xd[0], xd[1]
    L = xf - x0

    x = [0.0]*nnodes
    for i in range(nnodes):
        xi = i/(nnodes-1)
        gi = xi**p
        # gi = 0.5*(1.0-math.cos(p*math.pi*xi))
        # gi = math.tanh(p*(2*xi-1))/(2*math.tanh(p)) + 0.5
        x[i] = x0 + L*gi

    return x

def newton_raphson_bvp_1d(f: Callable[[float, List[float], List[float]], List[float]],
    x: List[float], bc: Dict[str, Dict[str, Any]],
    eps: float = 1e-8, k_max: int = 1000, r: float = 1.0) -> List[float]:
    """Returns the solution of the non-linear system of algebraic equations
    that results from FD approximation of a second order 1D-BVP
    (system of ODEs) problem using the Newton-Raphson method"""

    nnodes = len(x)
    neq = len(bc['left']['value'])

    nunknowns = nnodes*neq

    u = [1.0]*nunknowns

    for k in range(1, k_max+1):

        res = residual(f, x, u)

        jac = jacobian(res, f, x, u)

        for j in range(neq):

            h = x[1]-x[0]
            res[j], jac[j] = left_boundary_condition(u, j, bc['left'], h)

            row = (nnodes-1)*neq+j
            h = x[-1]-x[-2]
            res[row], jac[row] = right_boundary_condition(u, j, bc['right'], h)

        res_norm = vector_operations.norm2(res)
        # print(matrix_operations.condition_number(jac))

        b = [-val for val in res]

        du = direct_solvers.lu_decomposition_solve(jac, b)
        # du = iterative_solvers.sor_solver(jac, b, x0=None, w=0.5, output=False)

        cor_norm = vector_operations.norm2(du)

        print(f'k = {k}, Res Norm: {res_norm:.4e}, Cor Norm: {cor_norm:.4e}')
        if cor_norm < eps and res_norm < eps:
            break

        u = [ u[i]+r*du[i] for i in range(nunknowns) ]

    return u

def residual(f: Callable[[float, List[float], List[float]], List[float]],
    x: List[float], u: List[float]) -> List[float]:
    """Returns the residual of a second order 1D-BVP
    (system of ODEs) problem."""

    nnodes, nunknowns = len(x), len(u)

    neq = nunknowns//nnodes

    res = [0.0]*nunknowns
    for i in range(1,nnodes-1):

        hf = x[i+1] - x[i]
        hb = x[i] - x[i-1]

        row1 = i*neq

        ui = u[row1:row1+neq]

        t1_hb = 1.0/hb
        t1_hf = 1.0/hf
        t2_hf_hb = 2.0/(hf+hb)

        dui = [0.0]*neq
        for jeq in range(neq):
            row = row1 + jeq
            dui[jeq] = (u[row+neq] - u[row])*(t1_hf)
            # dui[jeq] = (u[row+neq] - u[row-neq])/(hb+hf)

        fi = f(x[i],ui,dui)

        for jeq in range(neq):
            row = row1 + jeq
            term2 = (u[row+neq] - u[row])*t1_hf
            term2 -= (u[row] - u[row-neq])*t1_hb
            term2 *= t2_hf_hb
            res[row] = term2 + fi[jeq]

    return res

def jacobian(res: List[float],
    f: Callable[[float, List[float], List[float]], List[float]],
    x: float, u: List[float],
    e: float = 1e-8) -> List[List[float]]:
    """Returns the Jacobian of a system of non-linear algebraic equations
    around the values, u, of the unknowns."""

    n, m = len(u), len(res)

    jac = [ [0.0 for _ in range(n)] for _ in range(m) ]
    for j in range(n):

        e = 1e-8*(abs(u[j]) + 1.0)

        u_f = u[:]
        u_f[j] += e
        res_f = residual(f, x, u_f)

        # u_b = u[:]
        # u_b[j] -= e
        # #res_b = residual(f, x, u_b)

        for i in range(m):
            jac[i][j] = (res_f[i]-res[i])/e
            # jac[i][j] = (res_f[i]-res_b[i])/(2*e)

    return jac

def left_boundary_condition(u: List[float], ieq: int,
    bc: Dict[str, Any], hf: float) -> Tuple[float, List[float]]:
    """Returns the residual and Jacobian entries after applying
    a boundary condition on the left boundary of a 1D domain."""

    neq, nunknowns = len(bc['value']), len(u)

    res_bc = 0.0
    jac_bc = [0.0]*nunknowns

    node1, node2 = ieq, ieq+neq

    bc_type = bc['type']
    val = bc['value'][ieq]

    if bc_type == 'dirichlet':
        jac_bc[node1] = 1.0
        res_bc = u[node1] - val
    elif bc_type == 'neumann':
        jac_bc[node1] = -1.0/hf
        jac_bc[node2] = +1.0/hf
        res_bc = (u[node2] - u[node1])/hf - val
    elif bc_type == 'robin':
        a0 = bc.get('a_robin', [0.0]*neq)[ieq]
        jac_bc[node1] = a0 - 1.0/hf
        jac_bc[node2] = +1.0/hf
        res_bc = (u[node2] - u[node1])/hf +a0*u[node1] - val

    return res_bc, jac_bc

def right_boundary_condition(u: List[float], ieq: int,
    bc: Dict[str, Any], hb: float) -> Tuple[float, List[float]]:
    """Returns the residual and Jacobian entries after applying
    a boundary condition on the right boundary of a 1D domain."""

    neq, nunknowns = len(bc['value']), len(u)

    res_bc = 0.0
    jac_bc = [0.0]*nunknowns

    node1 = nunknowns-neq+ieq
    node2 = node1-neq

    bc_type = bc['type']
    val = bc['value'][ieq]

    if bc_type == 'dirichlet':
        jac_bc[node1] = 1.0
        res_bc = u[node1] - val
    elif bc_type == 'neumann':
        jac_bc[node1] = +1.0/hb
        jac_bc[node2] = -1.0/hb
        res_bc = (u[node1] - u[node2])/hb - val
    elif bc_type == 'robin':
        an = bc.get('a_robin', [0.0]*neq)[ieq]
        jac_bc[node1] = an+1.0/hb
        jac_bc[node2] = -1.0/hb
        res_bc = (u[node1] - u[node2])/hb +an*u[node1] - val

    return res_bc, jac_bc
