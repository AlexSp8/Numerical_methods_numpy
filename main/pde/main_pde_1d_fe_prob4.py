
import sys
from pathlib import Path

import numpy as np

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from linear_systems import direct_solver
from ode import plot
from pde import bvp_setup

# Problem 4
def bc_left(t: float, x: float, u: np.ndarray[tuple[int]], dudx: np.ndarray[tuple[int]]
    ) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 0.0     # Dirichlet
    # res[1] = dudx[1] - 0.0  # Neumann: does not work due to cross-coupling
    res[1] = u[1] - np.exp(-t)
    return res

def bc_right(t: float, x: float, u: np.ndarray[tuple[int]], dudx: np.ndarray[tuple[int]]
    ) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    # res[0] = dudx[0] - 0.0  # Neumann: does not work due to cross-coupling
    res[0] = u[0] - np.exp(-t)
    res[1] = u[1] - 0.0     # Dirichlet
    return res

def initial_condition(t0: float, x: np.ndarray[tuple[int]], neq: int
    ) -> np.ndarray[tuple[int]]:

    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)

    u0[0::neq] = np.sin(x)
    u0[1::neq] = np.cos(x)

    return u0

def f_weak(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], w: float, dwdx: float) -> np.ndarray[tuple[int]]:

    neq = u.shape[0]
    res = np.zeros(neq)
    
    v, D, s = physical_quantities(t, x, u)

    res[0] = (dudt[0] - s[0])*w + D[0]*dudx[1]*dwdx
    res[1] = (dudt[1] - s[1])*w - D[1]*dudx[0]*dwdx
    return res

def f_strong(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], d2udx2: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

    neq = u.shape[0]
    res = np.zeros(neq)
    
    v, D, s = physical_quantities(t, x, u)

    res[0] = dudt[0] - (D[0]*d2udx2[1] + s[0])
    res[1] = dudt[1] - (-D[1]*d2udx2[0] + s[1])
    return res

def physical_quantities(t: float, x: float, u: np.ndarray[tuple[int]]):

    neq = u.shape[0]
    D = np.ones(neq)

    v = np.zeros(neq)

    s = np.zeros(neq)    

    ext = np.exp(-t)
    cos_x = np.cos(x)
    sin_x = np.sin(x)

    s[0] = ext*(cos_x - sin_x)
    s[1] = -ext*(cos_x + sin_x)

    return v, D, s

def stabilization_operator(t: float, x: float, u: np.ndarray[tuple[int]],
    w: float, dwdt: np.ndarray[tuple[int]], dwdx: float,
    d2wdx2: float) -> np.ndarray[tuple[int]]:

    neq = u.shape[0]
    p_mat = np.zeros((neq, neq))

    return p_mat

def u_analytical(t: float, x: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    nnodes = x.shape[0]
    neq = 2
    u_exact = np.zeros(nnodes*neq)
    u_exact[0::neq] = np.exp(-t)*np.sin(x)
    u_exact[1::neq] = np.exp(-t)*np.cos(x)
    return u_exact

def main():

    neq = 2

    bvp_solver = bvp_setup.BVPSetupFE(neq)

    xd = np.array([0.0, np.pi/2])
    nel, nbf = 5, 3
    nnodes = (nbf-1)*nel+1
    bvp_solver.create_mesh(nnodes=nnodes, xd=xd, p=1.0)
    x = bvp_solver.x

    D, v = 1.0, 0.0

    bvp_solver.check_dx(D, v)

    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 0.5, 0.01, 1e-3, 0.1, 1e-3, 1e-2
    # dt = bvp_solver.check_dt(D, v, dt)
    
    u0 = initial_condition(t0, x, neq)

    bvp_solver.set_nr_solver(u0, k_max=100, tol=1e-8, r=1.0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)
    
    ng = 3
    bvp_solver.set_fe_gauss_int(ng=ng, nbf=nbf)
    
    bvp_solver.set_fe_stabilization(f_quantities=physical_quantities, f_strong=f_strong,
                                    f_operator=stabilization_operator, order=2)

    bc = {'left': bc_left, 'right': bc_right}

    bvp_solver.set_fe_problem(f_weak, bc, theta=1.0)

    u = bvp_solver.solve(dtw=0.2)

    u_exact = u_analytical(tf, x)
    err = np.max(np.abs(u_exact-u))
    print(f'err = {err:.4e}')

    for ieq in range(neq):
        ui_num, ui_ex = u[ieq::neq].copy(), u_exact[ieq::neq].copy()
        plot.plot_ode_system_evolution(x, [ui_num, ui_ex])

if __name__ == '__main__':
    main()
