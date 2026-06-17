
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

# Time-dependent BCs
def bc_left(t: float, x_b: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    # Exact value at x=0 is sin(t)*(0) = 0
    return u_b - 0.0

def bc_right(t: float, x_b: float, u_b: np.ndarray[tuple[int]],
    dudx_b: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    # Exact value at x=2 is sin(t)*(4 + 2) = 6*sin(t)
    return u_b - 6.0*np.sin(t)

def initial_condition(t0: float, x: np.ndarray[tuple[int]], neq: int
    ) -> np.ndarray[tuple[int]]:
    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)
    return u0

def f_res(t: float, x: float, u: np.ndarray[tuple[int]],
    dudt: np.ndarray[tuple[int]], dudx: np.ndarray[tuple[int]],
    d2udx2: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    D, v = 1.0, 2.0
    S = np.cos(t)*(x**2 + x) + 4.0*x*np.sin(t)
    return dudt - (D*d2udx2 - v*dudx + S)

def u_analytical(t: float, x: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    return np.sin(t)*(x**2+x)

def main():

    neq = 1
    
    bvp_solver = bvp_setup.BVPSetupFD(neq)

    xd = np.array([0.0, 2.0])
    bvp_solver.create_mesh(nnodes=31, xd=xd, p=1.0)
    x = bvp_solver.x

    D, v = 1.0, 2.0

    bvp_solver.check_dx(D, v)

    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 2.0, 0.1, 1e-2, 0.5, 1e-2, 1e-1
    # dt = bvp_solver.check_dt(D, v, dt)

    u0 = initial_condition(t0, x, neq)
    
    bvp_solver.set_nr_solver(u0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)

    bc = {'left': bc_left, 'right': bc_right}

    bvp_solver.set_fd_problem(f_res, bc)

    u = bvp_solver.solve(dtw=1.0)

    u_exact = u_analytical(tf, x)
    err = np.max(np.abs(u_exact-u))
    print(f'err = {err:.4e}')

    for ieq in range(neq):
        ui_num, ui_ex = u[ieq::neq].copy(), u_exact[ieq::neq].copy()
        plot.plot_ode_system_evolution(x, [ui_num, ui_ex])

if __name__ == '__main__':
    main()
