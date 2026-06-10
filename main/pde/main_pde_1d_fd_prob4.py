
import sys
from pathlib import Path

import numpy as np

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from linear_systems import direct_solver
from ode import plot
from pde import bvp_fd

# Problem 4
def bc_left(t: float, x: float, u: np.ndarray[tuple[int]], dudx: np.ndarray[tuple[int]]
    ) -> np.ndarray[tuple[int]]:
    res0 = u[0] - 0.0     # Dirichlet
    res1 = dudx[1] - 0.0  # Neumann
    return np.array([res0, res1])

def bc_right(t: float, x: float, u: np.ndarray[tuple[int]], dudx: np.ndarray[tuple[int]]
    ) -> np.ndarray[tuple[int]]:
    res0 = dudx[0] - 0.0  # Neumann
    res1 = u[1] - 0.0     # Dirichlet
    return np.array([res0, res1])

def initial_condition(t0: float, x: np.ndarray[tuple[int]], neq: int
    ) -> np.ndarray[tuple[int]]:

    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)

    u0[0::neq] = np.sin(x)
    u0[1::neq] = np.cos(x)

    return u0

def f_res(t: float, x: float, u: np.ndarray[tuple[int]], dudt: np.ndarray[tuple[int]],
    dudx: np.ndarray[tuple[int]], d2udx2: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

    ext = np.exp(-t)
    cos_x = np.cos(x)
    sin_x = np.sin(x)

    S0 = ext*(cos_x - sin_x)
    S1 = -ext*(cos_x + sin_x)

    res0 = dudt[0] - (d2udx2[1] + S0)
    res1 = dudt[1] - (-d2udx2[0] + S1)
    return np.array([res0, res1])

def u_analytical(t: float, x: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
    nnodes = x.shape[0]
    neq = 2
    u_exact = np.zeros(nnodes*neq)
    u_exact[0::neq] = np.exp(-t)*np.sin(x)
    u_exact[1::neq] = np.exp(-t)*np.cos(x)
    return u_exact

def main():
    
    neq = 2
    
    bvp_solver = bvp_fd.BVPSolver()

    xd = np.array([0.0, np.pi/2])
    bvp_solver.create_mesh(nnodes=11, xd=xd, p=1.0)
    x = bvp_solver.x

    D, v = 1.0, 0.0

    bvp_solver.check_dx(D, v)

    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 0.5, 0.01, 1e-3, 0.1, 1e-3, 1e-2
    # dt = bvp_solver.check_dt(D, v, dt)

    u0 = initial_condition(t0, x, neq)
    
    bvp_solver.set_nr_solver(u0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)

    bc = {'left': bc_left, 'right': bc_right}

    bvp_solver.set_problem(f_res, bc, theta=1.0)

    u = bvp_solver.solve(dtw=0.2)

    u_exact = u_analytical(tf, x)
    err = np.max(np.abs(u_exact-u))
    print(f'err = {err:.4e}')

    for ieq in range(neq):
        ui_num, ui_ex = u[ieq::neq].copy(), u_exact[ieq::neq].copy()
        plot.plot_ode_system_evolution(x, [ui_num, ui_ex])

if __name__ == '__main__':
    main()
