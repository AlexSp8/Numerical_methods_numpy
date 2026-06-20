
import sys
from pathlib import Path

import numpy as np

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from linear_systems import direct_solver
from ode import plot
from pde import bvp_setup, mesh_discretization, boundary_conditions

# Heat Transfer
def bc_x0(t: float, x: float, u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0] - 40.0
    return res

def bc_xf(t: float, x: float, u: np.ndarray[tuple[int]],
    grad_u: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    res[0] = u[0]- 200.0
    # res[0] = grad_u[0,0]- 0.0
    return res

def initial_condition(t0: float, x: np.ndarray[tuple[int, int]], neq: int
    ) -> np.ndarray[tuple[int]]:

    nnodes = x.shape[0]

    u0 = np.zeros(nnodes*neq)

    u_left = np.zeros(neq)
    u_right = np.zeros(neq)

    ieq = 0
    u0[ieq::neq] = 0.0 # x

    u_left[ieq] = 40.0
    u_right[ieq] = 200.0

    u0[:neq] = u_left
    u0[-neq:] = u_right

    return u0

def f_res(t: float, x: np.ndarray[tuple[int]], u: np.ndarray[tuple[int]],
    dudt: np.ndarray[tuple[int]], grad_u: np.ndarray[tuple[int, int]],
    hess_u: np.ndarray[tuple[int, int, int]]) -> np.ndarray[tuple[int]]:
    neq = u.shape[0]
    res = np.zeros(neq)
    d2udx2 = hess_u[0,0,0]
    h, Ta = 0.01, 20.0
    res[0] = (d2udx2 + h*(Ta-u[0]))
    return res

def u_analytical(t: float, x: np.ndarray[tuple[int, int]]) -> np.ndarray[tuple[int]]:

    xm = x[:,0]

    Ta, L = 20.0, 10.0
    c = (180-20*np.cosh(0.1*L))/np.sinh(0.1*L)
    T_ss = Ta*np.cosh(0.1*xm) + c*np.sinh(0.1*xm) + Ta
    n_tot = 5
    L2 = L**2
    T_tot = T_ss
    for n in range(1, n_tot):
        n_p = np.pi*n
        c1 = 2.0/(n_p*(0.01*L2+n_p**2))
        c2 = ((-1)**n)*(0.2*L2+200*(n_p)**2) - (0.2*L2+40*n_p)
        an = c1*c2
        T_tot += an*np.sin(n_p*xm/L)*np.exp(-(n_p**2+1)*t/L2)
    return T_tot

def main():

    neq = 1

    x0 = np.array([0.0])
    xf = np.array([10.0])
    nodes_dim = np.array([11])
    
    mesh = mesh_discretization.MeshDiscretization(x0, xf, nodes_dim)

    p = np.array([1])
    mesh.set_rectangular_mesh(p, major_order='row')

    bc = [bc_x0, bc_xf]
    boundary = boundary_conditions.BoundaryConditions(neq, bc)

    bvp_solver = bvp_setup.BVPSetupFD(neq, mesh, boundary, f_res)
    
    ls_solver = direct_solver.LUSolver()
    bvp_solver.set_ls_solver(ls_solver)

    t0, tf, dt, dt_min, dt_max, atol, rtol = 0.0, 250.0, 1e-1, 1e-1, 1.0, 1e-3, 1e-2
    dt = tf
    
    u0 = initial_condition(t0, mesh.x_mesh, neq)

    bvp_solver.set_nr_solver(u0, k_max=100, tol=1e-8, r=1.0)

    bvp_solver.set_time_int(t0, tf, dt, dt_min, dt_max, atol, rtol, u0)

    bvp_solver.set_fd_problem(theta=1.0)

    u = bvp_solver.solve(dtw=tf)

    x = mesh.x_mesh
    
    u_exact = u_analytical(tf, x)
    err = np.max(np.abs(u_exact-u))
    print(f'err = {err:.4e}')

    for ieq in range(neq):
        ui_num, ui_ex = u[ieq::neq].copy(), u_exact[ieq::neq].copy()
        plot.plot_ode_system_evolution(x, [ui_num, ui_ex])

if __name__ == '__main__':
    main()
