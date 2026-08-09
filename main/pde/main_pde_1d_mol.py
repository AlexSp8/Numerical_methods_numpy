
import numpy as np

from ode import bvp_1d_fd, plot, ivp
from pde import mol_1d

def bc_left_main(t: float, x_b: float, u_b: NDArray[np.float64],
    dudx_b: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u_b.shape[0]
    res = np.zeros(neq)
    res[0] = u_b[0] - 40.0
    return res

def bc_right_main(t: float, x_b: float, u_b: NDArray[np.float64],
    dudx_b: NDArray[np.float64]) -> NDArray[np.float64]:
    neq = u_b.shape[0]
    res = np.zeros(neq)
    res[0] = u_b[0]- 200.0
    # res[0] = dudx_b[0]- 0.0
    return res

def initial_condition_main(t0: float, x: NDArray[np.float64], neq: int
    ) -> NDArray[np.float64]:

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

def f_main(t: float, x: float, u: NDArray[np.float64], dudx: NDArray[np.float64],
    d2udx2: NDArray[np.float64]) -> NDArray[np.float64]:
    T = u[0]
    h, Ta = 0.01, 20.0
    return d2udx2 + h*(Ta-T)

def u_exact_main(t: float, x: NDArray[np.float64]) -> NDArray[np.float64]:

    Ta, L = 20.0, 10.0
    c = (180-20*np.cosh(0.1*L))/np.sinh(0.1*L)
    T_ss = Ta*np.cosh(0.1*x) + c*np.sinh(0.1*x) + Ta
    n_tot = 5
    L2 = L**2
    T_tot = T_ss
    for n in range(1, n_tot):
        n_p = np.pi*n
        c1 = 2.0/(n_p*(0.01*L2+n_p**2))
        c2 = ((-1)**n)*(0.2*L2+200*(n_p)**2) - (0.2*L2+40*n_p)
        an = c1*c2
        T_tot += an*np.sin(n_p*x/L)*np.exp(-(n_p**2+1)*t/L2)
    return T_tot

def setup_main():

    xd = np.array([0.0, 10.0])
    x = bvp_1d_fd.create_mesh(nnodes=11, xd=xd, p=1.0)

    neq = 1

    bc = {'left': bc_left_main, 'right': bc_right_main}
    problem = mol_1d.FiniteDifferencesMOL1D(f_main, x, neq, bc)

    D, v = 1.0, 0.0

    problem.check_dx(x, D, v)

    u0 = initial_condition_main(t0=0.0, x=x, neq=neq)

    t0, tf, method, dt0 = 0.0, 70.0, 'euler', 0.4
    dt = problem.check_dt(x, D, v, dt0)
    ivp_solver = ivp.ExplicitIVPSolver(problem.g, t0, tf, u0, method, h=dt)


    # t0, tf, method = 0.0, 70.0, 'fehlberg'
    # h_min, h_max, atol, rtol = 1e-6, 1.0, 1e-8, 1e-6
    # ivp_solver = ivp.ExplicitAdaptiveIVPSolver(problem.g, t0, tf, u0, method,
    #     h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)

    # t0, tf, method, dt = 0.0, 70.0, 'euler', 1.0
    # ivp_solver = ivp.ImplicitIVPSolver(problem.g, t0, tf, u0, method, h=dt)

    # t0, tf, method = 0.0, 70.0, 'radau-iia-5-adaptive'
    # h_min, h_max, atol, rtol = 1e-1, 10.0, 1e-4, 1e-2
    # ivp_solver = ivp.ImplicitAdaptiveIVPSolver(problem.g, t0, tf, u0, method,
    #     h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)

    return x, ivp_solver


# Time-dependent BCs
def bc_left_prob1(t: float, x: float, u: NDArray[np.float64],
    dudx: NDArray[np.float64]) -> NDArray[np.float64]:
    # Exact value at x=0 is sin(t)*(0) = 0
    return u - 0.0

def bc_right_prob1(t: float, x: float, u: NDArray[np.float64],
    dudx: NDArray[np.float64]) -> NDArray[np.float64]:
    # Exact value at x=2 is sin(t)*(4 + 2) = 6*sin(t)
    return u - 6.0*np.sin(t)

def initial_condition_time(t0: float, x: NDArray[np.float64], neq: int
    ) -> NDArray[np.float64]:
    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)
    return u0

def f_prob1(t: float, x: float, u: NDArray[np.float64],
    dudx: NDArray[np.float64], d2udx2: NDArray[np.float64]
    ) -> NDArray[np.float64]:
    D, v = 1.0, 2.0
    S = np.cos(t)*(x**2 + x) + 4.0*x*np.sin(t)
    return D*d2udx2 - v*dudx + S

def u_exact_time(t: float, x: NDArray[np.float64]) -> NDArray[np.float64]:
    return np.sin(t)*(x**2+x)

def setup_time():

    xd = np.array([0.0, 2.0])
    nnodes = 31
    x = bvp_1d_fd.create_mesh(nnodes=nnodes, xd=xd, p=1.0)

    neq = 1

    bc = {'left': bc_left_prob1, 'right': bc_right_prob1}
    problem = mol_1d.FiniteDifferencesMOL1D(f_prob1, x, neq, bc)

    D, v = 1.0, 2.0

    problem.check_dx(x, D, v)

    u0 = initial_condition_time(t0=0.0, x=x, neq=neq)

    t0, tf, method, dt0 = 0.0, 1.0, 'euler', 0.1
    dt = problem.check_dt(x, D, v, dt0)
    ivp_solver = ivp.ExplicitIVPSolver(problem.g, t0, tf, u0, method, h=dt)

    # t0, tf, method = 0.0, 2.0, 'fehlberg'
    # h_min, h_max, atol, rtol = 1e-6, 1.0, 1e-8, 1e-6
    # ivp_solver = ivp.ExplicitAdaptiveIVPSolver(problem.g, t0, tf, u0, method,
    #     h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)

    # t0, tf, method, dt = 0.0, 2.0, 'euler', 0.1
    # ivp_solver = ivp.ImplicitIVPSolver(problem.g, t0, tf, u0, method, h=dt)

    # t0, tf, method = 0.0, 70.0, 'radau-iia-5-adaptive'
    # h_min, h_max, atol, rtol = 1e-1, 10.0, 1e-4, 1e-2
    # ivp_solver = ivp.ImplicitAdaptiveIVPSolver(problem.g, t0, tf, u0, method,
    #     h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)

    return x, ivp_solver

# Coupled
def bc_left_prob2(t: float, x: float, u: NDArray[np.float64],
    dudx: NDArray[np.float64]) -> NDArray[np.float64]:
    res0 = u[0] - 0.0                      # Dirichlet
    res1 = dudx[1] - (t + 1.0)             # Neumann
    return np.array([res0, res1])

def bc_right_prob2(t: float, x: float, u: NDArray[np.float64],
    dudx: NDArray[np.float64]) -> NDArray[np.float64]:
    res0 = dudx[0] - 2.0*(u[1] - 1.0)      # Cross-coupled Robin
    # res0 = u[0] - t
    res1 = u[1] - (t + 1.0)                # Dirichlet
    return np.array([res0, res1])

def initial_condition_coupled(t0: float, x: NDArray[np.float64], neq: int
    ) -> NDArray[np.float64]:

    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)

    ieq = 1
    u0[ieq::neq] = x

    return u0

def f_prob2(t: float, x: float, u: NDArray[np.float64],
    dudx: NDArray[np.float64], d2udx2: NDArray[np.float64]
    ) -> NDArray[np.float64]:

    S0 = x**2 - 2.0*t - t*(t + 1.0)*(x**2)
    S1 = x - 2.0*t*(t + 1.0)*(x**2)

    dudt0 = d2udx2[0] + u[0]*dudx[1] + S0
    dudt1 = d2udx2[1] + u[1]*dudx[0] + S1
    return np.array([dudt0, dudt1])

def u_exact_coupled(t: float, x: NDArray[np.float64]) -> NDArray[np.float64]:
    nnodes = x.shape[0]
    neq = 2
    u_exact = np.zeros(nnodes*neq)
    u_exact[0::neq] = t*(x**2)
    u_exact[1::neq] = (t + 1.0)*x
    return u_exact

def setup_coupled():

    xd = np.array([0.0, 1.0])
    nnodes = 11
    x = bvp_1d_fd.create_mesh(nnodes=nnodes, xd=xd, p=1.0)

    neq = 2

    bc = {'left': bc_left_prob2, 'right': bc_right_prob2}
    problem = mol_1d.FiniteDifferencesMOL1D(f_prob2, x, neq, bc)

    u0 = initial_condition_coupled(t0=0.0, x=x, neq=neq)

    D, v = 1.0, np.max(u0)

    problem.check_dx(x, D, v)

    # t0, tf, method, dt0 = 0.0, 5.0, 'euler', 0.1
    # dt = problem.check_dt(x, D, v, dt0)
    # ivp_solver = ivp.ExplicitIVPSolver(problem.g, t0, tf, u0, method, h=dt)

    # t0, tf, method = 0.0, 5.0, 'fehlberg'
    # h_min, h_max, atol, rtol = 1e-6, 0.1, 1e-8, 1e-6
    # ivp_solver = ivp.ExplicitAdaptiveIVPSolver(problem.g, t0, tf, u0, method,
    #     h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)

    t0, tf, method, dt = 0.0, 5.0, 'euler', 0.1
    ivp_solver = ivp.ImplicitIVPSolver(problem.g, t0, tf, u0, method, h=dt)

    # t0, tf, method = 0.0, 70.0, 'radau-iia-5-adaptive'
    # h_min, h_max, atol, rtol = 1e-1, 10.0, 1e-4, 1e-2
    # ivp_solver = ivp.ImplicitAdaptiveIVPSolver(problem.g, t0, tf, u0, method,
    #     h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)

    return x, ivp_solver

# Problem 3
def f_prob3(t, x, u, dudx, d2udx2):
    S0 = np.exp(-t)*(np.cos(x) - np.sin(x))
    S1 = np.exp(-t)*(np.sin(x) - np.cos(x))
    dudt0 = d2udx2[0] + S0
    dudt1 = d2udx2[1] + S1
    return np.array([dudt0, dudt1])

def bc_left_prob3(t, x, u, dudx):
    res0 = u[0] - 0.0     # Dirichlet
    res1 = dudx[1] - 0.0  # Neumann
    return np.array([res0, res1])

def bc_right_prob3(t, x, u, dudx):
    res0 = dudx[0] - 0.0  # Neumann
    res1 = u[1] - 0.0     # Dirichlet
    return np.array([res0, res1])

def initial_condition_3(t0: float, x: NDArray[np.float64], neq: int
    ) -> NDArray[np.float64]:

    nnodes = x.shape[0]
    u0 = np.zeros(nnodes*neq)

    u0[0::neq] = np.sin(x)
    u0[1::neq] = np.cos(x)

    return u0

def u_exact_3(t: float, x: NDArray[np.float64]) -> NDArray[np.float64]:
    nnodes = x.shape[0]
    neq = 2
    u_exact = np.zeros(nnodes*neq)
    u_exact[0::neq] = np.exp(-t)*np.sin(x)
    u_exact[1::neq] = np.exp(-t)*np.cos(x)
    return u_exact

def setup_3():

    xd = np.array([0.0, np.pi/2.0])
    nnodes = 11
    x = bvp_1d_fd.create_mesh(nnodes=nnodes, xd=xd, p=1.0)

    neq = 2

    bc = {'left': bc_left_prob3, 'right': bc_right_prob3}
    problem = mol_1d.FiniteDifferencesMOL1D(f_prob3, x, neq, bc)

    u0 = initial_condition_3(t0=0.0, x=x, neq=neq)

    D, v = 1.0, 0.0

    problem.check_dx(x, D, v)

    t0, tf, method, dt0 = 0.0, 5.0, 'euler', 0.1
    dt = problem.check_dt(x, D, v, dt0)
    ivp_solver = ivp.ExplicitIVPSolver(problem.g, t0, tf, u0, method, h=dt)

    # t0, tf, method = 0.0, 5.0, 'fehlberg'
    # h_min, h_max, atol, rtol = 1e-6, 0.1, 1e-8, 1e-6
    # ivp_solver = ivp.ExplicitAdaptiveIVPSolver(problem.g, t0, tf, u0, method,
    #     h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)

    # t0, tf, method, dt = 0.0, 5.0, 'euler', 0.1
    # ivp_solver = ivp.ImplicitIVPSolver(problem.g, t0, tf, u0, method, h=dt)

    # t0, tf, method = 0.0, 70.0, 'radau-iia-5-adaptive'
    # h_min, h_max, atol, rtol = 1e-1, 10.0, 1e-4, 1e-2
    # ivp_solver = ivp.ImplicitAdaptiveIVPSolver(problem.g, t0, tf, u0, method,
    #     h_min=h_min, h_max=h_max, atol=atol, rtol=rtol)

    return x, ivp_solver

def main():

    # pde = setup_main
    # pde = setup_time
    # pde = setup_coupled
    pde = setup_3

    x, ivp_solver = pde()
    t, u = ivp_solver.solve()

    nnodes, nunknowns = x.shape[0], u.shape[1]
    neq = round(nunknowns/nnodes)

    nt = t.shape[0]

    # print(nnodes, nunknowns, neq, nt, round(nt/5))

    for i in range(0, t.shape[0], round(nt/5)):
        # u_exact = u_exact_main(t[i], x)
        # u_exact = u_exact_time(t[i], x)
        # u_exact = u_exact_coupled(t[i], x)
        u_exact = u_exact_3(t[i], x)
        u_num = u[i,:]
        err = np.max(np.abs(u_exact-u_num))
        print(f't = {t[i]:.2e}, err: {err:.4e}')

        for ieq in range(neq):
            ui_num, ui_ex = u_num[ieq::neq].copy(), u_exact[ieq::neq].copy()
            plot.plot_ode_system_evolution(x, [ui_num, ui_ex])

    # for i in range(0, t.shape[0], 50):
    #     u_str = "[" + ", ".join([f"{val:.2e}" for val in u[i,:]]) + "]"
    #     print(f't = {t[i]:.2e}, u: {u_str}')
    #     plot.plot_ode_evolution(x, u[i,:])
    # plot.plot_ode_evolution(t, u)

if __name__ == '__main__':
    main()
