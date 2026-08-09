
from typing import Callable
import numpy as np
from numpy.typing import NDArray

class FEStabilization():

    def __init__(self, f_quantities: Callable, f_strong: Callable,
        f_operator: Callable, order: int = 1):

        self.f_quantities = f_quantities
        self.f_strong = f_strong
        self.f_operator = f_operator
        self.order = order

    def tau_stabilization(self, t: float, x: NDArray[np.float64],
        u: NDArray[np.float64], h_el: float, dt: float) -> NDArray[np.float64]:

        neq = u.shape[0]

        p = self.order

        v, D, s = self.f_quantities(t, x, u)

        tau = np.zeros(neq)

        for ieq in range(neq):

            v_i = np.abs(v[ieq])
            v_i = max(1e-12, v_i)

            D_i = np.abs(D[ieq])
            D_i = max(1e-12, D_i)

            Pe_i = v_i*h_el/(2.0*D_i)

            s_i = s[ieq]

            # tau[ieq] = (h_el/(2.0*v_i))*(1.0/np.tanh(Pe_i) - 1.0/Pe_i) # linear elements: convection-diffusion
            # tau[ieq] = (h_el/(2.0*v_i))*((1.0 + 1.0/Pe_i + h_el*s_i/(2.0*v_i))**(-1)) # Codina (2nd order)
            # tau[ieq] = (h_el/(2.0*v_i))*((1.0 + 9.0/(Pe_i**2) + (h_el*s_i/(2.0*v_i))**2)**(-0.5)) # Shakib (4th order)
            tau[ieq] = ((2.0/dt)**2 + (2.0*v_i*p/h_el)**2 + (6.0*D_i*(p**2)/(h_el**2))**2 + s_i**2)**(-0.5)

        return tau

    def stabilization_terms(self, t: float, x: NDArray[np.float64],
        u: NDArray[np.float64], r_strong: NDArray[np.float64],
        w: float, dwdt: float, grad_w: NDArray[np.float64],
        hess_w: NDArray[np.float64], tau: NDArray[np.float64])-> NDArray[np.float64]:

        p_mat = self.f_operator(t, x, u, w, dwdt, grad_w, hess_w)

        terms = np.dot(p_mat, tau*r_strong)
        # neq = u.shape[0]
        # terms = np.zeros(neq)
        # for ieq in range(neq):
        #     terms[ieq] = np.dot(p_mat[ieq,:], tau*r_strong)

        return terms