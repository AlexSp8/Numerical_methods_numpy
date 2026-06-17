
from typing import Callable
import numpy as np

class FEStabilization():

    def __init__(self, f_quantities: Callable, f_strong: Callable,
        f_operator: Callable, order: int = 1):

        self.f_quantities = f_quantities
        self.f_strong = f_strong
        self.f_operator = f_operator
        self.order = order

    def tau_stabilization(self, t: float, x: float,
        u: np.ndarray[tuple[int]], h: float, dt: float):

        neq = u.shape[0]

        p = self.order

        v, D, s = self.f_quantities(t, x, u)

        tau = np.zeros(neq)

        for ieq in range(neq):

            v_i = np.abs(v[ieq])
            v_i = max(1e-12, v_i)

            D_i = np.abs(D[ieq])
            D_i = max(1e-12, D_i)

            Pe_i = v_i*h/(2.0*D_i)

            s_i = s[ieq]

            tau[ieq] = (h/(2.0*v_i))*(1.0/np.tanh(Pe_i) - 1.0/Pe_i) # linear elements: convection-diffusion
            # tau[ieq] = (h/(2.0*v_i))*((1.0 + 1.0/Pe_i + h*s_i/(2.0*v_i))**(-1)) # Codina (2nd order)
            # tau[ieq] = (h/(2.0*v_i))*((1.0 + 9.0/(Pe_i**2) + (h*s_i/(2.0*v_i))**2)**(-0.5)) # Shakib (4th order)
            tau[ieq] = ((2.0/dt)**2 + (2.0*v_i*p/h)**2 + (6.0*D_i*(p**2)/(h**2))**2 + s_i**2)**(-0.5)

        return tau
    
    def stabilization_terms(self, t: float, x: float,
        u: np.ndarray[tuple[int]], r_strong: np.ndarray[tuple[int]],
        w: float, dwdt: float, dwdx: float, d2wdx2: float,
        tau: np.ndarray[tuple[int]]):

        neq = u.shape[0]

        p_mat = self.f_operator(t, x, u, w, dwdt, dwdx, d2wdx2)
        
        terms = np.zeros(neq)
        for ieq in range(neq):
            s = 0.0
            for jeq in range(neq):
                s += p_mat[ieq,jeq]*tau[jeq]*r_strong[jeq]
            
            terms[ieq] = s

        return terms