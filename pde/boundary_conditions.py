
from typing import Callable
import numpy as np

class BoundaryConditions():

    def __init__(self, neq: int, bc: list[Callable]):

        self.neq = neq
        self.bc = bc

        n_bnd = len(bc)
        self.n_bnd = n_bnd

        self.is_Dirichlet = np.zeros((n_bnd, neq), dtype=bool)

        for i in range(n_bnd):
            self.is_Dirichlet[i,:] = self.classify_bc(bc[i])
        
    def classify_bc(self, bc: Callable):
        """Checks which boundary conditions involve derivatives of the unknowns
        and which are of Dirichlet type.

        Args:
            self: the BoundaryConditions class object
            bc: the boundary condition function
        Returns:
            A boolean array for the BCs on the boundary.
            True indicates a Dirichlet type BC and False indicates a derivative type BC."""

        neq = self.neq
        ndim = 3

        t0_dum = 1.0
        x_b_dum = np.ones(ndim)
        u_b_dum = np.ones(neq)
        dudx_b_dum = np.ones((ndim,neq))

        res = bc(t0_dum, x_b_dum, u_b_dum, dudx_b_dum)

        is_Dirichlet = np.zeros(neq, dtype=bool)

        res_p = bc(t0_dum, x_b_dum, u_b_dum, dudx_b_dum+1e-8)
        for i in range(neq):
            if np.abs(res_p[i] - res[i]) < 1e-12:
                is_Dirichlet[i] = True

        return is_Dirichlet
