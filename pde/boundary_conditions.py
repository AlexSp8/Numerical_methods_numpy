
from typing import Callable
import numpy as np


class BoundaryConditions():

    def __init__(self, f_bc: list[Callable], bc_type: list[list[str]]):

        self.f_bc = f_bc

        self.n_bnd = len(f_bc)
        self.neq = len(bc_type[1])

        self.is_Dirichlet = np.zeros((self.n_bnd, self.neq), dtype=bool)

        for ibd in range(self.n_bnd):
            bc_type_bnd = bc_type[ibd]
            for ieq in range(self.neq):
                self.is_Dirichlet[ibd,ieq] = (bc_type_bnd[ieq] == 'dirichlet')

        # for i in range(n_bnd):
        #     self.is_Dirichlet[i,:] = self.classify_bc(f_bc[i])

    def classify_bc(self, f_bc: Callable):
        """Checks which boundary conditions involve derivatives of the unknowns
        and which are of Dirichlet type.

        Args:
            self: the BoundaryConditions class object
            f_bc: the boundary condition function
        Returns:
            A boolean array for the BCs on the boundary.
            True indicates a Dirichlet type BC and False indicates a derivative type BC."""

        neq = self.neq
        ndim = 3

        t0_dum = 1.0
        x_b_dum = np.ones(ndim)
        u_b_dum = np.ones(neq)
        grad_u_dum = np.ones((ndim,neq))

        res = f_bc(t0_dum, x_b_dum, u_b_dum, grad_u_dum)

        is_Dirichlet = np.zeros(neq, dtype=bool)

        res_p = f_bc(t0_dum, x_b_dum, u_b_dum, grad_u_dum+1e-8)
        for i in range(neq):
            if np.abs(res_p[i] - res[i]) < 1e-12:
                is_Dirichlet[i] = True

        return is_Dirichlet
