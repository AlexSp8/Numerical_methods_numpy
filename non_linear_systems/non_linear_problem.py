
from typing import Callable
from abc import ABC, abstractmethod
import numpy as np
from numpy.typing import NDArray

from ode.time_integration import TimeIntegration
from integration.fe_gauss_integration import FEGaussIntegration
from pde.fe_stabilization import FEStabilization
from pde.mesh_discretization import MeshDiscretization
from pde.boundary_conditions import BoundaryConditions

class NonlinearProblem(ABC):

    def __init__(self, f: Callable[[NDArray[np.float64]], NDArray[np.float64]]):
        self.f = f

    @abstractmethod
    def f_res(self, u: NDArray[np.float64]) -> NDArray[np.float64]:
        """Returns the residual vector of the non-linear system.

        Args:
            u: the vector of unknowns
        Returns:
            The residual vector at the values of the unknowns, u."""

        pass

    @abstractmethod
    def jacobian(self, u: NDArray[np.float64], res: NDArray[np.float64],
        h: float = 1e-8) -> NDArray[np.float64]:
        """Returns the Jacobian of a system of non-linear algebraic equations
        around the values, u, of the unknowns.

        Args:
            res: the residual vector at the values of the unknowns, u
            f_res: the function that calculates the residual
            u: the unknown vector u
        Returns:
            The Jacobian of the non-linear system around the unknown vector, u."""

        n, m = u.shape[0], res.shape[0]

        jac = np.zeros((m,n))

        for j in range(n):

            u_val = u[j]
            u[j] += h
            res_f = self.f_res(u)
            u[j] = u_val

            # u_val = u[j]
            # u[j] -= h
            # res_b = self.f_res(u)
            # u[j] = u_val

            jac[:,j] = (res_f-res)/h
            # jac[:,j] = (res_f-res_b)/(2*h)

        return jac


class Regular(NonlinearProblem):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def f_res(self, u: NDArray[np.float64]) -> NDArray[np.float64]:
        return self.f(u)

    def jacobian(self, u: NDArray[np.float64], res: NDArray[np.float64],
        h: float = 1e-8) -> NDArray[np.float64]:
        return super().jacobian(u, res, h)


class Optimization(NonlinearProblem):

    def __init__(self, *args,
        df: Callable[[Callable[[float], float], float, float, float], float],
        grad_f: Callable[[Callable[[Callable[[float], float], float, float, float], float],
        Callable[[NDArray[np.float64]], float], NDArray[np.float64],
        float], NDArray[np.float64]],
        hessian_f: Callable[[Callable[[Callable[[float], float], float, float, float], float],
        Callable[[NDArray[np.float64]], float], NDArray[np.float64],
        float], NDArray[np.float64]], **kwargs):

        super().__init__(*args, **kwargs)
        self.df = df
        self.grad_f = grad_f
        self.hessian_f = hessian_f

    def f_res(self, x: NDArray[np.float64]
        ) -> NDArray[np.float64]:

        return self.grad_f(self.df, self.f, x, h=1e-6)

    def jacobian(self, x: NDArray[np.float64],
        res: NDArray[np.float64] = None, h: float = 1e-6
        ) -> NDArray[np.float64]:

        # return self.hessian_f(self.df, self.f, x, h=1e-4)
        return super().jacobian(x, res, h)


class LagrangeMultipliers(NonlinearProblem):

    def __init__(self, *args, nx: int,
        g: Callable[[NDArray[np.float64]], NDArray[np.float64]],
        df: Callable[[Callable[[float], float], float, float, float], float],
        grad_f: Callable[[Callable[[Callable[[float], float], float, float, float], float],
        Callable[[NDArray[np.float64]], float], NDArray[np.float64],
        float], NDArray[np.float64]], **kwargs):

        super().__init__(*args, **kwargs)
        self.nx = nx
        self.g = g
        self.df = df
        self.grad_f = grad_f

    def f_res(self, u: NDArray[np.float64]) -> NDArray[np.float64]:

        nx = self.nx

        x = u[:nx]
        l_m = u[nx:]

        r2 = self.g(x)

        grad_f = self.grad_f(self.df, self.f, x)

        jac_g = self.jacobian_g(x, r2)

        r1 =  grad_f + np.dot(l_m, jac_g)

        # return np.concatenate([r1, r2])

        r_tot = np.zeros(nx+l_m.shape[0])
        r_tot[:nx] = r1.copy()
        r_tot[nx:] = r2.copy()

        return r_tot


    def jacobian(self, u: NDArray[np.float64], res: NDArray[np.float64],
        h: float = 1e-8) -> NDArray[np.float64]:

        return super().jacobian(u, res, h)

    def jacobian_g(self, u: NDArray[np.float64], res: NDArray[np.float64],
        h: float = 1e-4) -> NDArray[np.float64]:
        """Returns the Jacobian of the constraints, g, around the unknowns values, u."""

        n, m = u.shape[0], res.shape[0]

        jac = np.zeros((m,n))
        for j in range(n):

            u_val = u[j]
            u[j] += h
            res_f = self.g(u)
            u[j] = u_val

            # u_val = u[j]
            # u[j] -= h
            # res_b = self.g(u)
            # u[j] = u_val

            jac[:,j] = (res_f-res)/h
            # jac[:,j] = (res_f-res_b)/(2*h)

        return jac


class Regression(NonlinearProblem):

    def __init__(self, *args, xi: NDArray[np.float64],
        yi: NDArray[np.float64], **kwargs):

        super().__init__(*args, **kwargs)
        self.xi = xi
        self.yi = yi

    def f_res(self, u: NDArray[np.float64]) -> NDArray[np.float64]:
        return self.yi - self.f(u, self.xi)

    def jacobian(self, u: NDArray[np.float64], res: NDArray[np.float64],
        h: float = 1e-8) -> NDArray[np.float64]:

        return super().jacobian(u, res, h)


class ImplicitODE(NonlinearProblem):

    def __init__(self, *args, yi: NDArray[np.float64], ti: float, h: float,
        method_dict: dict, **kwargs):

        super().__init__(*args, **kwargs)
        self.yi = yi
        self.ti = ti
        self.h = h
        self.method_dict = method_dict

    def update(self, yi: NDArray[np.float64], ti: float, h: float = None):

        self.yi = yi
        self.ti = ti
        if h is not None:
            self.h = h

    def f_res(self, u: NDArray[np.float64]) -> NDArray[np.float64]:

        neq = (self.yi).shape[0]
        ns = self.method_dict['ns']
        q = self.method_dict['q']
        p = self.method_dict['p']

        ti, yi, h = self.ti, self.yi, self.h

        # k_mat = u.reshape(ns, neq)
        k_mat = np.zeros((ns, neq))
        for j in range(ns):
            k_mat[j,:] = u[j*neq:(j+1)*neq]

        res = np.zeros(ns*neq)
        for j in range(ns):
            yj = yi[:] + h*np.dot(q[j,:], k_mat[:,:])
            fj = self.f(ti+p[j]*h, yj)
            res[j*neq:(j+1)*neq] = k_mat[j,:] - fj[:]

        return res

    def jacobian(self, u: NDArray[np.float64], res: NDArray[np.float64],
        h: float = 1e-8) -> NDArray[np.float64]:

        return super().jacobian(u, res, h)


class AdamsMoultonODE(NonlinearProblem):

    def __init__(self, *args, yi: NDArray[np.float64], ti: float, h: float,
        b0: float, **kwargs):

        super().__init__(*args, **kwargs)
        self.yi = yi
        self.ti = ti
        self.h = h
        self.b0 = b0

    def update(self, yi: NDArray[np.float64], ti: float):

        self.yi = yi
        self.ti = ti

    def f_res(self, yn: NDArray[np.float64]) -> NDArray[np.float64]:

        ti, h, yi, b0 = self.ti, self.h, self.yi, self.b0
        f_val = self.f(ti+h, yn)
        return yn - yi - h*b0*f_val

    def jacobian(self, u: NDArray[np.float64], res: NDArray[np.float64],
        h: float = 1e-8) -> NDArray[np.float64]:

        return super().jacobian(u, res, h)


class FiniteDifferencesBVP(NonlinearProblem):

    def __init__(self, *args, neq: int, mesh: MeshDiscretization,
        boundary: BoundaryConditions, time_int: TimeIntegration,
        theta: float = 1.0, **kwargs):

        super().__init__(*args, **kwargs)
        self.neq = neq
        self.mesh = mesh
        self.boundary = boundary
        self.time_int = time_int
        self.theta = theta

    def f_res(self, u: NDArray[np.float64]) -> NDArray[np.float64]:

        f = self.f

        neq = self.neq
        theta = self.theta

        x = self.mesh.x_mesh

        dt = self.time_int.dt
        t = self.time_int.t
        up = self.time_int.up

        ti = theta*t + (1.0-theta)*(t-dt)

        nunknowns = u.shape[0]

        res = np.zeros(nunknowns)

        bulk_nodes = self.mesh.bulk_nodes
        for i in range(len(bulk_nodes)):

            gnod = bulk_nodes[i]

            row1 = gnod*neq
            rowf = row1+neq

            ui = u[row1:rowf]
            ui_p = up[row1:rowf]
            u_val = theta*ui + (1.0-theta)*ui_p

            dudt, grad_u, hess_u = self.calculate_derivatives(u, gnod)

            res[row1:rowf] = f(ti, x[gnod], u_val, dudt, grad_u, hess_u)

        self.set_boundary_conditions(u, res)

        return res

    def calculate_derivatives(self, u: NDArray[np.float64], gnod: int):

        neq = self.neq
        theta = self.theta

        up = self.time_int.up
        dt = self.time_int.dt
        dtp = self.time_int.dtp

        x = self.mesh.x_mesh
        connectivity = self.mesh.connectivity
        ndim = self.mesh.ndim

        row1 = gnod*neq
        rowf = row1+neq

        ui = u[row1:rowf]
        ui_p = up[row1:rowf]
        ui_2p = self.time_int.up2[row1:rowf]
        if self.time_int.inc == 1:
            ui_2p = 2*ui_p - ui

        o1_h2_p_i = (2*dt+dtp)/(dt*(dt+dtp))  # 1st order (h2) backward, ui
        o1_h2_p_p = -(dt+dtp)/(dt*dtp)        # 1st order (h2) backward, up
        o1_h2_p_2p = dt/(dtp*(dt+dtp))        # 1st order (h2) backward, u2p

        # dudt = (ui - ui_p)/dt # O(h) backward
        # dudt = (3.0*ui - 4.0*ui_p + ui_2p)/(2.0*dt) # O(h2) backward, dt: constant
        dudt = (o1_h2_p_i*ui + o1_h2_p_p*ui_p + o1_h2_p_2p*ui_2p) # O(h2) backward

        grad_u = np.zeros((ndim, neq))
        hess_u = np.zeros((ndim, ndim, neq))

        for d in range(ndim):

            idx_b, idx_f = connectivity[gnod,d,0], connectivity[gnod,d,1]

            hf = x[idx_f,d] - x[gnod,d]
            hb = x[gnod,d] - x[idx_b,d]

            ub = u[idx_b*neq : idx_b*neq+neq]
            uf = u[idx_f*neq : idx_f*neq+neq]
            ub_p = up[idx_b*neq : idx_b*neq+neq]
            uf_p = up[idx_f*neq : idx_f*neq+neq]

            # idx_2b = connectivity[idx_b,d,0]
            # idx_2f = connectivity[idx_f,d,1]
            # if idx_f != -1 and idx_2f != -1:
            #     # Forward FD O(h2)
            #     h2f = x[idx_2f,d] - x[idx_f,d]
            #     o1_h2_f_i = -(2*hf+h2f)/(hf*(hf+h2f)) # 1st order (h2) forward, ui
            #     o1_h2_f_f = (hf+h2f)/(hf*h2f)         # 1st order (h2) forward, uf
            #     o1_h2_f_2f = -hf/(h2f*(hf+h2f))       # 1st order (h2) forward, u2f
            #     u2f = u[idx_2f*neq : idx_2f*neq+neq]
            #     u2f_p = up[idx_2f*neq : idx_2f*neq+neq]
            #     dudx = o1_h2_f_i*ui + o1_h2_f_f*uf + o1_h2_f_2f*u2f # O(h2) forward
            #     dudx_p = o1_h2_f_i*ui_p + o1_h2_f_f*uf_p + o1_h2_f_2f*u2f_p # O(h2) forward
            # elif idx_b != -1 and idx_2b != -1:
            #     # Backward FD O(h2)
            #     h2b = x[idx_b,d] - x[idx_2b,d]
            #     o1_h2_b_i = (2*hb+h2b)/(hb*(hb+h2b)) # 1st order (h2) backward, ui
            #     o1_h2_b_b = -(hb+h2b)/(hb*h2b)       # 1st order (h2) backward, ub
            #     o1_h2_b_2b = hb/(h2b*(hb+h2b))       # 1st order (h2) backward, u2b
            #     u2b = u[idx_2b*neq : idx_2b*neq + neq]
            #     u2b_p = up[idx_2b*neq : idx_2b*neq + neq]
            #     dudx = o1_h2_b_i*ui + o1_h2_b_b*ub + o1_h2_b_2b*u2b # O(h2) backward
            #     dudx_p = o1_h2_b_i*ui_p + o1_h2_b_b*ub_p + o1_h2_b_2b*u2b_p # O(h2) backward
            # elif idx_f != -1:
            #     # Forward FD O(h)
            #     dudx = (uf - ui)/hf # O(h) forward
            #     dudx_p = (uf_p - ui_p)/hf # O(h) forward
            # elif idx_b != -1:
            #     # Backward FD O(h)
            #     dudx = (ui - ub)/hb # O(h) backward
            #     dudx_p = (ui_p - ub_p)/hb # O(h) backward

            # Central FD
            o1_h2_c = 1.0/(hf*hb*(hf+hb)) # 1st order (h2) central
            # dudx = (uf - ub)/(hf+hb) # O(h) central
            # dudx_p = (uf_p - ub_p)/(hf+hb) # O(h) central
            dudx = ((hb**2)*uf + (hf**2-hb**2)*ui - (hf**2)*ub)*o1_h2_c # O(h2) central
            dudx_p = ((hb**2)*uf_p + (hf**2-hb**2)*ui_p - (hf**2)*ub_p)*o1_h2_c # O(h2) central

            o2_h2_c = 2.0/(hf+hb) # 2nd order (h2) central
            d2udx2 = ((uf-ui)/hf - (ui-ub)/hb)*o2_h2_c # O(h2) central
            d2udx2_p = ((uf_p-ui_p)/hf - (ui_p-ub_p)/hb)*o2_h2_c # O(h2) central

            grad_u[d,:] = theta*dudx + (1.0-theta)*dudx_p
            hess_u[d,d,:] = theta*d2udx2 + (1.0-theta)*d2udx2_p

            for dj in range(d+1,ndim):

                idx_b_j, idx_f_j = connectivity[gnod,dj,0], connectivity[gnod,dj,1]

                hf_j = x[idx_f_j,dj] - x[gnod,dj]
                hb_j = x[gnod,dj] - x[idx_b_j,dj]

                idx_ff = connectivity[idx_f, dj, 1]
                uff = u[idx_ff*neq : idx_ff*neq + neq]
                uff_p = up[idx_ff*neq : idx_ff*neq + neq]

                idx_fb = connectivity[idx_f, dj, 0]
                ufb = u[idx_fb*neq : idx_fb*neq + neq]
                ufb_p = up[idx_fb*neq : idx_fb*neq + neq]

                idx_bf = connectivity[idx_b, dj, 1]
                ubf = u[idx_bf*neq : idx_bf*neq + neq]
                ubf_p = up[idx_bf*neq : idx_bf*neq + neq]

                idx_bb = connectivity[idx_b, dj, 0]
                ubb = u[idx_bb*neq : idx_bb*neq + neq]
                ubb_p = up[idx_bb*neq : idx_bb*neq + neq]

                denom = (hf+hb)*(hf_j+hb_j)

                d2u_mixed = (uff-ufb-ubf+ubb)/denom
                d2u_mixed_p = (uff_p-ufb_p-ubf_p+ubb_p)/denom

                hess_u[d,dj,:] = theta*d2u_mixed + (1.0-theta)*d2u_mixed_p
                hess_u[dj,d,:] = hess_u[d,dj,:]

        return dudt, grad_u, hess_u

    def set_boundary_conditions(self, u: NDArray[np.float64], res: NDArray[np.float64]):

        neq = self.neq

        boundary_nodes = self.mesh.boundary_nodes

        f_bc = self.boundary.f_bc
        is_Dirichlet = self.boundary.is_Dirichlet

        for i in range(len(boundary_nodes)):

            gnod, bnds = boundary_nodes[i]

            row1 = gnod*neq
            rowf = row1+neq

            node_res = np.zeros(neq)

            for ieq in range(neq):

                dirichlet_res = []
                flux_res = []

                for j in range(len(bnds)):

                    bnd = bnds[j]
                    f = f_bc[bnd]
                    res_bc = self.bc_residual(u, gnod, f)

                    if is_Dirichlet[bnd, ieq]:
                        dirichlet_res.append(res_bc[ieq])
                    else:
                        flux_res.append(res_bc[ieq])

                if len(dirichlet_res) > 0:
                    node_res[ieq] = sum(dirichlet_res)/len(dirichlet_res)
                else:
                    node_res[ieq] = sum(flux_res)/len(flux_res)
                    # node_res[ieq] = flux_res[0]

            res[row1:rowf] = node_res

    def bc_residual(self, u: NDArray[np.float64],
        gnod: int, bc: Callable) -> NDArray[np.float64]:
        """Returns the residual of a boundary condition along a boundary.

        Args:
            u: the vector of unknowns
            gnod: the node of the boundary
            bc: the function of the boundary condition residual
        Returns:
            The residual of a boundary condition along a boundary."""

        neq = self.neq
        ndim = self.mesh.ndim
        x = self.mesh.x_mesh
        t = self.time_int.t

        connectivity = self.mesh.connectivity

        row1, rowf = gnod*neq, gnod*neq+neq

        ui = u[row1:rowf]

        grad_u = np.zeros((ndim, neq))

        for d in range(ndim):

            idx_b = connectivity[gnod,d,0]
            idx_f = connectivity[gnod,d,1]

            hf = x[idx_f,d] - x[gnod,d]
            hb = x[gnod,d] - x[idx_b,d]

            idx_2b = connectivity[idx_b,d,0]
            idx_2f = connectivity[idx_f,d,1]

            ub = u[idx_b*neq : idx_b*neq+neq]
            uf = u[idx_f*neq : idx_f*neq+neq]

            if idx_f != -1 and idx_2f != -1:
                # Forward FD O(h2)
                h2f = x[idx_2f,d] - x[idx_f,d]
                o1_h2_f_i = -(2*hf+h2f)/(hf*(hf+h2f)) # 1st order (h2) forward, ui
                o1_h2_f_f = (hf+h2f)/(hf*h2f)         # 1st order (h2) forward, uf
                o1_h2_f_2f = -hf/(h2f*(hf+h2f))       # 1st order (h2) forward, u2f
                u2f = u[idx_2f*neq : idx_2f*neq+neq]
                grad_u[d,:] = o1_h2_f_i*ui + o1_h2_f_f*uf + o1_h2_f_2f*u2f # O(h2) forward
            elif idx_b != -1 and idx_2b != -1:
                # Backward FD O(h2)
                h2b = x[idx_b,d] - x[idx_2b,d]
                o1_h2_b_i = (2*hb+h2b)/(hb*(hb+h2b)) # 1st order (h2) backward, ui
                o1_h2_b_b = -(hb+h2b)/(hb*h2b)       # 1st order (h2) backward, ub
                o1_h2_b_2b = hb/(h2b*(hb+h2b))       # 1st order (h2) backward, u2b
                u2b = u[idx_2b*neq : idx_2b*neq + neq]
                grad_u[d,:] = o1_h2_b_i*ui + o1_h2_b_b*ub + o1_h2_b_2b*u2b # O(h2) backward
            elif idx_f != -1:
                # Forward FD O(h)
                grad_u[d,:] = (uf - ui)/hf # O(h) forward
            elif idx_b != -1:
                # Backward FD O(h)
                grad_u[d,:] = (ui - ub)/hb # O(h) backward

        return bc(t, x[gnod], ui, grad_u)

    def jacobian(self, u: NDArray[np.float64], res: NDArray[np.float64],
        h: float = 1e-8) -> NDArray[np.float64]:

        return super().jacobian(u, res, h)


class FiniteElementsBVP(NonlinearProblem):

    def __init__(self, *args, neq: int, mesh: MeshDiscretization, boundary: BoundaryConditions,
        time_int: TimeIntegration, gauss_int: FEGaussIntegration, fe_stab: FEStabilization,
        theta: float = 1.0, **kwargs):

        super().__init__(*args, **kwargs)

        self.neq = neq
        self.mesh = mesh
        self.boundary = boundary
        self.time_int = time_int
        self.gauss_int = gauss_int
        self.fe_stab = fe_stab
        self.theta = theta

    def f_res(self, u: NDArray[np.float64]) -> NDArray[np.float64]:

        f_weak = self.f

        neq = self.neq

        ndim = self.mesh.ndim
        nbf = self.mesh.nbf
        nel = self.mesh.nel
        connectivity = self.mesh.connectivity

        ng = self.gauss_int.ng
        w_g = self.gauss_int.w
        bfn = self.gauss_int.bfn
        grad_fc = self.gauss_int.grad_fc
        hess_fc = self.gauss_int.hess_fc

        fe_stab = self.fe_stab

        theta = self.theta

        t = self.time_int.t

        dt = self.time_int.dt
        dtp = self.time_int.dtp

        o1_h2_p_i = (2*dt+dtp)/(dt*(dt+dtp))  # 1st order (h2) backward, ui
        o1_h2_p_p = -(dt+dtp)/(dt*dtp)        # 1st order (h2) backward, up
        o1_h2_p_2p = dt/(dtp*(dt+dtp))        # 1st order (h2) backward, u2p

        t_val = theta*t + (1.0-theta)*(t-dt)

        nunknowns = u.shape[0]

        res = np.zeros(nunknowns)

        for iel in range(nel):

            gnodes = connectivity[iel,:]

            x_el, u_el, up_el, u2p_el = self.elemental_arrays(u, gnodes)

            for ig in range(ng):

                bfn_ig = bfn[:,ig]
                grad_fc_ig = grad_fc[:,:,ig]
                hess_fc_ig = hess_fc[:,:,:,ig]

                det_jac_ig, grad_fx_ig, hess_fx_ig = self.det_jac(x_el, grad_fc_ig, hess_fc_ig)

                w_ig = w_g[ig]*det_jac_ig

                x_ig = np.dot(bfn_ig, x_el)

                u_ig = np.dot(bfn_ig, u_el)
                up_ig = np.dot(bfn_ig, up_el)
                u2p_ig = np.dot(bfn_ig, u2p_el)

                grad_u_ig = np.dot(grad_fx_ig, u_el)
                grad_up_ig = np.dot(grad_fx_ig, up_el)

                hess_u_ig = np.dot(hess_fx_ig, u_el)
                hess_up_ig = np.dot(hess_fx_ig, up_el)

                # dudt_ig = (u_ig-up_ig)/dt # O(h) backward
                # dudt_ig = (3.0*u_ig - 4.0*up_ig + u2p_ig)/(2.0*dt) # O(h2) backward, dt: constant
                dudt_ig = (o1_h2_p_i*u_ig + o1_h2_p_p*up_ig + o1_h2_p_2p*u2p_ig) # O(h2) backward

                u_val = theta*u_ig + (1.0-theta)*up_ig
                grad_u_val = theta*grad_u_ig + (1.0-theta)*grad_up_ig
                hess_u_val = theta*hess_u_ig + (1.0-theta)*hess_up_ig

                h_el = det_jac_ig**(1.0/ndim)
                h_el *= 2

                tau_ig = fe_stab.tau_stabilization(t_val, x_ig, u_ig, h_el, dt)

                r_strong = fe_stab.f_strong(t_val, x_ig, u_val, dudt_ig, grad_u_val, hess_u_val)

                # res_ig = np.zeros(nbf*neq)
                for inod in range(nbf):

                    gnod = gnodes[inod]
                    row1 = neq*gnod
                    rowf = row1+neq

                    w = bfn_ig[inod]
                    grad_w = grad_fx_ig[:,inod]
                    hess_w = hess_fx_ig[:,:,inod]
                    dwdt = 0.0

                    res_ig = f_weak(t_val, x_ig, u_val, dudt_ig, grad_u_val, w, grad_w)

                    res_ig += fe_stab.stabilization_terms(t_val, x_ig, u_val, r_strong,
                                                        w, dwdt, grad_w, hess_w, tau_ig)

                    res[row1:rowf] += res_ig*w_ig

        self.set_boundary_conditions(u, res)

        return res

    def elemental_arrays(self, u: NDArray[np.float64], gnodes: NDArray[np.float64]):
        """Returns the arrays with the values of the mesh and the unknowns at the
        local element.

        Args:
            self: the FiniteElementsBVP class object
            u: the global array of unknowns
            gnodes: the global enumeration of the nodes of the element
        Returns:
            The arrays of values (x, u, up, u2p) at the element nodes."""

        neq = self.neq

        up = self.time_int.up
        u2p = self.time_int.up2

        nbf = self.gauss_int.nbf

        x_el = self.mesh.x_mesh[gnodes,:]

        u_el = np.zeros((nbf, neq))
        up_el = np.zeros((nbf, neq))
        u2p_el = np.zeros((nbf, neq))
        for inod in range(nbf):
            gnod = gnodes[inod]
            row1 = gnod*neq
            rowf = row1+neq
            u_el[inod,:] = u[row1:rowf]
            up_el[inod,:] = up[row1:rowf]
            u2p_el[inod,:] = u2p[row1:rowf]

        return x_el, u_el, up_el, u2p_el

    def det_jac(self, x_el: NDArray[np.float64], grad_fc: NDArray[np.float64],
        hess_fc: NDArray[np.float64]) -> float:
        """Returns the determinant of the transformation jacobian
        from the parent to the physical element, and the
        derivatives with respect to the physical at the Gauss point.

        Args:
            self: the FiniteElementsBVP class object
            x_el: the array of mesh values at the element nodes
            grad_fc: the gradient of the basis functions with
                respect to the parent coordinates at the Gauss point (ndim, nbf).
            hess_fc: the hessian of the basis functions with
                respect to the parent coordinates at the Gauss point (ndim, ndim, nbf).
        Returns:
            The determinant of the transformation jacobian and the
            physical derivatives at the Gauss point."""

        jac = np.dot(grad_fc, x_el)

        det_jac = np.abs(np.linalg.det(jac))

        inv_jac = np.linalg.inv(jac)

        inv_jacT = np.transpose(inv_jac)

        grad_fx = np.dot(inv_jacT, grad_fc)

        nbf, ndim = x_el.shape

        hess_xc = np.dot(hess_fc, x_el)

        hess_fx = np.zeros((ndim, ndim, nbf))
        for inod in range(nbf):
            hess_fc_mod = hess_fc[:, :, inod] - np.dot(hess_xc, grad_fx[:,inod])
            # J^-1 @ H @ J^-T
            hess_fx[:, :, inod] = np.matmul(np.matmul(inv_jac, hess_fc_mod), inv_jacT)

        return det_jac, grad_fx, hess_fx

    def set_boundary_conditions(self, u: NDArray[np.float64], res: NDArray[np.float64]):

        neq = self.neq

        ndim = self.mesh.ndim

        boundary_nodes = self.mesh.boundary_nodes

        f_bc = self.boundary.f_bc
        is_Dirichlet = self.boundary.is_Dirichlet

        t = self.time_int.t

        grad_u = np.zeros((ndim, neq))

        for i in range(len(boundary_nodes)):

            gnod, bnds = boundary_nodes[i]

            x = self.mesh.x_mesh[gnod,:]

            row1 = gnod*neq
            rowf = row1+neq

            ui = u[row1:rowf]

            for ieq in range(neq):

                dirichlet_res = []
                flux_res = []

                for j in range(len(bnds)):

                    bnd = bnds[j]
                    if ndim == 1:
                        grad_u = self.calculate_boundary_derivatives_1d(bnd, u)
                    f = f_bc[bnd]
                    res_bc = f(t, x, ui, grad_u)

                    if is_Dirichlet[bnd, ieq]:
                        dirichlet_res.append(res_bc[ieq])
                    else:
                        flux_res.append(res_bc[ieq])

                irow = row1+ieq
                if len(dirichlet_res) > 0:
                    res[irow] = sum(dirichlet_res)/len(dirichlet_res)
                else:
                    res[irow] += sum(flux_res)/len(flux_res)
                    # res[irow] = flux_res[0]

    def calculate_boundary_derivatives_1d(self, bnd: int, u: NDArray[np.float64]):

        neq = self.neq

        ndim = self.mesh.ndim
        fe_order = self.mesh.fe_order
        x = self.mesh.x_mesh
        nnodes = x.shape[0]

        if bnd == 0:

            gnod = 0
            row1 = gnod*neq

            r = row1
            ui = u[r:r+neq]
            r += neq+1
            uf = u[r:r+neq]
            r += neq+1
            u2f = u[r:r+neq]
            # forward fd O(h2)
            hf = x[1,0]-x[0,0]
            h2f = x[2,0]-x[1,0]

            o1_h2_f_i = -(2*hf+h2f)/(hf*(hf+h2f)) # 1st order (h2) forward, ui
            o1_h2_f_f = (hf+h2f)/(hf*h2f)         # 1st order (h2) forward, uf
            o1_h2_f_2f = -hf/(h2f*(hf+h2f))       # 1st order (h2) forward, u2f
        elif bnd == 1:

            gnod = nnodes-1
            row1 = gnod*neq

            r = row1
            ui = u[r:r+neq]
            r -= neq
            ub = u[r:r+neq]
            r -= neq
            u2b = u[r:r+neq]

            # backward fd O(h2)
            hb = x[gnod,0]-x[gnod-1,0]
            h2b = x[gnod-1,0]-x[gnod-2,0]

            o1_h2_b_i = (2*hb+h2b)/(hb*(hb+h2b)) # 1st order (h2) backward, ui
            o1_h2_b_b = -(hb+h2b)/(hb*h2b)       # 1st order (h2) backward, ub
            o1_h2_b_2b = hb/(h2b*(hb+h2b))       # 1st order (h2) backward, u2b
        else:
            raise ValueError(f'Wrong bnd = {bnd} in calculate_boundary_derivatives_1d!')

        grad_u = np.zeros((ndim,neq))

        if fe_order == 'linear':

            if bnd == 0:
                # x0
                # dfdx = np.array([-1.0/hf, +1.0/hf])
                # grad_u[0,:] = dfdx[0]*ui + dfdx[1]*uf
                grad_u[0,:] = o1_h2_f_i*ui + o1_h2_f_f*uf + o1_h2_f_2f*u2f
            elif bnd == 1:
                # xf
                # dfdx = np.array([-1.0/hb, +1.0/hb])
                # grad_u[0,:] = dfdx[0]*ub + dfdx[1]*ui
                grad_u[0,:] = o1_h2_b_i*ui + o1_h2_b_b*ub + o1_h2_b_2b*u2b

        elif fe_order == 'quadratic':

            if bnd == 0:
                # x0
                x0, x1, x2 = x[0,0], x[1,0], x[2,0]
                dfdx = np.array([(2*x0 - x1 - x2)/((x0 - x1)*(x0 - x2)),
                                (x0 - x2)/((x1 - x0)*(x1 - x2)),
                                (x0 - x1)/((x2 - x0)*(x2 - x1))])
                grad_u[0,:] = dfdx[0]*ui + dfdx[1]*uf + dfdx[2]*u2f
            elif bnd == 1:
                # xf
                x0, x1, x2 = x[gnod-2,0], x[gnod-1,0], x[gnod,0]  # x2 is the boundary coordinate xf
                dfdx = np.array([(x2 - x1)/((x0 - x1)*(x0 - x2)),
                                (x2 - x0)/((x1 - x0)*(x1 - x2)),
                                (2*x2 - x0 - x1)/((x2 - x0)*(x2 - x1))])
                grad_u[0,:] = dfdx[0]*u2b + dfdx[1]*ub + dfdx[2]*ui

        return grad_u

    def jacobian(self, u: NDArray[np.float64], res: NDArray[np.float64],
        h: float = 1e-8) -> NDArray[np.float64]:

        return super().jacobian(u, res, h)
