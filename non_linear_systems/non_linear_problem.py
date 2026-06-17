
from typing import Callable
from abc import ABC, abstractmethod
import numpy as np
import numpy.typing as npt

from ode.time_integration import TimeIntegration
from integration.fe_gauss_integration import FEGaussIntegration
from pde.fe_stabilization import FEStabilization
from pde.mesh_discretization import MeshDiscretization

class NonlinearProblem(ABC):

    def __init__(self,
        f: Callable[[np.ndarray[tuple[int]]], np.ndarray[tuple[int]]]):

        self.f = f

    @abstractmethod
    def f_res(self, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
        """Returns the residual vector of the non-linear system.

        Args:
            u: the vector of unknowns
        Returns:
            The residual vector at the values of the unknowns, u."""

        pass

    @abstractmethod
    def jacobian(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-8) -> np.ndarray[tuple[int, int]]:
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

    def f_res(self, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

        return self.f(u)

    def jacobian(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-8) -> np.ndarray[tuple[int, int]]:

        return super().jacobian(u, res, h)


class Optimization(NonlinearProblem):

    def __init__(self, *args,
        df: Callable[[Callable[[float], float], float, float, float], float],
        grad_f: Callable[[Callable[[Callable[[float], float], float, float, float], float],
        Callable[[npt.NDArray[np.float64]], float], npt.NDArray[np.float64],
        float], npt.NDArray[np.float64]],
        hessian_f: Callable[[Callable[[Callable[[float], float], float, float, float], float],
        Callable[[npt.NDArray[np.float64]], float], npt.NDArray[np.float64],
        float], npt.NDArray[np.float64]], **kwargs):

        super().__init__(*args, **kwargs)
        self.df = df
        self.grad_f = grad_f
        self.hessian_f = hessian_f

    def f_res(self, x: np.ndarray[tuple[int]]
        ) -> np.ndarray[tuple[int]]:

        return self.grad_f(self.df, self.f, x, h=1e-6)

    def jacobian(self, x: np.ndarray[tuple[int]],
        res: np.ndarray[tuple[int]] = None, h: float = 1e-6
        ) -> np.ndarray[tuple[int, int]]:

        # return self.hessian_f(self.df, self.f, x, h=1e-4)
        return super().jacobian(x, res, h)


class LagrangeMultipliers(NonlinearProblem):

    def __init__(self, *args, nx: int,
        g: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
        df: Callable[[Callable[[float], float], float, float, float], float],
        grad_f: Callable[[Callable[[Callable[[float], float], float, float, float], float],
        Callable[[npt.NDArray[np.float64]], float], npt.NDArray[np.float64],
        float], npt.NDArray[np.float64]], **kwargs):

        super().__init__(*args, **kwargs)
        self.nx = nx
        self.g = g
        self.df = df
        self.grad_f = grad_f

    def f_res(self, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

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


    def jacobian(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-8) -> np.ndarray[tuple[int, int]]:

        return super().jacobian(u, res, h)

    def jacobian_g(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-4) -> np.ndarray[tuple[int, int]]:
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

    def __init__(self, *args, xi: np.ndarray[tuple[int, int]],
        yi: np.ndarray[tuple[int]], **kwargs):

        super().__init__(*args, **kwargs)
        self.xi = xi
        self.yi = yi

    def f_res(self, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:
        return self.yi - self.f(u, self.xi)

    def jacobian(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-8) -> np.ndarray[tuple[int, int]]:

        return super().jacobian(u, res, h)


class ImplicitODE(NonlinearProblem):

    def __init__(self, *args, yi: np.ndarray[tuple[int]], ti: float, h: float,
        method_dict: dict, **kwargs):

        super().__init__(*args, **kwargs)
        self.yi = yi
        self.ti = ti
        self.h = h
        self.method_dict = method_dict

    def update(self, yi: np.ndarray[tuple[int]], ti: float, h: float = None):

        self.yi = yi
        self.ti = ti
        if h is not None:
            self.h = h

    def f_res(self, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

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

    def jacobian(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-8) -> np.ndarray[tuple[int, int]]:

        return super().jacobian(u, res, h)


class AdamsMoultonODE(NonlinearProblem):

    def __init__(self, *args, yi: np.ndarray[tuple[int]], ti: float, h: float,
        b0: float, **kwargs):

        super().__init__(*args, **kwargs)
        self.yi = yi
        self.ti = ti
        self.h = h
        self.b0 = b0

    def update(self, yi: np.ndarray[tuple[int]], ti: float):

        self.yi = yi
        self.ti = ti

    def f_res(self, yn: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

        ti, h, yi, b0 = self.ti, self.h, self.yi, self.b0
        f_val = self.f(ti+h, yn)
        return yn - yi - h*b0*f_val

    def jacobian(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-8) -> np.ndarray[tuple[int, int]]:

        return super().jacobian(u, res, h)


class FiniteDifferencesBVP1D(NonlinearProblem):

    def __init__(self, *args, neq: int, x: np.ndarray[tuple[int]], bc: dict[str, Callable],
        time_int: TimeIntegration, theta: float = 1.0, **kwargs):

        super().__init__(*args, **kwargs)
        self.neq = neq
        self.x = x
        self.bc = bc
        self.time_int = time_int
        self.theta = theta

    def f_res(self, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

        f = self.f
        neq = self.neq
        x = self.x
        time_int = self.time_int
        theta = self.theta

        t = time_int.t

        up = time_int.up
        u2p = time_int.up2

        if time_int.inc == 1:
            u2p = 2*up - u

        dt = time_int.dt
        dtp = time_int.dtp

        o1_h2_p_i = (2*dt+dtp)/(dt*(dt+dtp))  # 1st order (h2) backward, ui
        o1_h2_p_p = -(dt+dtp)/(dt*dtp)        # 1st order (h2) backward, up
        o1_h2_p_2p = dt/(dtp*(dt+dtp))        # 1st order (h2) backward, u2p

        nnodes, nunknowns = x.shape[0], u.shape[0]

        res = np.zeros(nunknowns)

        for i in range(1, nnodes-1):

            hf = x[i+1] - x[i]
            hb = x[i] - x[i-1]

            o2_h2_c = 2.0/(hf+hb) # 2nd order (h2) central

            # o1_h2_c = 1.0/(hf*hb*(hf+hb)) # 1st order (h2) central

            # h2b = x[i-1] - x[i-2]
            # o1_h2_b_i = (2*hb+h2b)/(hb*(hb+h2b)) # 1st order (h2) backward, ui
            # o1_h2_b_b = -(hb+h2b)/(hb*h2b)       # 1st order (h2) backward, ub
            # o1_h2_b_2b = hb/(h2b*(hb+h2b))       # 1st order (h2) backward, u2b

            # h2f = x[i+2] - x[i+1]
            # o1_h2_f_i = -(2*hf+h2f)/(hf*(hf+h2f)) # 1st order (h2) forward, ui
            # o1_h2_f_f = (hf+h2f)/(hf*h2f)         # 1st order (h2) forward, uf
            # o1_h2_f_2f = -hf/(h2f*(hf+h2f))       # 1st order (h2) forward, u2f

            row1 = i*neq
            rowf = row1+neq

            ub_i = u[row1-neq:row1]
            u_i = u[row1:rowf]
            uf_i = u[rowf:rowf+neq]

            up_i = up[row1:rowf]
            u2p_i = u2p[row1:rowf]
            
            upb_i = up[row1-neq:row1]
            upf_i = up[rowf:rowf+neq]

            # dudt_i = (u_i - up_i)/dt # O(h) backward
            # dudt_i = (3.0*u_i - 4.0*up_i + u2p_i)/(2.0*dt) # O(h2) backward, dt: constant
            dudt_i = (o1_h2_p_i*u_i + o1_h2_p_p*up_i + o1_h2_p_2p*u2p_i) # O(h2) backward

            dudx_i = (uf_i - u_i)/hf # O(h) forward
            # dudx_i = (uf_i - ub_i)/(hf+hb) # O(h) central
            # dudx_i = ((hb**2)*uf_i + (hf**2-hb**2)*u_i + (hf**2)*ub_i)*o1_h2_c # O(h2) central
            
            dupdx_i = (upf_i - up_i)/hf # O(h) forward
            # dupdx_i = (upf_i - upb_i)/(hf+hb) # O(h) central
            # dupdx_i = ((hb**2)*upf_i + (hf**2-hb**2)*up_i + (hf**2)*upb_i)*o1_h2_c # O(h2) central

            d2udx2_i = ((uf_i-u_i)/hf - (u_i-ub_i)/hb)*o2_h2_c # O(h2) central
            
            d2updx2_i = ((upf_i-up_i)/hf - (up_i-upb_i)/hb)*o2_h2_c # O(h2) central
            
            t_val = theta*t + (1.0-theta)*(t-dt)
            u_val = theta*u_i + (1.0-theta)*up_i
            dudx_val = theta*dudx_i + (1.0-theta)*dupdx_i
            d2udx2_val = theta*d2udx2_i + (1.0-theta)*d2updx2_i

            res[row1:rowf] = f(t_val, x[i], u_val, dudt_i, dudx_val, d2udx2_val) # implicit form

        res[:neq] = self.left_bc_residual(u, t)
        res[-neq:] = self.right_bc_residual(u, t)

        return res

    def left_bc_residual(self, u: np.ndarray[tuple[int]],
        t: float) -> np.ndarray[tuple[int]]:
        """Returns the residual on the left boundary of a 1D domain.

        Args:
            u: the vector of unknowns
        Returns:
            The residual on the left boundary after applying the boundary condition."""

        neq = self.neq
        x = self.x
        bc = self.bc['left']

        hf = x[1]-x[0]
        h2f = x[2] - x[1]
        o1_h2_f_i = -(2*hf+h2f)/(hf*(hf+h2f)) # 1st order (h2) forward, ui
        o1_h2_f_b = (hf+h2f)/(hf*h2f)         # 1st order (h2) forward, uf
        o1_h2_f_2b = -hf/(h2f*(hf+h2f))       # 1st order (h2) forward, u2f

        u_b = u[0:neq]
        uf_b = u[neq:2*neq]
        u2f_b = u[2*neq:3*neq]

        # dudx_b = (uf_b - u_b)/hf
        # dudx_b = (-u2f_b+4.0*uf_b-3.0*u_b)/(2.0*hf)
        dudx_b = o1_h2_f_i*u_b + o1_h2_f_b*uf_b + o1_h2_f_2b*u2f_b

        return bc(t, x[0], u_b, dudx_b)

    def right_bc_residual(self, u: np.ndarray[tuple[int]],
        t: float) -> np.ndarray[tuple[int]]:
        """Returns the residual on the right boundary of a 1D domain.

        Args:
            u: the vector of unknowns
        Returns:
            The residual on the right boundary after applying the boundary condition."""

        neq = self.neq
        x = self.x
        bc = self.bc['right']

        nunknowns = u.shape[0]

        hb = x[-1]-x[-2]
        h2b = x[-2] - x[-3]
        o1_h2_b_i = (2*hb+h2b)/(hb*(hb+h2b)) # 1st order (h2) backward, ui
        o1_h2_b_b = -(hb+h2b)/(hb*h2b)       # 1st order (h2) backward, ub
        o1_h2_b_2b = hb/(h2b*(hb+h2b))       # 1st order (h2) backward, u2b

        u_b = u[nunknowns-neq:]
        ub_b = u[nunknowns-2*neq:nunknowns-neq]
        ub_2b = u[nunknowns-3*neq:nunknowns-2*neq]

        # dudx_b = (u_b - ub_b)/hb
        # dudx_b = (3.0*u_b - 4.0*ub_b + ub_2b)/(2.0*hb)
        dudx_b = o1_h2_b_i*u_b + o1_h2_b_b*ub_b + o1_h2_b_2b*ub_2b

        return bc(t, x[-1], u_b, dudx_b)

    def jacobian(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-8) -> np.ndarray[tuple[int, int]]:

        return super().jacobian(u, res, h)


class FiniteElementsBVP1D(NonlinearProblem):

    def __init__(self, *args, neq: int, x: np.ndarray[tuple[int]], bc: dict[str, Callable],
        time_int: TimeIntegration, gauss_int: FEGaussIntegration, fe_stab: FEStabilization,
        theta: float = 1.0, **kwargs):

        super().__init__(*args, **kwargs)

        self.neq = neq
        self.x = x
        self.bc = bc
        self.time_int = time_int
        self.gauss_int = gauss_int
        self.fe_stab = fe_stab
        self.theta = theta
        
        self.is_Dirichlet_left = self._classify_bc(x[0], bc['left'])
        self.is_Dirichlet_right = self._classify_bc(x[-1], bc['right'])

    def _classify_bc(self, x_b, bc):
        """Checks which boundary conditions involve derivatives of the unknowns
        and which are of Dirichlet type.

        Args:
            self: the FiniteElementsBVP1D class object
        Returns:
            Two boolean arrays for the BCs on the left and right boundaries, respectively.
            True indicates a Dirichlet type BC and False indicates a derivative type BC."""

        neq = self.neq

        t0_dum = 1.0
        u_b_dum = np.ones(neq)
        dudx_b_dum = np.ones(neq)

        res = bc(t0_dum, x_b, u_b_dum, dudx_b_dum)

        is_Dirichlet = np.zeros(neq, dtype=bool)

        res_p = bc(t0_dum, x_b, u_b_dum, dudx_b_dum+1e-8)
        for i in range(neq):
            if np.abs(res_p[i] - res[i]) < 1e-12:
                is_Dirichlet[i] = True

        return is_Dirichlet

    def f_res(self, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

        f = self.f

        neq = self.neq

        x = self.x

        gauss_int = self.gauss_int
        ng = gauss_int.ng
        nbf = gauss_int.nbf
        w_g = gauss_int.w
        bfn = gauss_int.bfn
        dfdc = gauss_int.dfdc
        d2fdc2 = gauss_int.d2fdc2

        fe_stab = self.fe_stab
        
        theta = self.theta

        time_int = self.time_int
        t = time_int.t

        up = time_int.up
        u2p = time_int.up2

        if time_int.inc == 1:
            u2p = 2*up - u

        dt = time_int.dt
        dtp = time_int.dtp

        o1_h2_p_i = (2*dt+dtp)/(dt*(dt+dtp))  # 1st order (h2) backward, ui
        o1_h2_p_p = -(dt+dtp)/(dt*dtp)        # 1st order (h2) backward, up
        o1_h2_p_2p = dt/(dtp*(dt+dtp))        # 1st order (h2) backward, u2p
        
        t_val = theta*t + (1.0-theta)*(t-dt)

        nnodes, nunknowns = x.shape[0], u.shape[0]

        nel = round((nnodes-1)/(nbf-1))
        # print(nnodes, neq, nel, nunknowns)

        res = np.zeros(nunknowns)

        for iel in range(nel):

            gnodes = gnodes = np.array([(nbf-1)*iel+inod for inod in range(nbf)])

            x_el, u_el, up_el, u2p_el = self.elemental_arrays(u, gnodes)

            h_el = np.abs(x_el[-1]-x_el[0])

            for ig in range(ng):

                bfn_ig = bfn[:,ig]
                dfdc_ig = dfdc[:,ig]
                d2fdc2_ig = d2fdc2[:,ig]

                det_jac_ig = self.det_jac(x_el, dfdc_ig)

                dfdx_ig = dfdc_ig/det_jac_ig
                d2fdx2_ig = d2fdc2_ig/(det_jac_ig**2)

                x_ig = np.dot(x_el, bfn_ig)

                u_ig = np.dot(bfn_ig, u_el)
                up_ig = np.dot(bfn_ig, up_el)
                u2p_ig = np.dot(bfn_ig, u2p_el)

                dudx_ig = np.dot(dfdx_ig, u_el)
                dupdx_ig = np.dot(dfdx_ig, up_el)

                d2udx2_ig = np.dot(d2fdx2_ig, u_el)
                d2updx2_ig = np.dot(d2fdx2_ig, up_el)

                # dudt_ig = (u_ig-up_ig)/dt # O(h) backward
                # dudt_ig = (3.0*u_ig - 4.0*up_ig + u2p_ig)/(2.0*dt) # O(h2) backward, dt: constant
                dudt_ig = (o1_h2_p_i*u_ig + o1_h2_p_p*up_ig + o1_h2_p_2p*u2p_ig) # O(h2) backward

                w_ig = w_g[ig]*det_jac_ig

                u_val = theta*u_ig + (1.0-theta)*up_ig
                dudx_val = theta*dudx_ig + (1.0-theta)*dupdx_ig
                d2udx2_val = theta*d2udx2_ig + (1.0-theta)*d2updx2_ig
                
                tau_ig = fe_stab.tau_stabilization(t, x_ig, u_ig, h_el, dt)

                r_strong = fe_stab.f_strong(t, x_ig, u_val, dudt_ig, dudx_val, d2udx2_val)

                # res_ig = np.zeros(nbf*neq)
                for inod in range(nbf):

                    gnod = gnodes[inod]
                    row1 = neq*gnod
                    rowf = row1+neq

                    w = bfn_ig[inod]
                    dwdx = dfdx_ig[inod]
                    d2wdx2 = d2fdx2_ig[inod]
                    dwdt = 0.0

                    res_ig = f(t_val, x_ig, u_val, dudt_ig, dudx_val, w, dwdx)

                    res_ig += fe_stab.stabilization_terms(t_val, x_ig, u_val, r_strong,
                                                        w, dwdt, dwdx, d2wdx2, tau_ig)

                    res[row1:rowf] += res_ig*w_ig

        # res[:neq] = self.bc['left'](t, x[0], u[0:neq], np.zeros(neq))
        # res[-neq:] = self.bc['right'](t, x[-1], u[-neq:], np.zeros(neq))

        left_bc = self.bc['left'](t, x[0], u[0:neq], np.zeros(neq))
        right_bc = self.bc['right'](t, x[-1], u[-neq:], np.zeros(neq))

        for ieq in range(neq):

            if self.is_Dirichlet_left[ieq]:
                res[ieq] = left_bc[ieq]
            else:
                res[ieq] -= left_bc[ieq]
                
            r_idx = -neq + ieq
            if self.is_Dirichlet_right[ieq]:
                res[r_idx] = right_bc[ieq]
            else:
                res[r_idx] += right_bc[ieq]

        return res
    
    def elemental_arrays(self, u: np.ndarray[tuple[int]],
        gnodes: np.ndarray[tuple[int]]):
        """Returns the arrays with the values of the mesh and the unknowns at the
        local element.

        Args:
            self: the FiniteElementsBVP1D class object
            u: the global array of unknowns
            gnodes: the global enumeration of the nodes of the element
        Returns:
            The arrays of values (x, u, up, u2p) at the element nodes."""

        neq = self.neq

        up = self.time_int.up
        u2p = self.time_int.up2

        nbf = self.gauss_int.nbf
        
        x_el = self.x[gnodes]

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

    def det_jac(self, x_el: np.ndarray[tuple[int]],
        dfdc: np.ndarray[tuple[int]]) -> float:
        """Returns the transformation jacobian from the parent to the
        physical element.

        Args:
            self: the FiniteElementsBVP1D class object
            x_el: the array of mesh values at the element nodes
            dfdc: the value of the derivative of the basis function with
            respect to the parent coordinate, c, at each node of the element.
        Returns:
            The arrays of values (x, u, up, u2p) at the element nodes."""
        
        dxdc_ig = np.dot(x_el, dfdc)
        
        jac_ig = dxdc_ig
        return jac_ig

    def jacobian(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-8) -> np.ndarray[tuple[int, int]]:

        return super().jacobian(u, res, h)

class FiniteDifferencesBVP(NonlinearProblem):

    def __init__(self, *args, neq: int, mesh: MeshDiscretization,
        bc: list[list[Callable]], time_int: TimeIntegration,
        theta: float = 1.0, **kwargs):

        super().__init__(*args, **kwargs)
        self.neq = neq
        self.mesh = mesh
        self.bc = bc
        self.time_int = time_int
        self.theta = theta

    def f_res(self, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

        f = self.f

        neq = self.neq
        theta = self.theta

        x = self.mesh.x_mesh
        connectivity = self.mesh.connectivity
        ndim = self.mesh.ndim
        nodes_total = self.mesh.nodes_total

        dt = self.time_int.dt
        t = self.time_int.t
        up = self.time_int.up
            
        ti = theta*t + (1.0-theta)*(t-dt)

        nunknowns = u.shape[0]

        res = np.zeros(nunknowns)

        for i in range(nodes_total):

            row1 = i*neq
            rowf = row1+neq

            # boundary_bcs = []
            # for d in range(ndim):
            #     for j in range(connectivity.shape[2]):
            #         if connectivity[i, d, j] == -1:
            #             boundary_bcs.append((j, d, self.bc[j][d]))
            #             is_boundary = True

            is_boundary = False
            # the highest dimension index dominates corners
            for d in range(ndim):
                for j in range(connectivity.shape[2]):
                    if connectivity[i, d, j] == -1:
                        bc = self.bc[j][d]
                        res[row1:rowf] = self.bc_residual(u, i, bc)
                        is_boundary = True
                        break
            
            if is_boundary:
                continue

            ui = u[row1:rowf]
            ui_p = up[row1:rowf]
            u_val = theta*ui + (1.0-theta)*ui_p
            
            dudt, grad_u, hess_u = self.calculate_derivatives(u, i)

            res[row1:rowf] = f(ti, x[i], u_val, dudt, grad_u, hess_u) # implicit form

        return res
    
    def calculate_derivatives(self, u: np.ndarray[tuple[int]], i: int):

        neq = self.neq
        theta = self.theta

        up = self.time_int.up
        dt = self.time_int.dt
        dtp = self.time_int.dtp

        x = self.mesh.x_mesh
        connectivity = self.mesh.connectivity
        ndim = self.mesh.ndim
        
        row1 = i*neq
        rowf = row1+neq
        
        ui = u[row1:rowf]
        ui_p = up[row1:rowf]
        ui_2p = self.time_int.up2[row1:rowf]
        if self.time_int.inc == 1:
            ui_2p = 2*up - u
        
        o1_h2_p_i = (2*dt+dtp)/(dt*(dt+dtp))  # 1st order (h2) backward, ui
        o1_h2_p_p = -(dt+dtp)/(dt*dtp)        # 1st order (h2) backward, up
        o1_h2_p_2p = dt/(dtp*(dt+dtp))        # 1st order (h2) backward, u2p

        # dudt = (ui - ui_p)/dt # O(h) backward
        # dudt = (3.0*ui - 4.0*ui_p + ui_2p)/(2.0*dt) # O(h2) backward, dt: constant
        dudt = (o1_h2_p_i*ui + o1_h2_p_p*ui_p + o1_h2_p_2p*ui_2p) # O(h2) backward

        grad_u = np.zeros((ndim, neq))
        hess_u = np.zeros((ndim, ndim, neq))
        
        for d in range(ndim):
        
            idx_b, idx_f = connectivity[i,d,0], connectivity[i,d,1]

            hf = x[idx_f,d] - x[i,d]
            hb = x[i,d] - x[idx_b,d]

            ub = u[idx_b*neq : idx_b*neq + neq]
            uf = u[idx_f*neq : idx_f*neq + neq]
            ub_p = up[idx_b*neq : idx_b*neq + neq]
            uf_p = up[idx_f*neq : idx_f*neq + neq]

            # idx_2b = connectivity[idx_b,d,0]
            # h2b = x[idx_b,d] - x[idx_2b,d]
            # o1_h2_b_i = (2*hb+h2b)/(hb*(hb+h2b)) # 1st order (h2) backward, ui
            # o1_h2_b_b = -(hb+h2b)/(hb*h2b)       # 1st order (h2) backward, ub
            # o1_h2_b_2b = hb/(h2b*(hb+h2b))       # 1st order (h2) backward, u2b
            # dudx = (ui - ub)/hb # O(h) backward
            # u2b = u[idx_2b*neq : idx_2b*neq + neq]
            # dudx = o1_h2_b_i*ui + o1_h2_b_b*ub + o1_h2_b_2b*u2b # O(h2) backward
            # dudx_p = (ui_p - ub_p)/hb # O(h) backward
            # u2b_p = up[idx_2b*neq : idx_2b*neq + neq]
            # dudx_p = o1_h2_b_i*ui_p + o1_h2_b_b*ub_p + o1_h2_b_2b*u2b_p # O(h2) backward

            # idx_2f = connectivity[idx_f,d,1]
            # h2f = x[idx_2f,d] - x[idx_f,d]
            # o1_h2_f_i = -(2*hf+h2f)/(hf*(hf+h2f)) # 1st order (h2) forward, ui
            # o1_h2_f_f = (hf+h2f)/(hf*h2f)         # 1st order (h2) forward, uf
            # o1_h2_f_2f = -hf/(h2f*(hf+h2f))       # 1st order (h2) forward, u2f
            dudx = (uf - ui)/hf # O(h) forward
            # u2f = u[idx_2f*neq : idx_2f*neq+neq]
            # dudx = o1_h2_f_i*ui + o1_h2_f_f*uf + o1_h2_f_2f*u2f # O(h2) forward
            dudx_p = (uf_p - ui_p)/hf # O(h) forward
            # u2f_p = up[idx_2f*neq : idx_2f*neq+neq]
            # dudx_p = o1_h2_f_i*ui_p + o1_h2_f_f*uf_p + o1_h2_f_2f*u2f_p # O(h2) forward

            # o1_h2_c = 1.0/(hf*hb*(hf+hb)) # 1st order (h2) central
            # dudx = (uf - ub)/(hf+hb) # O(h) central
            # dudx = ((hb**2)*uf + (hf**2-hb**2)*ui + (hf**2)*ub)*o1_h2_c # O(h2) central
            # dudx_p = (uf_p - ub_p)/(hf+hb) # O(h) central
            # dudx_p = ((hb**2)*uf_p + (hf**2-hb**2)*ui_p + (hf**2)*ub_p)*o1_h2_c # O(h2) central

            o2_h2_c = 2.0/(hf+hb) # 2nd order (h2) central
            d2udx2 = ((uf-ui)/hf - (ui-ub)/hb)*o2_h2_c # O(h2) central
            d2udx2_p = ((uf_p-ui_p)/hf - (ui_p-ub_p)/hb)*o2_h2_c # O(h2) central

            grad_u[d,:] = theta*dudx + (1.0-theta)*dudx_p
            hess_u[d,d,:] = theta*d2udx2 + (1.0-theta)*d2udx2_p
            
            for dj in range(d+1,ndim):

                idx_b_j, idx_f_j = connectivity[i,dj,0], connectivity[i,dj,1]

                hf_j = x[idx_f_j,dj] - x[i,dj]
                hb_j = x[i,dj] - x[idx_b_j,dj]

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
    
    def bc_residual(self, u: np.ndarray[tuple[int]],
        i: int, bc: Callable) -> np.ndarray[tuple[int]]:
        """Returns the residual of a boundary condition along a boundary.

        Args:
            u: the vector of unknowns
            i: the node of the boundary
            bc: the function of the boundary condition residual
        Returns:
            The residual of a boundary condition along a boundary."""
        
        neq = self.neq
        ndim = self.mesh.ndim
        x = self.mesh.x_mesh
        t = self.time_int.t

        connectivity = self.mesh.connectivity
        
        row1, rowf = i*neq, i*neq+neq

        ui = u[row1:rowf]

        grad_u = np.zeros((ndim, neq))

        for d in range(ndim):
            
            if connectivity[i,d,1] != -1:

                idx_f = connectivity[i,d,1]
                hf = x[idx_f,d] - x[i,d]
                
                uf = u[idx_f*neq : idx_f*neq+neq]

                idx_2f = connectivity[idx_f,d,1]
                if idx_2f == -1:
                    grad_u[d,:] = (uf - ui)/hf
                else:
                    h2f = x[idx_2f,d] - x[idx_f,d]

                    o1_h2_f_i = -(2*hf+h2f)/(hf*(hf+h2f)) # 1st order (h2) forward, ui
                    o1_h2_f_f = (hf+h2f)/(hf*h2f)         # 1st order (h2) forward, uf
                    o1_h2_f_2f = -hf/(h2f*(hf+h2f))       # 1st order (h2) forward, u2f

                    u2f = u[idx_2f*neq : idx_2f*neq+neq]

                    grad_u[d,:] = o1_h2_f_i*ui + o1_h2_f_f*uf + o1_h2_f_2f*u2f

            elif connectivity[i,d,0] != -1:

                idx_b = connectivity[i,d,0]
                hb = x[i,d] - x[idx_b,d]
                
                ub = u[idx_b*neq : idx_b*neq + neq]
                
                idx_2b = connectivity[idx_b,d,0]
                if idx_2b == -1:
                    grad_u[d, :] = (ui - ub)/hb
                else:
                    h2b = x[idx_b,d] - x[idx_2b,d]

                    o1_h2_b_i = (2*hb+h2b)/(hb*(hb+h2b)) # 1st order (h2) backward, ui
                    o1_h2_b_b = -(hb+h2b)/(hb*h2b)       # 1st order (h2) backward, ub
                    o1_h2_b_2b = hb/(h2b*(hb+h2b))       # 1st order (h2) backward, u2b

                    u2b = u[idx_2b*neq : idx_2b*neq + neq]

                    grad_u[d,:] = o1_h2_b_i*ui + o1_h2_b_b*ub + o1_h2_b_2b*u2b

        return bc(t, x[i], ui, grad_u)

    def jacobian(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-8) -> np.ndarray[tuple[int, int]]:

        return super().jacobian(u, res, h)
