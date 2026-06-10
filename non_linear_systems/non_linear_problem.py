
from typing import Callable
from abc import ABC, abstractmethod
import numpy as np
import numpy.typing as npt

from ode.time_integration import TimeIntegration

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


class TransientBVP1DFD(NonlinearProblem):

    def __init__(self, *args, x: np.ndarray[tuple[int]], bc: dict[str, Callable],
        time_int: TimeIntegration, theta: float = 1.0, **kwargs):

        super().__init__(*args, **kwargs)
        self.x = x
        self.bc = bc
        self.time_int = time_int
        self.theta = theta

    def f_res(self, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

        f = self.f
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

        neq = round(nunknowns/nnodes)

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
            row_f = row1+neq

            ub_i = u[row1-neq:row1]
            u_i = u[row1:row_f]
            uf_i = u[row_f:row_f+neq]

            up_i = up[row1:row_f]
            u2p_i = u2p[row1:row_f]
            
            upb_i = up[row1-neq:row1]
            upf_i = up[row_f:row_f+neq]

            # dudt_i = (u_i - up_i)/dt # O(h) backward
            # dudt_i = (3.0*u_i - 4.0*up_i + up2_i)/(2.0*dt) # O(h2) backward, dt: constant
            dudt_i = (o1_h2_p_i*u_i + o1_h2_p_p*up_i + o1_h2_p_2p*u2p_i) # O(h2) backward

            dudx_i = (uf_i - u_i)/hf # O(h) forward
            # dudx_i = (uf_i - ub_i)/(hf+hb) # O(h) central
            # dudx_i = ((hb**2)*uf_i + (hf**2-hb**2)*u_i + (hf**2)*ub_i)*o1_h2_c # O(h2) central
            dudx_i *= theta
            
            dupdx_i = (upf_i - up_i)/hf # O(h) forward
            # dupdx_i = (upf_i - upb_i)/(hf+hb) # O(h) central
            # dupdx_i = ((hb**2)*upf_i + (hf**2-hb**2)*up_i + (hf**2)*upb_i)*o1_h2_c # O(h2) central
            dudx_i += (1.0-theta)*dupdx_i

            d2udx2_i = ((uf_i-u_i)/hf - (u_i-ub_i)/hb)*o2_h2_c # O(h2) central
            d2udx2_i *= theta
            
            d2updx2_i = ((upf_i-up_i)/hf - (up_i-upb_i)/hb)*o2_h2_c # O(h2) central
            d2udx2_i += (1.0-theta)*d2updx2_i

            u_theta_i = theta*u_i + (1.0-theta)*up_i

            t_theta = theta*t + (1.0-theta)*(t-dt)

            res[row1:row_f] = f(t_theta, x[i], u_theta_i, dudt_i, dudx_i, d2udx2_i) # implicit form

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

        x = self.x
        bc = self.bc['left']

        nnodes, nunknowns = x.shape[0], u.shape[0]
        neq = round(nunknowns/nnodes)

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

        x = self.x
        bc = self.bc['right']

        nnodes, nunknowns = x.shape[0], u.shape[0]
        neq = round(nunknowns/nnodes)

        hb = x[-1]-x[-2]
        h2b = x[-2] - x[-3]
        o1_h2_b_i = (2*hb+h2b)/(hb*(hb+h2b)) # 1st order (h2) backward, ui
        o1_h2_b_b = -(hb+h2b)/(hb*h2b)       # 1st order (h2) backward, ub
        o1_h2_b_2b = hb/(h2b*(hb+h2b))       # 1st order (h2) backward, u2b

        ub = u[nunknowns-neq:]
        ub_b = u[nunknowns-2*neq:nunknowns-neq]
        ub_2b = u[nunknowns-3*neq:nunknowns-2*neq]

        # dudx_b = (ub - ub_b)/hb
        # dudx_b = (3.0*ub - 4.0*ub_b + ub_2b)/(2.0*hb)
        dudx_b = o1_h2_b_i*ub + o1_h2_b_b*ub_b + o1_h2_b_2b*ub_2b

        return bc(t, x[-1], ub, dudx_b)

    def jacobian(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-8) -> np.ndarray[tuple[int, int]]:

        return super().jacobian(u, res, h)
