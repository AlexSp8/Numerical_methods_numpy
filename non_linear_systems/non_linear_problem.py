
from typing import Callable
from abc import ABC, abstractmethod
import numpy as np
import numpy.typing as npt

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

    def f_res(self, yn: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

        neq = (self.yi).shape[0]
        ns = self.method_dict['ns']
        q = self.method_dict['q']
        p = self.method_dict['p']

        # Ys = yn.reshape((ns, neq))
        # k_j = np.zeros((ns,neq))
        # for j in range(ns):
        #     k_j[j] = self.f(self.ti + p[j]*self.h, Ys[j])
        # res_matrix = Ys - self.yi - self.h*(np.matmul(q,k_j))
        # return res_matrix.flatten()

        Ys = np.zeros((ns,neq))
        k_j = np.zeros((ns,neq))
        for j in range(ns):
            Ys[j] = yn[j*neq:(j+1)*neq]
            k_j[j] = self.f(self.ti + p[j]*self.h, Ys[j])

        res = np.zeros(ns*neq)
        for j in range(ns):
            for ieq in range(neq):
                row = j*neq + ieq
                s_sum = np.dot(q[j,:], k_j[:,ieq])
                res[row] = Ys[j,ieq] - self.yi[ieq] - self.h*s_sum

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

        f_val = self.f(self.ti+self.h, yn)
        return yn - self.yi - self.h*self.b0*f_val

    def jacobian(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-8) -> np.ndarray[tuple[int, int]]:

        return super().jacobian(u, res, h)


class BVP1D(NonlinearProblem):

    def __init__(self, *args, x: np.ndarray[tuple[int]],
        bc: dict, neq: int, **kwargs):

        super().__init__(*args, **kwargs)
        self.x = x
        self.bc = bc
        self.neq = neq

    def f_res(self, u: np.ndarray[tuple[int]]) -> np.ndarray[tuple[int]]:

        nnodes, nunknowns = (self.x).shape[0], u.shape[0]

        neq = self.neq

        res = np.zeros(nunknowns)

        for i in range(1, nnodes-1):

            hf = self.x[i+1] - self.x[i]
            hb = self.x[i] - self.x[i-1]

            row1 = i*neq

            ui = u[row1:row1+neq]

            t2_hf_hb = 2.0/(hf+hb)

            # 1st, 2nd derivatives
            dui = np.zeros(neq)
            d2ui = np.zeros(neq)
            for jeq in range(neq):

                row = row1 + jeq

                dui[jeq] = (u[row+neq] - u[row])/hf
                # dui[jeq] = (u[row+neq] - u[row-neq])/(hb+hf)

                d2ui[jeq] = (u[row+neq] - u[row])/hf
                d2ui[jeq] -= (u[row] - u[row-neq])/hb
                d2ui[jeq] *= t2_hf_hb

            # source term
            # fi = self.f(self.x[i], ui, dui) # explicit form
            fi = self.f(self.x[i], ui, dui, d2ui) # implicit form

            for jeq in range(neq):
                row = row1 + jeq
                # res[row] = d2ui[jeq] + fi[jeq] # explicit form
                res[row] = fi[jeq] # implicit form

        hf = self.x[1]-self.x[0]
        hb = self.x[-1]-self.x[-2]
        for j in range(neq):

            res[j] = self.left_boundary_condition(u, j, hf)

            row = (nnodes-1)*neq+j
            res[row] = self.right_boundary_condition(u, j, hb)

        return res

    def left_boundary_condition(self, u: np.ndarray[tuple[int]],
        ieq: int, hf: float) -> float:
        """Returns the residual on the left boundary of a 1D domain.

        Args:
            u: the vector of unknowns
            ieq: the equation on which the boundary condition is applied
            hf: the distance between the first two points of the domain (for derivative BCs)
        Returns:
            The residual on the left boundary after applying the boundary condition."""

        bc = self.bc['left']

        neq = self.neq

        res_bc = 0.0

        node1, node2 = ieq, ieq+neq

        bc_type = bc['type']
        val = bc['value'][ieq]

        if bc_type == 'dirichlet':
            res_bc = u[node1] - val
        elif bc_type == 'neumann':
            res_bc = (u[node2] - u[node1])/hf - val
        elif bc_type == 'robin':
            a0 = bc.get('a_robin', [0.0]*neq)[ieq]
            res_bc = (u[node2] - u[node1])/hf +a0*u[node1] - val

        return res_bc

    def right_boundary_condition(self, u: np.ndarray[tuple[int]],
        ieq: int, hb: float) -> float:
        """Returns the residual on the right boundary of a 1D domain.

        Args:
            u: the vector of unknowns
            ieq: the equation on which the boundary condition is applied
            hf: the distance between the first two points of the domain (for derivative BCs)
        Returns:
            The residual on the right boundary after applying the boundary condition."""

        bc = self.bc['right']

        neq = self.neq

        nunknowns = u.shape[0]

        node1 = nunknowns-neq+ieq
        node2 = node1-neq

        res_bc = 0.0

        bc_type = bc['type']
        val = bc['value'][ieq]

        if bc_type == 'dirichlet':
            res_bc = u[node1] - val
        elif bc_type == 'neumann':
            res_bc = (u[node1] - u[node2])/hb - val
        elif bc_type == 'robin':
            an = bc.get('a_robin', [0.0]*neq)[ieq]
            res_bc = (u[node1] - u[node2])/hb +an*u[node1] - val

        return res_bc

    def jacobian(self, u: np.ndarray[tuple[int]], res: np.ndarray[tuple[int]],
        h: float = 1e-8) -> np.ndarray[tuple[int, int]]:

        return super().jacobian(u, res, h)


class PDESystem(NonlinearProblem):

    def __init__(self, h, conductivity):
        self.h = h
        self.k = conductivity

    def f_res(self, u):
        return (u**2)/self.h + self.k
