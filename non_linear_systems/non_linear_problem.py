
from typing import Callable
from abc import ABC, abstractmethod
import numpy as np
import numpy.typing as npt

class NonlinearProblem(ABC):

    def __init__(self,
        f: Callable[[np.ndarray[tuple[int], np.dtype[np.float64]]], npt.NDArray[np.float64]]):

        self.f = f

    @abstractmethod
    def f_res(self, u: np.ndarray[tuple[int], np.dtype[np.float64]]
        ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:
        """Returns the residual vector of the non-linear system.

        Args:
            u: the vector of unknowns
        Returns:
            The residual vector at the values of the unknowns, u."""

        pass

    @abstractmethod
    def jacobian(self, u: np.ndarray[tuple[int], np.dtype[np.float64]],
        res: np.ndarray[tuple[int], np.dtype[np.float64]], h: float = 1e-8
        ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:
        """Returns the Jacobian of a system of non-linear algebraic equations
        around the values, u, of the unknowns.

        Args:
            res: the residual vector at the values of the unknowns, u
            f_res: the function that calculates the residual
            u: the unknown vector u
        Returns:
            The Jacobian of the non-linear system around the unknown vector, u."""

        n, m = len(u), len(res)

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

    def f_res(self, u: np.ndarray[tuple[int], np.dtype[np.float64]]
        ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:

        return self.f(u)

    def jacobian(self, u: np.ndarray[tuple[int], np.dtype[np.float64]],
        res: np.ndarray[tuple[int], np.dtype[np.float64]], h: float = 1e-8
        ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:

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

    def f_res(self, x: np.ndarray[tuple[int], np.dtype[np.float64]]
        ) -> np.ndarray[tuple[int], np.dtype[np.float64]]:

        return self.grad_f(self.df, self.f, x, h=1e-6)

    def jacobian(self, x: np.ndarray[tuple[int], np.dtype[np.float64]],
        res: np.ndarray[tuple[int], np.dtype[np.float64]] = None, h: float = 1e-6
        ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:

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

    def f_res(self, u: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:

        nx = self.nx

        x = u[:nx]
        l_m = u[nx:]

        r2 = self.g(x)

        grad_f = self.grad_f(self.df, self.f, x)

        jac_g = self.jacobian_g(x, r2)

        r1 =  grad_f + np.dot(l_m, jac_g)

        # return np.concatenate([r1, r2])

        r_tot = np.zeros(nx+len(l_m))
        r_tot[:nx] = r1.copy()
        r_tot[nx:] = r2.copy()
        return r_tot


    def jacobian(self, u: np.ndarray[tuple[int], np.dtype[np.float64]],
        res: np.ndarray[tuple[int], np.dtype[np.float64]], h: float = 1e-8
        ) -> np.ndarray[tuple[int, int], np.dtype[np.float64]]:

        return super().jacobian(u, res, h)

    def jacobian_g(self, u: npt.NDArray[np.float64], res: npt.NDArray[np.float64],
        h: float = 1e-4) -> npt.NDArray[np.float64]:
        """Returns the Jacobian of the constraints, g, around the unknowns values, u."""

        n, m = len(u), len(res)

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

    def __init__(self, *args, xi: npt.NDArray[np.float64],
        yi: npt.NDArray[np.float64], **kwargs):

        super().__init__(*args, **kwargs)
        self.xi = xi
        self.yi = yi

    def f_res(self, u: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return self.yi - self.f(u, self.xi)


class PDESystem(NonlinearProblem):

    def __init__(self, h, conductivity):
        self.h = h
        self.k = conductivity

    def f_res(self, u):
        return (u**2)/self.h + self.k
