
from typing import Callable
from abc import ABC, abstractmethod
import numpy as np
import numpy.typing as npt

class NonlinearProblem(ABC):

    def __init__(self, f: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]]):
        self.f = f

    @abstractmethod
    def f_res(self, u: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        pass


class Regular(NonlinearProblem):

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def f_res(self, u: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return self.f(u)


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

    def f_res(self, x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return self.grad_f(self.df, self.f, x)

    def jac(self, x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:
        return self.hessian_f(self.df, self.f, x)


class LagrangeMultiplier(NonlinearProblem):

    def __init__(self, *args, nu: int,
        g: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
        df: Callable[[Callable[[float], float], float, float, float], float],
        grad_f: Callable[[Callable[[Callable[[float], float], float, float, float], float],
        Callable[[npt.NDArray[np.float64]], float], npt.NDArray[np.float64],
        float], npt.NDArray[np.float64]], **kwargs):

        super().__init__(*args, **kwargs)
        self.nu = nu
        self.g = g
        self.df = df
        self.grad_f = grad_f

    def f_res(self, x: npt.NDArray[np.float64]) -> npt.NDArray[np.float64]:

        n = self.nu

        u = x[:n]
        l = x[n:]

        r2 = self.g(u)

        grad_f = self.grad_f(self.df, self.f, u)
        jac_g = self.jac_g(r2, self.g, u)
        r1 =  grad_f + np.dot(l, jac_g)

        return np.concatenate([r1, r2])

    @staticmethod
    def jac_g(res: npt.NDArray[np.float64],
        f_res: Callable[[npt.NDArray[np.float64]], npt.NDArray[np.float64]],
        u: npt.NDArray[np.float64], h: float = 1e-4) -> npt.NDArray[np.float64]:
        """Returns the Jacobian of the constraints, g, around the unknowns values, u."""

        n, m = len(u), len(res)

        jac = np.zeros((m,n))
        for j in range(n):

            u_val = u[j]
            u[j] += h
            res_f = f_res(u)
            u[j] = u_val

            # u_val = u[j]
            # u[j] -= h
            # res_b = f_res(u)
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
