
import numpy as np
from numpy.typing import NDArray
from abc import ABC, abstractmethod

def gauss_coefficients_1d(ng: int = 1, k_max = 100,
    ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Solves for Gauss-Legendre points and weights of degree n."""

    points = np.zeros(ng)
    weights = np.zeros(ng)

    # We only need to find half the roots due to symmetry
    m = (ng + 1)//2

    for i in range(m):
        # Initial guess for the i-th root (Chebyshev-like spacing)
        x = np.cos(np.pi*(i + 0.75)/(ng + 0.5))

        for k in range(1, k_max+1):
            p, dp = legendre_polynomial(x, ng)
            dx = p/dp
            x -= dx
            if abs(dx) < 1e-12:
                break

        if k >= k_max:
            print(f"Warning: Newton maximum iterations ({k}) reached without converging.")

        points[i] = -x
        points[ng-1-i] = x

        weight = 2.0/((1.0 - x**2)*(dp**2))
        weights[i] = weight
        weights[ng-1-i] = weight

    return points, weights

def legendre_polynomial(x: float, ng: int) -> tuple[float, float]:
    """Evaluates the n-th Legendre polynomial and its derivative at x
    using the recursive relation."""

    if ng == 0:
        return 1.0, 0.0

    p_prev, p = 0.0, 1.0

    for i in range(ng):
        p_next = ((2.0*i + 1.0)*x*p - i*p_prev)/(i + 1.0)
        p_prev = p
        p = p_next

    dp = ng*(x*p - p_prev)/(x**2 - 1)

    return p, dp


class FEGaussIntegration(ABC):

    def __init__(self, nbf: int):
        self.nbf = nbf

    @abstractmethod
    def set_linear_basis_functions(self):
        pass

    @abstractmethod
    def set_quadratic_basis_functions(self):
        pass


class FEGaussIntegrationLine(FEGaussIntegration):

    def __init__(self, ng: int, *args, **kwargs):

        super().__init__(*args, **kwargs)

        self.c, self.w = gauss_coefficients_1d(ng)

        self.ng = ng

        nbf = self.nbf
        ndim = 1

        self.bfn = np.zeros((nbf, ng))
        self.grad_fc = np.zeros((ndim, nbf, ng))
        self.hess_fc = np.zeros((ndim, ndim, nbf, ng))

        if nbf == 2:
            self.set_linear_basis_functions()
        elif nbf == 3:
            self.set_quadratic_basis_functions()
        else:
            raise ValueError(f"Unsupported nbf = {nbf} in FEGaussIntegrationLine!")

    def set_linear_basis_functions(self):

        c = self.c

        self.bfn[0,:] = (1.0-c)/2.0
        self.bfn[1,:] = (1.0+c)/2.0

        # dfdx
        self.grad_fc[0,0,:] = -1.0/2.0
        self.grad_fc[0,1,:] = +1.0/2.0

    def set_quadratic_basis_functions(self):

        c = self.c

        self.bfn[0,:] = -(1.0-c)*c/2.0
        self.bfn[1,:] = +(1.0-c)*(1.0+c) # mid node
        self.bfn[2,:] = +(1.0+c)*c/2.0

        # dfdx
        self.grad_fc[0,0,:] = -(1.0-2.0*c)/2.0
        self.grad_fc[0,1,:] = -2.0*c # mid node
        self.grad_fc[0,2,:] = +(1.0+2.0*c)/2.0

        # d2fdx2
        self.hess_fc[0,0,0,:] = +1.0
        self.hess_fc[0,0,1,:] = -2.0 # mid node
        self.hess_fc[0,0,2,:] = +1.0


class FEGaussIntegrationQuadrangle(FEGaussIntegration):

    def __init__(self, ng_1d: int, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.ng = ng_1d**2

        c_1d, w_1d = gauss_coefficients_1d(ng_1d)

        self.c, self.w = self.set_gauss_coefficients_2d(ng_1d, c_1d, w_1d)

        ng = self.ng
        nbf = self.nbf
        ndim = 2

        self.bfn = np.zeros((nbf, ng))
        self.grad_fc = np.zeros((ndim, nbf, ng))
        self.hess_fc = np.zeros((ndim, ndim, nbf, ng))

        if nbf == 4:
            self.set_linear_basis_functions()
        elif nbf == 9:
            self.set_quadratic_basis_functions()
        elif nbf == 8:
            self.set_serendipity_basis_functions()
        else:
            raise ValueError(f"Unsupported nbf = {nbf} in FEGaussIntegrationQuadrangle!")

    def set_gauss_coefficients_2d(self, ng_1d: int, c_1d: NDArray[np.float64],
        w_1d: NDArray[np.float64]) -> NDArray[np.float64]:

        c_2d = np.zeros((self.ng,2))
        w_2d = np.zeros(self.ng)

        ng = 0
        for i in range(ng_1d):
            for j in range(ng_1d):
                c_2d[ng,0] = c_1d[i]
                c_2d[ng,1] = c_1d[j]
                w_2d[ng] = w_1d[i]*w_1d[j]
                ng += 1

        # C, E = np.meshgrid(c_1d, c_1d, indexing='ij')
        # c_2d = np.stack((C.ravel(), E.ravel()), axis=-1)
        # w_2d np.outer(w_1d, w_1d).ravel()
        return c_2d, w_2d

    def set_linear_basis_functions(self):

        c = self.c[:,0]
        e = self.c[:,1]

        self.bfn[0,:] = (1.0-c)*(1.0-e)/4.0
        self.bfn[1,:] = (1.0+c)*(1.0-e)/4.0
        self.bfn[2,:] = (1.0+c)*(1.0+e)/4.0
        self.bfn[3,:] = (1.0-c)*(1.0+e)/4.0

        # dfdc
        self.grad_fc[0,0,:] = -(1.0-e)/4.0
        self.grad_fc[0,1,:] = +(1.0-e)/4.0
        self.grad_fc[0,2,:] = +(1.0+e)/4.0
        self.grad_fc[0,3,:] = -(1.0+e)/4.0

        # dfde
        self.grad_fc[1,0,:] = -(1.0-c)/4.0
        self.grad_fc[1,1,:] = -(1.0+c)/4.0
        self.grad_fc[1,2,:] = +(1.0+c)/4.0
        self.grad_fc[1,3,:] = +(1.0-c)/4.0

    def set_quadratic_basis_functions(self):

        c = self.c[:,0]
        e = self.c[:,1]

        self.bfn[0,:] = +(1.0-c)*c*(1.0-e)*e/4.0
        self.bfn[1,:] = -(1.0+c)*c*(1.0-e)*e/4.0
        self.bfn[2,:] = +(1.0+c)*c*(1.0+e)*e/4.0
        self.bfn[3,:] = -(1.0-c)*c*(1.0+e)*e/4.0
        self.bfn[4,:] = -(1.0+c)*(1.0-c)*(1.0-e)*e/2.0
        self.bfn[5,:] = +(1.0+c)*c*(1.0+e)*(1.0-e)/2.0
        self.bfn[6,:] = +(1.0+c)*(1.0-c)*(1.0+e)*e/2.0
        self.bfn[7,:] = -(1.0-c)*c*(1.0+e)*(1.0-e)/2.0
        self.bfn[8,:] = +(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)

        # dfdc
        self.grad_fc[0,0,:] = +(1.0-2.0*c)*(1.0-e)*e/4.0
        self.grad_fc[0,1,:] = -(1.0+2.0*c)*(1.0-e)*e/4.0
        self.grad_fc[0,2,:] = +(1.0+2.0*c)*(1.0+e)*e/4.0
        self.grad_fc[0,3,:] = -(1.0-2.0*c)*(1.0+e)*e/4.0
        self.grad_fc[0,4,:] = -(-2.0*c)*(1.0-e)*e/2.0
        self.grad_fc[0,5,:] = +(1.0+2.0*c)*(1.0+e)*(1.0-e)/2.0
        self.grad_fc[0,6,:] = +(-2.0*c)*(1.0+e)*e/2.0
        self.grad_fc[0,7,:] = -(1.0-2.0*c)*(1.0+e)*(1.0-e)/2.0
        self.grad_fc[0,8,:] = +(-2.0*c)*(1.0+e)*(1.0-e)

        # d2fdc2
        self.hess_fc[0,0,0,:] = +(-2.0)*(1.0-e)*e/4.0
        self.hess_fc[0,0,1,:] = -(+2.0)*(1.0-e)*e/4.0
        self.hess_fc[0,0,2,:] = +(+2.0)*(1.0+e)*e/4.0
        self.hess_fc[0,0,3,:] = -(-2.0)*(1.0+e)*e/4.0
        self.hess_fc[0,0,4,:] = -(-2.0)*(1.0-e)*e/2.0
        self.hess_fc[0,0,5,:] = +(+2.0)*(1.0+e)*(1.0-e)/2.0
        self.hess_fc[0,0,6,:] = +(-2.0)*(1.0+e)*e/2.0
        self.hess_fc[0,0,7,:] = -(-2.0)*(1.0+e)*(1.0-e)/2.0
        self.hess_fc[0,0,8,:] = +(-2.0)*(1.0+e)*(1.0-e)

        # d2fdce
        self.hess_fc[0,1,0,:] = +(1.0-2.0*c)*(1.0-2.0*e)/4.0
        self.hess_fc[0,1,1,:] = -(1.0+2.0*c)*(1.0-2.0*e)/4.0
        self.hess_fc[0,1,2,:] = +(1.0+2.0*c)*(1.0+2.0*e)/4.0
        self.hess_fc[0,1,3,:] = -(1.0-2.0*c)*(1.0+2.0*e)/4.0
        self.hess_fc[0,1,4,:] = -(-2.0*c)*(1.0-2.0*e)/2.0
        self.hess_fc[0,1,5,:] = +(1.0+2.0*c)*(-2.0*e)/2.0
        self.hess_fc[0,1,6,:] = +(-2.0*c)*(1.0+2.0*e)/2.0
        self.hess_fc[0,1,7,:] = -(1.0-2.0*c)*(-2.0*e)/2.0
        self.hess_fc[0,1,8,:] = +(-2.0*c)*(-2.0*e)
        self.hess_fc[1,0,:,:] = self.hess_fc[0,1,:,:]

        # dfde
        self.grad_fc[1,0,:] = +(1.0-c)*c*(1.0-2.0*e)/4.0
        self.grad_fc[1,1,:] = -(1.0+c)*c*(1.0-2.0*e)/4.0
        self.grad_fc[1,2,:] = +(1.0+c)*c*(1.0+2.0*e)/4.0
        self.grad_fc[1,3,:] = -(1.0-c)*c*(1.0+2.0*e)/4.0
        self.grad_fc[1,4,:] = -(1.0+c)*(1.0-c)*(1.0-2.0*e)/2.0
        self.grad_fc[1,5,:] = +(1.0+c)*c*(-2.0*e)/2.0
        self.grad_fc[1,6,:] = +(1.0+c)*(1.0-c)*(1.0+2.0*e)/2.0
        self.grad_fc[1,7,:] = -(1.0-c)*c*(-2.0*e)/2.0
        self.grad_fc[1,8,:] = +(1.0+c)*(1.0-c)*(-2.0*e)

        # d2fde2
        self.hess_fc[1,1,0,:] = +(1.0-c)*c*(-2.0)/4.0
        self.hess_fc[1,1,1,:] = -(1.0+c)*c*(-2.0)/4.0
        self.hess_fc[1,1,2,:] = +(1.0+c)*c*(+2.0)/4.0
        self.hess_fc[1,1,3,:] = -(1.0-c)*c*(+2.0)/4.0
        self.hess_fc[1,1,4,:] = -(1.0+c)*(1.0-c)*(-2.0)/2.0
        self.hess_fc[1,1,5,:] = +(1.0+c)*c*(-2.0)/2.0
        self.hess_fc[1,1,6,:] = +(1.0+c)*(1.0-c)*(+2.0)/2.0
        self.hess_fc[1,1,7,:] = -(1.0-c)*c*(-2.0)/2.0
        self.hess_fc[1,1,8,:] = +(1.0+c)*(1.0-c)*(-2.0)

    def set_serendipity_basis_functions(self):

        c = self.c[:,0]
        e = self.c[:,1]

        self.bfn[0,:] = -(1.0-c)*(1.0-e)*(1.0+c+e)/4.0
        self.bfn[1,:] = -(1.0+c)*(1.0-e)*(1.0-c+e)/4.0
        self.bfn[2,:] = -(1.0+c)*(1.0+e)*(1.0-c-e)/4.0
        self.bfn[3,:] = -(1.0-c)*(1.0+e)*(1.0+c-e)/4.0
        self.bfn[4,:] = +(1.0+c)*(1.0-c)*(1.0-e)/2.0
        self.bfn[5,:] = +(1.0+c)*(1.0+e)*(1.0-e)/2.0
        self.bfn[6,:] = +(1.0+c)*(1.0-c)*(1.0+e)/2.0
        self.bfn[7,:] = +(1.0-c)*(1.0+e)*(1.0-e)/2.0

        # dfdc
        self.grad_fc[0,0,:] = -((-1.0)*(1.0+c+e)+(1.0-c)*(+1.0))*(1.0-e)/4.0
        self.grad_fc[0,1,:] = -((+1.0)*(1.0-c+e)+(1.0+c)*(-1.0))*(1.0-e)/4.0
        self.grad_fc[0,2,:] = -((+1.0)*(1.0-c-e)+(1.0+c)*(-1.0))*(1.0+e)/4.0
        self.grad_fc[0,3,:] = -((-1.0)*(1.0+c-e)+(1.0-c)*(+1.0))*(1.0+e)/4.0
        self.grad_fc[0,4,:] = +(-2.0*c)*(1.0-e)/2.0
        self.grad_fc[0,5,:] = +(+1.0)*(1.0+e)*(1.0-e)/2.0
        self.grad_fc[0,6,:] = +(-2.0*c)*(1.0+e)/2.0
        self.grad_fc[0,7,:] = +(-1.0)*(1.0+e)*(1.0-e)/2.0

        # d2fdc2
        self.hess_fc[0,0,0,:] = -((-1.0)*(+1.0)+(-1.0)*(+1.0))*(1.0-e)/4.0
        self.hess_fc[0,0,1,:] = -((+1.0)*(-1.0)+(+1.0)*(-1.0))*(1.0-e)/4.0
        self.hess_fc[0,0,2,:] = -((+1.0)*(-1.0)+(+1.0)*(-1.0))*(1.0+e)/4.0
        self.hess_fc[0,0,3,:] = -((-1.0)*(+1.0)+(-1.0)*(+1.0))*(1.0+e)/4.0
        self.hess_fc[0,0,4,:] = +(-2.0)*(1.0-e)/2.0
        self.hess_fc[0,0,5,:] = +0.0
        self.hess_fc[0,0,6,:] = +(-2.0)*(1.0+e)/2.0
        self.hess_fc[0,0,7,:] = +0.0

        # d2fdce
        self.hess_fc[0,1,0,:] = -( -(+1.0)*(1.0-e)+(-(1.0+c+e)+(1.0-c))*(-1.0) )/4.0
        self.hess_fc[0,1,1,:] = -( +(+1.0)*(1.0-e)+(+(1.0-c+e)-(1.0+c))*(-1.0) )/4.0
        self.hess_fc[0,1,2,:] = -( +(-1.0)*(1.0+e)+(+(1.0-c-e)-(1.0+c))*(+1.0) )/4.0
        self.hess_fc[0,1,3,:] = -( -(-1.0)*(1.0+e)+(-(1.0+c-e)+(1.0-c))*(+1.0) )/4.0
        self.hess_fc[0,1,4,:] = +(-2.0*c)*(-1.0)/2.0
        self.hess_fc[0,1,5,:] = +(+1.0)*(-2.0*e)/2.0
        self.hess_fc[0,1,6,:] = +(-2.0*c)*(+1.0)/2.0
        self.hess_fc[0,1,7,:] = +(-1.0)*(-2.0*e)/2.0
        self.hess_fc[1,0,:,:] = self.hess_fc[0,1,:,:]

        # dfde
        self.grad_fc[1,0,:] = -(1.0-c)*((-1.0)*(1.0+c+e)+(1.0-e)*(+1.0))/4.0
        self.grad_fc[1,1,:] = -(1.0+c)*((-1.0)*(1.0-c+e)+(1.0-e)*(+1.0))/4.0
        self.grad_fc[1,2,:] = -(1.0+c)*((+1.0)*(1.0-c-e)+(1.0+e)*(-1.0))/4.0
        self.grad_fc[1,3,:] = -(1.0-c)*((+1.0)*(1.0+c-e)+(1.0+e)*(-1.0))/4.0
        self.grad_fc[1,4,:] = +(1.0+c)*(1.0-c)*(-1.0)/2.0
        self.grad_fc[1,5,:] = +(1.0+c)*(-2.0*e)/2.0
        self.grad_fc[1,6,:] = +(1.0+c)*(1.0-c)*(+1.0)/2.0
        self.grad_fc[1,7,:] = +(1.0-c)*(-2.0*e)/2.0

        # d2fde2
        self.hess_fc[1,1,0,:] = -(1.0-c)*(-(+1.0)+(-1.0))/4.0
        self.hess_fc[1,1,1,:] = -(1.0+c)*(-(+1.0)+(-1.0))/4.0
        self.hess_fc[1,1,2,:] = -(1.0+c)*(+(-1.0)-(+1.0))/4.0
        self.hess_fc[1,1,3,:] = -(1.0-c)*(+(-1.0)-(+1.0))/4.0
        self.hess_fc[1,1,4,:] = +0.0
        self.hess_fc[1,1,5,:] = +(1.0+c)*(-2.0)/2.0
        self.hess_fc[1,1,6,:] = +0.0
        self.hess_fc[1,1,7,:] = +(1.0-c)*(-2.0)/2.0


class FEGaussIntegrationHexahedron(FEGaussIntegration):

    def __init__(self, ng_1d: int, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.ng = ng_1d**3

        c_1d, w_1d = gauss_coefficients_1d(ng_1d)

        self.c, self.w = self.set_gauss_coefficients_3d(ng_1d, c_1d, w_1d)

        nbf = self.nbf
        ng = self.ng
        ndim = 3

        self.bfn = np.zeros((nbf, ng))
        self.grad_fc = np.zeros((ndim, nbf, ng))
        self.hess_fc = np.zeros((ndim, ndim, nbf, ng))

        if nbf == 8:
            self.set_linear_basis_functions()
        elif nbf == 27:
            self.set_quadratic_basis_functions()
        elif nbf == 20:
            self.set_serendipity_basis_functions()
        else:
            raise ValueError(f"Unsupported nbf = {nbf} in FEGaussIntegrationHexahedron!")

    def set_gauss_coefficients_3d(self, ng_1d: int, c_1d: NDArray[np.float64],
        w_1d: NDArray[np.float64]) -> NDArray[np.float64]:

        c_3d = np.zeros((self.ng,3))
        w_3d = np.zeros(self.ng)

        ng = 0
        for i in range(ng_1d):
            for j in range(ng_1d):
                for k in range(ng_1d):
                    c_3d[ng,0] = c_1d[i]
                    c_3d[ng,1] = c_1d[j]
                    c_3d[ng,2] = c_1d[k]
                    w_3d[ng] = w_1d[i]*w_1d[j]*w_1d[k]
                    ng += 1

        # C, E, S = np.meshgrid(c_1d, c_1d, c_1d, indexing='ij')
        # c_3d = np.stack((C.ravel(), E.ravel(), S.ravel()), axis=-1)
        # W1, W2, W3 = np.meshgrid(w_1d, w_1d, w_1d, indexing='ij')
        # w_3d = (W1*W2*W3).flatten()
        return c_3d, w_3d

    def set_linear_basis_functions(self):

        c = self.c[:,0]
        e = self.c[:,1]
        s = self.c[:,2]

        self.bfn[0,:] = (1.0-c)*(1.0-e)*(1.0-s)/8.0
        self.bfn[1,:] = (1.0+c)*(1.0-e)*(1.0-s)/8.0
        self.bfn[2,:] = (1.0+c)*(1.0+e)*(1.0-s)/8.0
        self.bfn[3,:] = (1.0-c)*(1.0+e)*(1.0-s)/8.0
        self.bfn[4,:] = (1.0-c)*(1.0-e)*(1.0+s)/8.0
        self.bfn[5,:] = (1.0+c)*(1.0-e)*(1.0+s)/8.0
        self.bfn[6,:] = (1.0+c)*(1.0+e)*(1.0+s)/8.0
        self.bfn[7,:] = (1.0-c)*(1.0+e)*(1.0+s)/8.0

        # dfdc
        self.grad_fc[0,0,:] = -(1.0-e)*(1.0-s)/8.0
        self.grad_fc[0,1,:] = +(1.0-e)*(1.0-s)/8.0
        self.grad_fc[0,2,:] = +(1.0+e)*(1.0-s)/8.0
        self.grad_fc[0,3,:] = -(1.0+e)*(1.0-s)/8.0
        self.grad_fc[0,4,:] = -(1.0-e)*(1.0+s)/8.0
        self.grad_fc[0,5,:] = +(1.0-e)*(1.0+s)/8.0
        self.grad_fc[0,6,:] = +(1.0+e)*(1.0+s)/8.0
        self.grad_fc[0,7,:] = -(1.0+e)*(1.0+s)/8.0

        # dfde
        self.grad_fc[1,0,:] = -(1.0-c)*(1.0-s)/8.0
        self.grad_fc[1,1,:] = -(1.0+c)*(1.0-s)/8.0
        self.grad_fc[1,2,:] = +(1.0+c)*(1.0-s)/8.0
        self.grad_fc[1,3,:] = +(1.0-c)*(1.0-s)/8.0
        self.grad_fc[1,4,:] = -(1.0-c)*(1.0+s)/8.0
        self.grad_fc[1,5,:] = -(1.0+c)*(1.0+s)/8.0
        self.grad_fc[1,6,:] = +(1.0+c)*(1.0+s)/8.0
        self.grad_fc[1,7,:] = +(1.0-c)*(1.0+s)/8.0

        # dfds
        self.grad_fc[2,0,:] = -(1.0-c)*(1.0-e)/8.0
        self.grad_fc[2,1,:] = -(1.0+c)*(1.0-e)/8.0
        self.grad_fc[2,2,:] = -(1.0+c)*(1.0+e)/8.0
        self.grad_fc[2,3,:] = -(1.0-c)*(1.0+e)/8.0
        self.grad_fc[2,4,:] = +(1.0-c)*(1.0-e)/8.0
        self.grad_fc[2,5,:] = +(1.0+c)*(1.0-e)/8.0
        self.grad_fc[2,6,:] = +(1.0+c)*(1.0+e)/8.0
        self.grad_fc[2,7,:] = +(1.0-c)*(1.0+e)/8.0

    def set_quadratic_basis_functions(self):

        c = self.c[:,0]
        e = self.c[:,1]
        s = self.c[:,2]

        # Corners
        self.bfn[0,:] = -(1.0-c)*c*(1.0-e)*e*(1.0-s)*s/8.0
        self.bfn[1,:] = +(1.0+c)*c*(1.0-e)*e*(1.0-s)*s/8.0
        self.bfn[2,:] = -(1.0+c)*c*(1.0+e)*e*(1.0-s)*s/8.0
        self.bfn[3,:] = +(1.0-c)*c*(1.0+e)*e*(1.0-s)*s/8.0
        self.bfn[4,:] = +(1.0-c)*c*(1.0-e)*e*(1.0+s)*s/8.0
        self.bfn[5,:] = -(1.0+c)*c*(1.0-e)*e*(1.0+s)*s/8.0
        self.bfn[6,:] = +(1.0+c)*c*(1.0+e)*e*(1.0+s)*s/8.0
        self.bfn[7,:] = -(1.0-c)*c*(1.0+e)*e*(1.0+s)*s/8.0
        # Mid-edge ξ = 0
        self.bfn[8,:]  = +(1.0+c)*(1.0-c)*(1.0-e)*e*(1.0-s)*s/4.0
        self.bfn[10,:] = -(1.0+c)*(1.0-c)*(1.0+e)*e*(1.0-s)*s/4.0
        self.bfn[12,:] = -(1.0+c)*(1.0-c)*(1.0-e)*e*(1.0+s)*s/4.0
        self.bfn[14,:] = +(1.0+c)*(1.0-c)*(1.0+e)*e*(1.0+s)*s/4.0
        # Mid-edge η = 0
        self.bfn[9,:] = -(1.0+c)*c*(1.0+e)*(1.0-e)*(1.0-s)*s/4.0
        self.bfn[11,:] = +(1.0-c)*c*(1.0+e)*(1.0-e)*(1.0-s)*s/4.0
        self.bfn[13,:] = +(1.0+c)*c*(1.0+e)*(1.0-e)*(1.0+s)*s/4.0
        self.bfn[15,:] = -(1.0-c)*c*(1.0+e)*(1.0-e)*(1.0+s)*s/4.0
        # Mid-edge ζ = 0
        self.bfn[16,:] = +(1.0-c)*c*(1.0-e)*e*(1.0+s)*(1.0-s)/4.0
        self.bfn[17,:] = -(1.0+c)*c*(1.0-e)*e*(1.0+s)*(1.0-s)/4.0
        self.bfn[18,:] = +(1.0+c)*c*(1.0+e)*e*(1.0+s)*(1.0-s)/4.0
        self.bfn[19,:] = -(1.0-c)*c*(1.0+e)*e*(1.0+s)*(1.0-s)/4.0
        # Mid-face ξ = +-1
        self.bfn[21,:] = +(1.0+c)*c*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)/2.0
        # self.bfn[22,:] = +(1.0+c)*c*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)/2.0  # Salome, ξ = +1
        self.bfn[23,:] = -(1.0-c)*c*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)/2.0
        # self.bfn[24,:] = -(1.0-c)*c*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)/2.0  # Salome, ξ = -1
        # Mid-face η = +-1
        self.bfn[20,:] = -(1.0+c)*(1.0-c)*(1.0-e)*e*(1.0+s)*(1.0-s)/2.0
        # self.bfn[21,:] = -(1.0+c)*(1.0-c)*(1.0-e)*e*(1.0+s)*(1.0-s)/2.0  # Salome, η = -1
        self.bfn[22,:] = +(1.0+c)*(1.0-c)*(1.0+e)*e*(1.0+s)*(1.0-s)/2.0
        # self.bfn[23,:] = +(1.0+c)*(1.0-c)*(1.0+e)*e*(1.0+s)*(1.0-s)/2.0  # Salome, η = +1
        # Mid-face ζ = +-1
        self.bfn[24,:] = -(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)*(1.0-s)*s/2.0
        # self.bfn[20,:] = -(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)*(1.0-s)*s/2.0  # Salome, ζ = -1
        self.bfn[25,:] = +(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)*(1.0+s)*s/2.0
        # Central
        self.bfn[26,:] = +(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)

        # dfdc
        # Corners
        self.grad_fc[0,0,:] = -(1.0-2.0*c)*(1.0-e)*e*(1.0-s)*s/8.0
        self.grad_fc[0,1,:] = +(1.0+2.0*c)*(1.0-e)*e*(1.0-s)*s/8.0
        self.grad_fc[0,2,:] = -(1.0+2.0*c)*(1.0+e)*e*(1.0-s)*s/8.0
        self.grad_fc[0,3,:] = +(1.0-2.0*c)*(1.0+e)*e*(1.0-s)*s/8.0
        self.grad_fc[0,4,:] = +(1.0-2.0*c)*(1.0-e)*e*(1.0+s)*s/8.0
        self.grad_fc[0,5,:] = -(1.0+2.0*c)*(1.0-e)*e*(1.0+s)*s/8.0
        self.grad_fc[0,6,:] = +(1.0+2.0*c)*(1.0+e)*e*(1.0+s)*s/8.0
        self.grad_fc[0,7,:] = -(1.0-2.0*c)*(1.0+e)*e*(1.0+s)*s/8.0
        # Mid-edge ξ = 0
        self.grad_fc[0,8,:]  = +(-2.0*c)*(1.0-e)*e*(1.0-s)*s/4.0
        self.grad_fc[0,10,:] = -(-2.0*c)*(1.0+e)*e*(1.0-s)*s/4.0
        self.grad_fc[0,12,:] = -(-2.0*c)*(1.0-e)*e*(1.0+s)*s/4.0
        self.grad_fc[0,14,:] = +(-2.0*c)*(1.0+e)*e*(1.0+s)*s/4.0
        # Mid-edge η = 0
        self.grad_fc[0,9,:] = -(1.0+2.0*c)*(1.0+e)*(1.0-e)*(1.0-s)*s/4.0
        self.grad_fc[0,11,:] = +(1.0-2.0*c)*(1.0+e)*(1.0-e)*(1.0-s)*s/4.0
        self.grad_fc[0,13,:] = +(1.0+2.0*c)*(1.0+e)*(1.0-e)*(1.0+s)*s/4.0
        self.grad_fc[0,15,:] = -(1.0-2.0*c)*(1.0+e)*(1.0-e)*(1.0+s)*s/4.0
        # Mid-edge ζ = 0
        self.grad_fc[0,16,:] = +(1.0-2.0*c)*(1.0-e)*e*(1.0+s)*(1.0-s)/4.0
        self.grad_fc[0,17,:] = -(1.0+2.0*c)*(1.0-e)*e*(1.0+s)*(1.0-s)/4.0
        self.grad_fc[0,18,:] = +(1.0+2.0*c)*(1.0+e)*e*(1.0+s)*(1.0-s)/4.0
        self.grad_fc[0,19,:] = -(1.0-2.0*c)*(1.0+e)*e*(1.0+s)*(1.0-s)/4.0
        # Mid-face ξ = +-1
        self.grad_fc[0,21,:] = +(1.0+2.0*c)*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)/2.0
        # self.grad_fc[0,22,:] = +(1.0+2.0*c)*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)/2.0  # Salome, ξ = +1
        self.grad_fc[0,23,:] = -(1.0-2.0*c)*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)/2.0
        # self.grad_fc[0,24,:] = -(1.0-2.0*c)*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)/2.0  # Salome, ξ = -1
        # Mid-face η = +-1
        self.grad_fc[0,20,:] = -(-2.0*c)*(1.0-e)*e*(1.0+s)*(1.0-s)/2.0
        # self.grad_fc[0,21,:] = -(-2.0*c)*(1.0-e)*e*(1.0+s)*(1.0-s)/2.0  # Salome, η = -1
        self.grad_fc[0,22,:] = +(-2.0*c)*(1.0+e)*e*(1.0+s)*(1.0-s)/2.0
        # self.grad_fc[0,23,:] = +(-2.0*c)*(1.0+e)*e*(1.0+s)*(1.0-s)/2.0  # Salome, η = +1
        # Mid-face ζ = +-1
        self.grad_fc[0,24,:] = -(-2.0*c)*(1.0+e)*(1.0-e)*(1.0-s)*s/2.0
        # self.grad_fc[0,20,:] = -(-2.0*c)*(1.0+e)*(1.0-e)*(1.0-s)*s/2.0  # Salome, ζ = -1
        self.grad_fc[0,25,:] = +(-2.0*c)*(1.0+e)*(1.0-e)*(1.0+s)*s/2.0
        # Central
        self.grad_fc[0,26,:] = +(-2.0*c)*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)

        # d2fdc2
        # Corners
        self.hess_fc[0,0,0,:]  = -(-2.0)*(1.0-e)*e*(1.0-s)*s/8.0
        self.hess_fc[0,0,1,:]  = +(+2.0)*(1.0-e)*e*(1.0-s)*s/8.0
        self.hess_fc[0,0,2,:]  = -(+2.0)*(1.0+e)*e*(1.0-s)*s/8.0
        self.hess_fc[0,0,3,:]  = +(-2.0)*(1.0+e)*e*(1.0-s)*s/8.0
        self.hess_fc[0,0,4,:]  = +(-2.0)*(1.0-e)*e*(1.0+s)*s/8.0
        self.hess_fc[0,0,5,:]  = -(+2.0)*(1.0-e)*e*(1.0+s)*s/8.0
        self.hess_fc[0,0,6,:]  = +(+2.0)*(1.0+e)*e*(1.0+s)*s/8.0
        self.hess_fc[0,0,7,:]  = -(-2.0)*(1.0+e)*e*(1.0+s)*s/8.0
        # Mid-edge ξ = 0
        self.hess_fc[0,0,8,:]  = +(-2.0)*(1.0-e)*e*(1.0-s)*s/4.0
        self.hess_fc[0,0,10,:] = -(-2.0)*(1.0+e)*e*(1.0-s)*s/4.0
        self.hess_fc[0,0,12,:] = -(-2.0)*(1.0-e)*e*(1.0+s)*s/4.0
        self.hess_fc[0,0,14,:] = +(-2.0)*(1.0+e)*e*(1.0+s)*s/4.0
        # Mid-edge η = 0
        self.hess_fc[0,0,9,:] = -(+2.0)*(1.0+e)*(1.0-e)*(1.0-s)*s/4.0
        self.hess_fc[0,0,11,:] = +(-2.0)*(1.0+e)*(1.0-e)*(1.0-s)*s/4.0
        self.hess_fc[0,0,13,:] = +(+2.0)*(1.0+e)*(1.0-e)*(1.0+s)*s/4.0
        self.hess_fc[0,0,15,:] = -(-2.0)*(1.0+e)*(1.0-e)*(1.0+s)*s/4.0
        # Mid-edge ζ = 0
        self.hess_fc[0,0,16,:] = +(-2.0)*(1.0-e)*e*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[0,0,17,:] = -(+2.0)*(1.0-e)*e*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[0,0,18,:] = +(+2.0)*(1.0+e)*e*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[0,0,19,:] = -(-2.0)*(1.0+e)*e*(1.0+s)*(1.0-s)/4.0
        # Mid-face ξ = +-1
        self.hess_fc[0,0,21,:] = +(+2.0)*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)/2.0
        # self.hess_fc[0,0,22,:] = +(+2.0)*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)/2.0  # Salome, ξ = +1
        self.hess_fc[0,0,23,:] = -(-2.0)*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)/2.0
        # self.hess_fc[0,0,24,:] = -(-2.0)*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)/2.0  # Salome, ξ = -1
        # Mid-face η = +-1
        self.hess_fc[0,0,20,:] = -(-2.0)*(1.0-e)*e*(1.0+s)*(1.0-s)/2.0
        # self.hess_fc[0,0,21,:] = -(-2.0)*(1.0-e)*e*(1.0+s)*(1.0-s)/2.0  # Salome, η = -1
        self.hess_fc[0,0,22,:] = +(-2.0)*(1.0+e)*e*(1.0+s)*(1.0-s)/2.0
        # self.hess_fc[0,0,23,:] = +(-2.0)*(1.0+e)*e*(1.0+s)*(1.0-s)/2.0  # Salome, η = +1
        # Mid-face ζ = +-1
        self.hess_fc[0,0,24,:] = -(-2.0)*(1.0+e)*(1.0-e)*(1.0-s)*s/2.0
        # self.hess_fc[0,0,20,:] = -(-2.0)*(1.0+e)*(1.0-e)*(1.0-s)*s/2.0  # Salome, ζ = -1
        self.hess_fc[0,0,25,:] = +(-2.0)*(1.0+e)*(1.0-e)*(1.0+s)*s/2.0
        # Central
        self.hess_fc[0,0,26,:] = +(-2.0)*(1.0+e)*(1.0-e)*(1.0+s)*(1.0-s)

        # d2fdce
        # Corners
        self.hess_fc[0,1,0,:]  = -(1.0-2.0*c)*(1.0-2.0*e)*(1.0-s)*s/8.0
        self.hess_fc[0,1,1,:]  = +(1.0+2.0*c)*(1.0-2.0*e)*(1.0-s)*s/8.0
        self.hess_fc[0,1,2,:]  = -(1.0+2.0*c)*(1.0+2.0*e)*(1.0-s)*s/8.0
        self.hess_fc[0,1,3,:]  = +(1.0-2.0*c)*(1.0+2.0*e)*(1.0-s)*s/8.0
        self.hess_fc[0,1,4,:]  = +(1.0-2.0*c)*(1.0-2.0*e)*(1.0+s)*s/8.0
        self.hess_fc[0,1,5,:]  = -(1.0+2.0*c)*(1.0-2.0*e)*(1.0+s)*s/8.0
        self.hess_fc[0,1,6,:]  = +(1.0+2.0*c)*(1.0+2.0*e)*(1.0+s)*s/8.0
        self.hess_fc[0,1,7,:]  = -(1.0-2.0*c)*(1.0+2.0*e)*(1.0+s)*s/8.0
        # Mid-edge ξ = 0
        self.hess_fc[0,1,8,:]  = +(-2.0*c)*(1.0-2.0*e)*(1.0-s)*s/4.0
        self.hess_fc[0,1,10,:] = -(-2.0*c)*(1.0+2.0*e)*(1.0-s)*s/4.0
        self.hess_fc[0,1,12,:] = -(-2.0*c)*(1.0-2.0*e)*(1.0+s)*s/4.0
        self.hess_fc[0,1,14,:] = +(-2.0*c)*(1.0+2.0*e)*(1.0+s)*s/4.0
        # Mid-edge η = 0
        self.hess_fc[0,1,9,:] = -(1.0+2.0*c)*(-2.0*e)*(1.0-s)*s/4.0
        self.hess_fc[0,1,11,:] = +(1.0-2.0*c)*(-2.0*e)*(1.0-s)*s/4.0
        self.hess_fc[0,1,13,:] = +(1.0+2.0*c)*(-2.0*e)*(1.0+s)*s/4.0
        self.hess_fc[0,1,15,:] = -(1.0-2.0*c)*(-2.0*e)*(1.0+s)*s/4.0
        # Mid-edge ζ = 0
        self.hess_fc[0,1,16,:] = +(1.0-2.0*c)*(1.0-2.0*e)*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[0,1,17,:] = -(1.0+2.0*c)*(1.0-2.0*e)*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[0,1,18,:] = +(1.0+2.0*c)*(1.0+2.0*e)*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[0,1,19,:] = -(1.0-2.0*c)*(1.0+2.0*e)*(1.0+s)*(1.0-s)/4.0
        # Mid-face ξ = +-1
        self.hess_fc[0,1,21,:] = +(1.0+2.0*c)*(-2.0*e)*(1.0+s)*(1.0-s)/2.0
        # self.hess_fc[0,1,22,:] = +(1.0+2.0*c)*(-2.0*e)*(1.0+s)*(1.0-s)/2.0  # Salome, ξ = +1
        self.hess_fc[0,1,23,:] = -(1.0-2.0*c)*(-2.0*e)*(1.0+s)*(1.0-s)/2.0
        # self.hess_fc[0,1,24,:] = -(1.0-2.0*c)*(-2.0*e)*(1.0+s)*(1.0-s)/2.0  # Salome, ξ = -1
        # Mid-face η = +-1
        self.hess_fc[0,1,20,:] = -(-2.0*c)*(1.0-2.0*e)*(1.0+s)*(1.0-s)/2.0
        # self.hess_fc[0,1,21,:] = -(-2.0*c)*(1.0-2.0*e)*(1.0+s)*(1.0-s)/2.0  # Salome, η = -1
        self.hess_fc[0,1,22,:] = +(-2.0*c)*(1.0+2.0*e)*(1.0+s)*(1.0-s)/2.0
        # self.hess_fc[0,1,23,:] = +(-2.0*c)*(1.0+2.0*e)*(1.0+s)*(1.0-s)/2.0  # Salome, η = +1
        # Mid-face ζ = +-1
        self.hess_fc[0,1,24,:] = -(-2.0*c)*(-2.0*e)*(1.0-s)*s/2.0
        # self.hess_fc[0,1,20,:] = -(-2.0*c)*(-2.0*e)*(1.0-s)*s/2.0  # Salome, ζ = -1
        self.hess_fc[0,1,25,:] = +(-2.0*c)*(-2.0*e)*(1.0+s)*s/2.0
        # Central
        self.hess_fc[0,1,26,:] = +(-2.0*c)*(-2.0*e)*(1.0+s)*(1.0-s)
        self.hess_fc[1,0,:,:] = self.hess_fc[0,1,:,:]

        # d2fdcs
        # Corners
        self.hess_fc[0,2,0,:]  = -(1.0-2.0*c)*(1.0-e)*e*(1.0-2.0*s)/8.0
        self.hess_fc[0,2,1,:]  = +(1.0+2.0*c)*(1.0-e)*e*(1.0-2.0*s)/8.0
        self.hess_fc[0,2,2,:]  = -(1.0+2.0*c)*(1.0+e)*e*(1.0-2.0*s)/8.0
        self.hess_fc[0,2,3,:]  = +(1.0-2.0*c)*(1.0+e)*e*(1.0-2.0*s)/8.0
        self.hess_fc[0,2,4,:]  = +(1.0-2.0*c)*(1.0-e)*e*(1.0+2.0*s)/8.0
        self.hess_fc[0,2,5,:]  = -(1.0+2.0*c)*(1.0-e)*e*(1.0+2.0*s)/8.0
        self.hess_fc[0,2,6,:]  = +(1.0+2.0*c)*(1.0+e)*e*(1.0+2.0*s)/8.0
        self.hess_fc[0,2,7,:]  = -(1.0-2.0*c)*(1.0+e)*e*(1.0+2.0*s)/8.0
        # Mid-edge ξ = 0
        self.hess_fc[0,2,8,:]  = +(-2.0*c)*(1.0-e)*e*(1.0-2.0*s)/4.0
        self.hess_fc[0,2,10,:] = -(-2.0*c)*(1.0+e)*e*(1.0-2.0*s)/4.0
        self.hess_fc[0,2,12,:] = -(-2.0*c)*(1.0-e)*e*(1.0+2.0*s)/4.0
        self.hess_fc[0,2,14,:] = +(-2.0*c)*(1.0+e)*e*(1.0+2.0*s)/4.0
        # Mid-edge η = 0
        self.hess_fc[0,2,9,:] = -(1.0+2.0*c)*(1.0+e)*(1.0-e)*(1.0-2.0*s)/4.0
        self.hess_fc[0,2,11,:] = +(1.0-2.0*c)*(1.0+e)*(1.0-e)*(1.0-2.0*s)/4.0
        self.hess_fc[0,2,13,:] = +(1.0+2.0*c)*(1.0+e)*(1.0-e)*(1.0+2.0*s)/4.0
        self.hess_fc[0,2,15,:] = -(1.0-2.0*c)*(1.0+e)*(1.0-e)*(1.0+2.0*s)/4.0
        # Mid-edge ζ = 0
        self.hess_fc[0,2,16,:] = +(1.0-2.0*c)*(1.0-e)*e*(-2.0*s)/4.0
        self.hess_fc[0,2,17,:] = -(1.0+2.0*c)*(1.0-e)*e*(-2.0*s)/4.0
        self.hess_fc[0,2,18,:] = +(1.0+2.0*c)*(1.0+e)*e*(-2.0*s)/4.0
        self.hess_fc[0,2,19,:] = -(1.0-2.0*c)*(1.0+e)*e*(-2.0*s)/4.0
        # Mid-face ξ = +-1
        self.hess_fc[0,2,21,:] = +(1.0+2.0*c)*(1.0+e)*(1.0-e)*(-2.0*s)/2.0
        # self.hess_fc[0,2,22,:] = +(1.0+2.0*c)*(1.0+e)*(1.0-e)*(-2.0*s)/2.0  # Salome, ξ = +1
        self.hess_fc[0,2,23,:] = -(1.0-2.0*c)*(1.0+e)*(1.0-e)*(-2.0*s)/2.0
        # self.hess_fc[0,2,24,:] = -(1.0-2.0*c)*(1.0+e)*(1.0-e)*(-2.0*s)/2.0  # Salome, ξ = -1
        # Mid-face η = +-1
        self.hess_fc[0,2,20,:] = -(-2.0*c)*(1.0-e)*e*(-2.0*s)/2.0
        # self.hess_fc[0,2,21,:] = -(-2.0*c)*(1.0-e)*e*(-2.0*s)/2.0  # Salome, η = -1
        self.hess_fc[0,2,22,:] = +(-2.0*c)*(1.0+e)*e*(-2.0*s)/2.0
        # self.hess_fc[0,2,23,:] = +(-2.0*c)*(1.0+e)*e*(-2.0*s)/2.0  # Salome, η = +1
        # Mid-face ζ = +-1
        self.hess_fc[0,2,24,:] = -(-2.0*c)*(1.0+e)*(1.0-e)*(1.0-2.0*s)/2.0
        # self.hess_fc[0,2,20,:] = -(-2.0*c)*(1.0+e)*(1.0-e)*(1.0-2.0*s)/2.0  # Salome, ζ = -1
        self.hess_fc[0,2,25,:] = +(-2.0*c)*(1.0+e)*(1.0-e)*(1.0+2.0*s)/2.0
        # Central
        self.hess_fc[0,2,26,:] = +(-2.0*c)*(1.0+e)*(1.0-e)*(-2.0*s)
        self.hess_fc[2,0,:,:] = self.hess_fc[0,2,:,:]

        # dfde
        # Corners
        self.grad_fc[1,0,:]  = -(1.0-c)*c*(1.0-2.0*e)*(1.0-s)*s/8.0
        self.grad_fc[1,1,:]  = +(1.0+c)*c*(1.0-2.0*e)*(1.0-s)*s/8.0
        self.grad_fc[1,2,:]  = -(1.0+c)*c*(1.0+2.0*e)*(1.0-s)*s/8.0
        self.grad_fc[1,3,:]  = +(1.0-c)*c*(1.0+2.0*e)*(1.0-s)*s/8.0
        self.grad_fc[1,4,:]  = +(1.0-c)*c*(1.0-2.0*e)*(1.0+s)*s/8.0
        self.grad_fc[1,5,:]  = -(1.0+c)*c*(1.0-2.0*e)*(1.0+s)*s/8.0
        self.grad_fc[1,6,:]  = +(1.0+c)*c*(1.0+2.0*e)*(1.0+s)*s/8.0
        self.grad_fc[1,7,:]  = -(1.0-c)*c*(1.0+2.0*e)*(1.0+s)*s/8.0
        # Mid-edge ξ = 0
        self.grad_fc[1,8,:]  = +(1.0+c)*(1.0-c)*(1.0-2.0*e)*(1.0-s)*s/4.0
        self.grad_fc[1,10,:] = -(1.0+c)*(1.0-c)*(1.0+2.0*e)*(1.0-s)*s/4.0
        self.grad_fc[1,12,:] = -(1.0+c)*(1.0-c)*(1.0-2.0*e)*(1.0+s)*s/4.0
        self.grad_fc[1,14,:] = +(1.0+c)*(1.0-c)*(1.0+2.0*e)*(1.0+s)*s/4.0
        # Mid-edge η = 0
        self.grad_fc[1,9,:] = -(1.0+c)*c*(-2.0*e)*(1.0-s)*s/4.0
        self.grad_fc[1,11,:] = +(1.0-c)*c*(-2.0*e)*(1.0-s)*s/4.0
        self.grad_fc[1,13,:] = +(1.0+c)*c*(-2.0*e)*(1.0+s)*s/4.0
        self.grad_fc[1,15,:] = -(1.0-c)*c*(-2.0*e)*(1.0+s)*s/4.0
        # Mid-edge ζ = 0
        self.grad_fc[1,16,:] = +(1.0-c)*c*(1.0-2.0*e)*(1.0+s)*(1.0-s)/4.0
        self.grad_fc[1,17,:] = -(1.0+c)*c*(1.0-2.0*e)*(1.0+s)*(1.0-s)/4.0
        self.grad_fc[1,18,:] = +(1.0+c)*c*(1.0+2.0*e)*(1.0+s)*(1.0-s)/4.0
        self.grad_fc[1,19,:] = -(1.0-c)*c*(1.0+2.0*e)*(1.0+s)*(1.0-s)/4.0
        # Mid-face ξ = +-1
        self.grad_fc[1,21,:] = +(1.0+c)*c*(-2.0*e)*(1.0+s)*(1.0-s)/2.0
        # self.grad_fc[1,22,:] = +(1.0+c)*c*(-2.0*e)*(1.0+s)*(1.0-s)/2.0  # Salome, ξ = +1
        self.grad_fc[1,23,:] = -(1.0-c)*c*(-2.0*e)*(1.0+s)*(1.0-s)/2.0
        # self.grad_fc[1,24,:] = -(1.0-c)*c*(-2.0*e)*(1.0+s)*(1.0-s)/2.0  # Salome, ξ = -1
        # Mid-face η = +-1
        self.grad_fc[1,20,:] = -(1.0+c)*(1.0-c)*(1.0-2.0*e)*(1.0+s)*(1.0-s)/2.0
        # self.grad_fc[1,21,:] = -(1.0+c)*(1.0-c)*(1.0-2.0*e)*(1.0+s)*(1.0-s)/2.0  # Salome, η = -1
        self.grad_fc[1,22,:] = +(1.0+c)*(1.0-c)*(1.0+2.0*e)*(1.0+s)*(1.0-s)/2.0
        # self.grad_fc[1,23,:] = +(1.0+c)*(1.0-c)*(1.0+2.0*e)*(1.0+s)*(1.0-s)/2.0  # Salome, η = +1
        # Mid-face ζ = +-1
        self.grad_fc[1,24,:] = -(1.0+c)*(1.0-c)*(-2.0*e)*(1.0-s)*s/2.0
        # self.grad_fc[1,20,:] = -(1.0+c)*(1.0-c)*(-2.0*e)*(1.0-s)*s/2.0  # Salome, ζ = -1
        self.grad_fc[1,25,:] = +(1.0+c)*(1.0-c)*(-2.0*e)*(1.0+s)*s/2.0
        # Central
        self.grad_fc[1,26,:] = +(1.0+c)*(1.0-c)*(-2.0*e)*(1.0+s)*(1.0-s)

        # d2fdce2
        # Corners
        self.hess_fc[1,1,0,:]  = -(1.0-c)*c*(-2.0)*(1.0-s)*s/8.0
        self.hess_fc[1,1,1,:]  = +(1.0+c)*c*(-2.0)*(1.0-s)*s/8.0
        self.hess_fc[1,1,2,:]  = -(1.0+c)*c*(+2.0)*(1.0-s)*s/8.0
        self.hess_fc[1,1,3,:]  = +(1.0-c)*c*(+2.0)*(1.0-s)*s/8.0
        self.hess_fc[1,1,4,:]  = +(1.0-c)*c*(-2.0)*(1.0+s)*s/8.0
        self.hess_fc[1,1,5,:]  = -(1.0+c)*c*(-2.0)*(1.0+s)*s/8.0
        self.hess_fc[1,1,6,:]  = +(1.0+c)*c*(+2.0)*(1.0+s)*s/8.0
        self.hess_fc[1,1,7,:]  = -(1.0-c)*c*(+2.0)*(1.0+s)*s/8.0
        # Mid-edge ξ = 0
        self.hess_fc[1,1,8,:]  = +(1.0+c)*(1.0-c)*(-2.0)*(1.0-s)*s/4.0
        self.hess_fc[1,1,10,:] = -(1.0+c)*(1.0-c)*(+2.0)*(1.0-s)*s/4.0
        self.hess_fc[1,1,12,:] = -(1.0+c)*(1.0-c)*(-2.0)*(1.0+s)*s/4.0
        self.hess_fc[1,1,14,:] = +(1.0+c)*(1.0-c)*(+2.0)*(1.0+s)*s/4.0
        # Mid-edge η = 0
        self.hess_fc[1,1,9,:] = -(1.0+c)*c*(-2.0)*(1.0-s)*s/4.0
        self.hess_fc[1,1,11,:] = +(1.0-c)*c*(-2.0)*(1.0-s)*s/4.0
        self.hess_fc[1,1,13,:] = +(1.0+c)*c*(-2.0)*(1.0+s)*s/4.0
        self.hess_fc[1,1,15,:] = -(1.0-c)*c*(-2.0)*(1.0+s)*s/4.0
        # Mid-edge ζ = 0
        self.hess_fc[1,1,16,:] = +(1.0-c)*c*(-2.0)*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[1,1,17,:] = -(1.0+c)*c*(-2.0)*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[1,1,18,:] = +(1.0+c)*c*(+2.0)*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[1,1,19,:] = -(1.0-c)*c*(+2.0)*(1.0+s)*(1.0-s)/4.0
        # Mid-face ξ = +-1
        self.hess_fc[1,1,21,:] = +(1.0+c)*c*(-2.0)*(1.0+s)*(1.0-s)/2.0
        # self.hess_fc[1,1,22,:] = +(1.0+c)*c*(-2.0)*(1.0+s)*(1.0-s)/2.0  # Salome, ξ = +1
        self.hess_fc[1,1,23,:] = -(1.0-c)*c*(-2.0)*(1.0+s)*(1.0-s)/2.0
        # self.hess_fc[1,1,24,:] = -(1.0-c)*c*(-2.0)*(1.0+s)*(1.0-s)/2.0  # Salome, ξ = -1
        # Mid-face η = +-1
        self.hess_fc[1,1,20,:] = -(1.0+c)*(1.0-c)*(-2.0)*(1.0+s)*(1.0-s)/2.0
        # self.hess_fc[1,1,21,:] = -(1.0+c)*(1.0-c)*(-2.0)*(1.0+s)*(1.0-s)/2.0  # Salome, η = -1
        self.hess_fc[1,1,22,:] = +(1.0+c)*(1.0-c)*(+2.0)*(1.0+s)*(1.0-s)/2.0
        # self.hess_fc[1,1,23,:] = +(1.0+c)*(1.0-c)*(+2.0)*(1.0+s)*(1.0-s)/2.0  # Salome, η = +1
        # Mid-face ζ = +-1
        self.hess_fc[1,1,24,:] = -(1.0+c)*(1.0-c)*(-2.0)*(1.0-s)*s/2.0
        # self.hess_fc[1,1,20,:] = -(1.0+c)*(1.0-c)*(-2.0)*(1.0-s)*s/2.0  # Salome, ζ = -1
        self.hess_fc[1,1,25,:] = +(1.0+c)*(1.0-c)*(-2.0)*(1.0+s)*s/2.0
        # Central
        self.hess_fc[1,1,26,:] = +(1.0+c)*(1.0-c)*(-2.0)*(1.0+s)*(1.0-s)

        # d2fdes
        # Corners
        self.hess_fc[1,2,0,:]  = -(1.0-c)*c*(1.0-2.0*e)*(1.0-2.0*s)/8.0
        self.hess_fc[1,2,1,:]  = +(1.0+c)*c*(1.0-2.0*e)*(1.0-2.0*s)/8.0
        self.hess_fc[1,2,2,:]  = -(1.0+c)*c*(1.0+2.0*e)*(1.0-2.0*s)/8.0
        self.hess_fc[1,2,3,:]  = +(1.0-c)*c*(1.0+2.0*e)*(1.0-2.0*s)/8.0
        self.hess_fc[1,2,4,:]  = +(1.0-c)*c*(1.0-2.0*e)*(1.0+2.0*s)/8.0
        self.hess_fc[1,2,5,:]  = -(1.0+c)*c*(1.0-2.0*e)*(1.0+2.0*s)/8.0
        self.hess_fc[1,2,6,:]  = +(1.0+c)*c*(1.0+2.0*e)*(1.0+2.0*s)/8.0
        self.hess_fc[1,2,7,:]  = -(1.0-c)*c*(1.0+2.0*e)*(1.0+2.0*s)/8.0
        # Mid-edge ξ = 0
        self.hess_fc[1,2,8,:]  = +(1.0+c)*(1.0-c)*(1.0-2.0*e)*(1.0-2.0*s)/4.0
        self.hess_fc[1,2,10,:] = -(1.0+c)*(1.0-c)*(1.0+2.0*e)*(1.0-2.0*s)/4.0
        self.hess_fc[1,2,12,:] = -(1.0+c)*(1.0-c)*(1.0-2.0*e)*(1.0+2.0*s)/4.0
        self.hess_fc[1,2,14,:] = +(1.0+c)*(1.0-c)*(1.0+2.0*e)*(1.0+2.0*s)/4.0
        # Mid-edge η = 0
        self.hess_fc[1,2,9,:] = -(1.0+c)*c*(-2.0*e)*(1.0-2.0*s)/4.0
        self.hess_fc[1,2,11,:] = +(1.0-c)*c*(-2.0*e)*(1.0-2.0*s)/4.0
        self.hess_fc[1,2,13,:] = +(1.0+c)*c*(-2.0*e)*(1.0+2.0*s)/4.0
        self.hess_fc[1,2,15,:] = -(1.0-c)*c*(-2.0*e)*(1.0+2.0*s)/4.0
        # Mid-edge ζ = 0
        self.hess_fc[1,2,16,:] = +(1.0-c)*c*(1.0-2.0*e)*(-2.0*s)/4.0
        self.hess_fc[1,2,17,:] = -(1.0+c)*c*(1.0-2.0*e)*(-2.0*s)/4.0
        self.hess_fc[1,2,18,:] = +(1.0+c)*c*(1.0+2.0*e)*(-2.0*s)/4.0
        self.hess_fc[1,2,19,:] = -(1.0-c)*c*(1.0+2.0*e)*(-2.0*s)/4.0
        # Mid-face ξ = +-1
        self.hess_fc[1,2,21,:] = +(1.0+c)*c*(-2.0*e)*(-2.0*s)/2.0
        # self.hess_fc[1,2,22,:] = +(1.0+c)*c*(-2.0*e)*(-2.0*s)/2.0  # Salome, ξ = +1
        self.hess_fc[1,2,23,:] = -(1.0-c)*c*(-2.0*e)*(-2.0*s)/2.0
        # self.hess_fc[1,2,24,:] = -(1.0-c)*c*(-2.0*e)*(-2.0*s)/2.0  # Salome, ξ = -1
        # Mid-face η = +-1
        self.hess_fc[1,2,20,:] = -(1.0+c)*(1.0-c)*(1.0-2.0*e)*(-2.0*s)/2.0
        # self.hess_fc[1,2,21,:] = -(1.0+c)*(1.0-c)*(1.0-2.0*e)*(-2.0*s)/2.0  # Salome, η = -1
        self.hess_fc[1,2,22,:] = +(1.0+c)*(1.0-c)*(1.0+2.0*e)*(-2.0*s)/2.0
        # self.hess_fc[1,2,23,:] = +(1.0+c)*(1.0-c)*(1.0+2.0*e)*(-2.0*s)/2.0  # Salome, η = +1
        # Mid-face ζ = +-1
        self.hess_fc[1,2,24,:] = -(1.0+c)*(1.0-c)*(-2.0*e)*(1.0-2.0*s)/2.0
        # self.hess_fc[1,2,20,:] = -(1.0+c)*(1.0-c)*(-2.0*e)*(1.0-2.0*s)/2.0  # Salome, ζ = -1
        self.hess_fc[1,2,25,:] = +(1.0+c)*(1.0-c)*(-2.0*e)*(1.0+2.0*s)/2.0
        # Central
        self.hess_fc[1,2,26,:] = +(1.0+c)*(1.0-c)*(-2.0*e)*(-2.0*s)
        self.hess_fc[2,1,:,:] = self.hess_fc[1,2,:,:]

        # dfds
        # Corners
        self.grad_fc[2,0,:]  = -(1.0-c)*c*(1.0-e)*e*(1.0-2.0*s)/8.0
        self.grad_fc[2,1,:]  = +(1.0+c)*c*(1.0-e)*e*(1.0-2.0*s)/8.0
        self.grad_fc[2,2,:]  = -(1.0+c)*c*(1.0+e)*e*(1.0-2.0*s)/8.0
        self.grad_fc[2,3,:]  = +(1.0-c)*c*(1.0+e)*e*(1.0-2.0*s)/8.0
        self.grad_fc[2,4,:]  = +(1.0-c)*c*(1.0-e)*e*(1.0+2.0*s)/8.0
        self.grad_fc[2,5,:]  = -(1.0+c)*c*(1.0-e)*e*(1.0+2.0*s)/8.0
        self.grad_fc[2,6,:]  = +(1.0+c)*c*(1.0+e)*e*(1.0+2.0*s)/8.0
        self.grad_fc[2,7,:]  = -(1.0-c)*c*(1.0+e)*e*(1.0+2.0*s)/8.0
        # Mid-edge ξ = 0
        self.grad_fc[2,8,:]  = +(1.0+c)*(1.0-c)*(1.0-e)*e*(1.0-2.0*s)/4.0
        self.grad_fc[2,10,:] = -(1.0+c)*(1.0-c)*(1.0+e)*e*(1.0-2.0*s)/4.0
        self.grad_fc[2,12,:] = -(1.0+c)*(1.0-c)*(1.0-e)*e*(1.0+2.0*s)/4.0
        self.grad_fc[2,14,:] = +(1.0+c)*(1.0-c)*(1.0+e)*e*(1.0+2.0*s)/4.0
        # Mid-edge η = 0
        self.grad_fc[2,9,:] = -(1.0+c)*c*(1.0+e)*(1.0-e)*(1.0-2.0*s)/4.0
        self.grad_fc[2,11,:] = +(1.0-c)*c*(1.0+e)*(1.0-e)*(1.0-2.0*s)/4.0
        self.grad_fc[2,13,:] = +(1.0+c)*c*(1.0+e)*(1.0-e)*(1.0+2.0*s)/4.0
        self.grad_fc[2,15,:] = -(1.0-c)*c*(1.0+e)*(1.0-e)*(1.0+2.0*s)/4.0
        # Mid-edge ζ = 0
        self.grad_fc[2,16,:] = +(1.0-c)*c*(1.0-e)*e*(-2.0*s)/4.0
        self.grad_fc[2,17,:] = -(1.0+c)*c*(1.0-e)*e*(-2.0*s)/4.0
        self.grad_fc[2,18,:] = +(1.0+c)*c*(1.0+e)*e*(-2.0*s)/4.0
        self.grad_fc[2,19,:] = -(1.0-c)*c*(1.0+e)*e*(-2.0*s)/4.0
        # Mid-face ξ = +-1
        self.grad_fc[2,21,:] = +(1.0+c)*c*(1.0+e)*(1.0-e)*(-2.0*s)/2.0
        # self.grad_fc[2,22,:] = +(1.0+c)*c*(1.0+e)*(1.0-e)*(-2.0*s)/2.0  # Salome, ξ = +1
        self.grad_fc[2,23,:] = -(1.0-c)*c*(1.0+e)*(1.0-e)*(-2.0*s)/2.0
        # self.grad_fc[2,24,:] = -(1.0-c)*c*(1.0+e)*(1.0-e)*(-2.0*s)/2.0  # Salome, ξ = -1
        # Mid-face η = +-1
        self.grad_fc[2,20,:] = -(1.0+c)*(1.0-c)*(1.0-e)*e*(-2.0*s)/2.0
        # self.grad_fc[2,21,:] = -(1.0+c)*(1.0-c)*(1.0-e)*e*(-2.0*s)/2.0  # Salome, η = -1
        self.grad_fc[2,22,:] = +(1.0+c)*(1.0-c)*(1.0+e)*e*(-2.0*s)/2.0
        # self.grad_fc[2,23,:] = +(1.0+c)*(1.0-c)*(1.0+e)*e*(-2.0*s)/2.0  # Salome, η = +1
        # Mid-face ζ = +-1
        self.grad_fc[2,24,:] = -(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)*(1.0-2.0*s)/2.0
        # self.grad_fc[2,20,:] = -(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)*(1.0-2.0*s)/2.0  # Salome, ζ = -1
        self.grad_fc[2,25,:] = +(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)*(1.0+2.0*s)/2.0
        # Central
        self.grad_fc[2,26,:] = +(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)*(-2.0*s)

        # d2fds2
        # Corners
        self.hess_fc[2,2,0,:]  = -(1.0-c)*c*(1.0-e)*e*(-2.0)/8.0
        self.hess_fc[2,2,1,:]  = +(1.0+c)*c*(1.0-e)*e*(-2.0)/8.0
        self.hess_fc[2,2,2,:]  = -(1.0+c)*c*(1.0+e)*e*(-2.0)/8.0
        self.hess_fc[2,2,3,:]  = +(1.0-c)*c*(1.0+e)*e*(-2.0)/8.0
        self.hess_fc[2,2,4,:]  = +(1.0-c)*c*(1.0-e)*e*(+2.0)/8.0
        self.hess_fc[2,2,5,:]  = -(1.0+c)*c*(1.0-e)*e*(+2.0)/8.0
        self.hess_fc[2,2,6,:]  = +(1.0+c)*c*(1.0+e)*e*(+2.0)/8.0
        self.hess_fc[2,2,7,:]  = -(1.0-c)*c*(1.0+e)*e*(+2.0)/8.0
        # Mid-edge ξ = 0
        self.hess_fc[2,2,8,:]  = +(1.0+c)*(1.0-c)*(1.0-e)*e*(-2.0)/4.0
        self.hess_fc[2,2,10,:] = -(1.0+c)*(1.0-c)*(1.0+e)*e*(-2.0)/4.0
        self.hess_fc[2,2,12,:] = -(1.0+c)*(1.0-c)*(1.0-e)*e*(+2.0)/4.0
        self.hess_fc[2,2,14,:] = +(1.0+c)*(1.0-c)*(1.0+e)*e*(+2.0)/4.0
        # Mid-edge η = 0
        self.hess_fc[2,2,9,:] = -(1.0+c)*c*(1.0+e)*(1.0-e)*(-2.0)/4.0
        self.hess_fc[2,2,11,:] = +(1.0-c)*c*(1.0+e)*(1.0-e)*(-2.0)/4.0
        self.hess_fc[2,2,13,:] = +(1.0+c)*c*(1.0+e)*(1.0-e)*(+2.0)/4.0
        self.hess_fc[2,2,15,:] = -(1.0-c)*c*(1.0+e)*(1.0-e)*(+2.0)/4.0
        # Mid-edge ζ = 0
        self.hess_fc[2,2,16,:] = +(1.0-c)*c*(1.0-e)*e*(-2.0)/4.0
        self.hess_fc[2,2,17,:] = -(1.0+c)*c*(1.0-e)*e*(-2.0)/4.0
        self.hess_fc[2,2,18,:] = +(1.0+c)*c*(1.0+e)*e*(-2.0)/4.0
        self.hess_fc[2,2,19,:] = -(1.0-c)*c*(1.0+e)*e*(-2.0)/4.0
        # Mid-face ξ = +-1
        self.hess_fc[2,2,21,:] = +(1.0+c)*c*(1.0+e)*(1.0-e)*(-2.0)/2.0
        # self.hess_fc[2,2,22,:] = +(1.0+c)*c*(1.0+e)*(1.0-e)*(-2.0)/2.0  # Salome, ξ = +1
        self.hess_fc[2,2,23,:] = -(1.0-c)*c*(1.0+e)*(1.0-e)*(-2.0)/2.0
        # self.hess_fc[2,2,24,:] = -(1.0-c)*c*(1.0+e)*(1.0-e)*(-2.0)/2.0  # Salome, ξ = -1
        # Mid-face η = +-1
        self.hess_fc[2,2,20,:] = -(1.0+c)*(1.0-c)*(1.0-e)*e*(-2.0)/2.0
        # self.hess_fc[2,2,21,:] = -(1.0+c)*(1.0-c)*(1.0-e)*e*(-2.0)/2.0  # Salome, η = -1
        self.hess_fc[2,2,22,:] = +(1.0+c)*(1.0-c)*(1.0+e)*e*(-2.0)/2.0
        # self.hess_fc[2,2,23,:] = +(1.0+c)*(1.0-c)*(1.0+e)*e*(-2.0)/2.0  # Salome, η = +1
        # Mid-face ζ = +-1
        self.hess_fc[2,2,24,:] = -(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)*(-2.0)/2.0
        # self.hess_fc[2,2,20,:] = -(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)*(-2.0)/2.0  # Salome, ζ = -1
        self.hess_fc[2,2,25,:] = +(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)*(+2.0)/2.0
        # Central
        self.hess_fc[2,2,26,:] = +(1.0+c)*(1.0-c)*(1.0+e)*(1.0-e)*(-2.0)

    def set_serendipity_basis_functions(self):

        c = self.c[:,0]
        e = self.c[:,1]
        s = self.c[:,2]

        # Corners
        self.bfn[0,:]  = -(1.0-c)*(1.0-e)*(1.0-s)*(2.0+c+e+s)/8.0
        self.bfn[1,:]  = -(1.0+c)*(1.0-e)*(1.0-s)*(2.0-c+e+s)/8.0
        self.bfn[2,:]  = -(1.0+c)*(1.0+e)*(1.0-s)*(2.0-c-e+s)/8.0
        self.bfn[3,:]  = -(1.0-c)*(1.0+e)*(1.0-s)*(2.0+c-e+s)/8.0
        self.bfn[4,:]  = -(1.0-c)*(1.0-e)*(1.0+s)*(2.0+c+e-s)/8.0
        self.bfn[5,:]  = -(1.0+c)*(1.0-e)*(1.0+s)*(2.0-c+e-s)/8.0
        self.bfn[6,:]  = -(1.0+c)*(1.0+e)*(1.0+s)*(2.0-c-e-s)/8.0
        self.bfn[7,:]  = -(1.0-c)*(1.0+e)*(1.0+s)*(2.0+c-e-s)/8.0
        # Mid-edge ξ = 0
        self.bfn[8,:]  = +(1.0+c)*(1.0-c)*(1.0-e)*(1.0-s)/4.0
        self.bfn[10,:] = +(1.0+c)*(1.0-c)*(1.0+e)*(1.0-s)/4.0
        self.bfn[12,:] = +(1.0+c)*(1.0-c)*(1.0-e)*(1.0+s)/4.0
        self.bfn[14,:] = +(1.0+c)*(1.0-c)*(1.0+e)*(1.0+s)/4.0
        # Mid-edge η = 0
        self.bfn[9,:] = +(1.0+c)*(1.0+e)*(1.0-e)*(1.0-s)/4.0
        self.bfn[11,:] = +(1.0-c)*(1.0+e)*(1.0-e)*(1.0-s)/4.0
        self.bfn[13,:] = +(1.0+c)*(1.0+e)*(1.0-e)*(1.0+s)/4.0
        self.bfn[15,:] = +(1.0-c)*(1.0+e)*(1.0-e)*(1.0+s)/4.0
        # Mid-edge ζ = 0
        self.bfn[16,:] = +(1.0-c)*(1.0-e)*(1.0+s)*(1.0-s)/4.0
        self.bfn[17,:] = +(1.0+c)*(1.0-e)*(1.0+s)*(1.0-s)/4.0
        self.bfn[18,:] = +(1.0+c)*(1.0+e)*(1.0+s)*(1.0-s)/4.0
        self.bfn[19,:] = +(1.0-c)*(1.0+e)*(1.0+s)*(1.0-s)/4.0

        # dfdc
        # Corners
        self.grad_fc[0,0,:]  = -((-1.0)*(2.0+c+e+s)+(1.0-c)*(+1.0))*(1.0-e)*(1.0-s)/8.0
        self.grad_fc[0,1,:]  = -((+1.0)*(2.0-c+e+s)+(1.0+c)*(-1.0))*(1.0-e)*(1.0-s)/8.0
        self.grad_fc[0,2,:]  = -((+1.0)*(2.0-c-e+s)+(1.0+c)*(-1.0))*(1.0+e)*(1.0-s)/8.0
        self.grad_fc[0,3,:]  = -((-1.0)*(2.0+c-e+s)+(1.0-c)*(+1.0))*(1.0+e)*(1.0-s)/8.0
        self.grad_fc[0,4,:]  = -((-1.0)*(2.0+c+e-s)+(1.0-c)*(+1.0))*(1.0-e)*(1.0+s)/8.0
        self.grad_fc[0,5,:]  = -((+1.0)*(2.0-c+e-s)+(1.0+c)*(-1.0))*(1.0-e)*(1.0+s)/8.0
        self.grad_fc[0,6,:]  = -((+1.0)*(2.0-c-e-s)+(1.0+c)*(-1.0))*(1.0+e)*(1.0+s)/8.0
        self.grad_fc[0,7,:]  = -((-1.0)*(2.0+c-e-s)+(1.0-c)*(+1.0))*(1.0+e)*(1.0+s)/8.0
        # Mid-edge ξ = 0
        self.grad_fc[0,8,:]  = +(-2.0*c)*(1.0-e)*(1.0-s)/4.0
        self.grad_fc[0,10,:] = +(-2.0*c)*(1.0+e)*(1.0-s)/4.0
        self.grad_fc[0,12,:] = +(-2.0*c)*(1.0-e)*(1.0+s)/4.0
        self.grad_fc[0,14,:] = +(-2.0*c)*(1.0+e)*(1.0+s)/4.0
        # Mid-edge η = 0
        self.grad_fc[0,9,:] = +(+1.0)*(1.0+e)*(1.0-e)*(1.0-s)/4.0
        self.grad_fc[0,11,:] = +(-1.0)*(1.0+e)*(1.0-e)*(1.0-s)/4.0
        self.grad_fc[0,13,:] = +(+1.0)*(1.0+e)*(1.0-e)*(1.0+s)/4.0
        self.grad_fc[0,15,:] = +(-1.0)*(1.0+e)*(1.0-e)*(1.0+s)/4.0
        # Mid-edge ζ = 0
        self.grad_fc[0,16,:] = +(-1.0)*(1.0-e)*(1.0+s)*(1.0-s)/4.0
        self.grad_fc[0,17,:] = +(+1.0)*(1.0-e)*(1.0+s)*(1.0-s)/4.0
        self.grad_fc[0,18,:] = +(+1.0)*(1.0+e)*(1.0+s)*(1.0-s)/4.0
        self.grad_fc[0,19,:] = +(-1.0)*(1.0+e)*(1.0+s)*(1.0-s)/4.0

        # d2fdc2
        # Corners
        self.hess_fc[0,0,0,:]  = -((-1.0)*(+1.0)+(-1.0)*(+1.0))*(1.0-e)*(1.0-s)/8.0
        self.hess_fc[0,0,1,:]  = -((+1.0)*(-1.0)+(+1.0)*(-1.0))*(1.0-e)*(1.0-s)/8.0
        self.hess_fc[0,0,2,:]  = -((+1.0)*(-1.0)+(+1.0)*(-1.0))*(1.0+e)*(1.0-s)/8.0
        self.hess_fc[0,0,3,:]  = -((-1.0)*(+1.0)+(-1.0)*(+1.0))*(1.0+e)*(1.0-s)/8.0
        self.hess_fc[0,0,4,:]  = -((-1.0)*(+1.0)+(-1.0)*(+1.0))*(1.0-e)*(1.0+s)/8.0
        self.hess_fc[0,0,5,:]  = -((+1.0)*(-1.0)+(+1.0)*(-1.0))*(1.0-e)*(1.0+s)/8.0
        self.hess_fc[0,0,6,:]  = -((+1.0)*(-1.0)+(+1.0)*(-1.0))*(1.0+e)*(1.0+s)/8.0
        self.hess_fc[0,0,7,:]  = -((-1.0)*(+1.0)+(-1.0)*(+1.0))*(1.0+e)*(1.0+s)/8.0
        # Mid-edge ξ = 0
        self.hess_fc[0,0,8,:]  = +(-2.0)*(1.0-e)*(1.0-s)/4.0
        self.hess_fc[0,0,10,:] = +(-2.0)*(1.0+e)*(1.0-s)/4.0
        self.hess_fc[0,0,12,:] = +(-2.0)*(1.0-e)*(1.0+s)/4.0
        self.hess_fc[0,0,14,:] = +(-2.0)*(1.0+e)*(1.0+s)/4.0
        # Mid-edge η = 0
        self.hess_fc[0,0,9,:] = +0.0
        self.hess_fc[0,0,11,:] = +0.0
        self.hess_fc[0,0,13,:] = +0.0
        self.hess_fc[0,0,15,:] = +0.0
        # Mid-edge ζ = 0
        self.hess_fc[0,0,16,:] = +0.0
        self.hess_fc[0,0,17,:] = +0.0
        self.hess_fc[0,0,18,:] = +0.0
        self.hess_fc[0,0,19,:] = +0.0

        # d2fdce
        # Corners
        self.hess_fc[0,1,0,:]  = -( ((-1.0)*(+1.0)+0.0)*(1.0-e) + ((-1.0)*(2.0+c+e+s)+(1.0-c)*(+1.0))*(-1.0) )*(1.0-s)/8.0
        self.hess_fc[0,1,1,:]  = -( ((+1.0)*(+1.0)+0.0)*(1.0-e) + ((+1.0)*(2.0-c+e+s)+(1.0+c)*(-1.0))*(-1.0) )*(1.0-s)/8.0
        self.hess_fc[0,1,2,:]  = -( ((+1.0)*(-1.0)+0.0)*(1.0+e) + ((+1.0)*(2.0-c-e+s)+(1.0+c)*(-1.0))*(+1.0) )*(1.0-s)/8.0
        self.hess_fc[0,1,3,:]  = -( ((-1.0)*(-1.0)+0.0)*(1.0+e) + ((-1.0)*(2.0+c-e+s)+(1.0-c)*(+1.0))*(+1.0) )*(1.0-s)/8.0
        self.hess_fc[0,1,4,:]  = -( ((-1.0)*(+1.0)+0.0)*(1.0-e) + ((-1.0)*(2.0+c+e-s)+(1.0-c)*(+1.0))*(-1.0) )*(1.0+s)/8.0
        self.hess_fc[0,1,5,:]  = -( ((+1.0)*(+1.0)+0.0)*(1.0-e) + ((+1.0)*(2.0-c+e-s)+(1.0+c)*(-1.0))*(-1.0) )*(1.0+s)/8.0
        self.hess_fc[0,1,6,:]  = -( ((+1.0)*(-1.0)+0.0)*(1.0+e) + ((+1.0)*(2.0-c-e-s)+(1.0+c)*(-1.0))*(+1.0) )*(1.0+s)/8.0
        self.hess_fc[0,1,7,:]  = -( ((-1.0)*(-1.0)+0.0)*(1.0+e) + ((-1.0)*(2.0+c-e-s)+(1.0-c)*(+1.0))*(+1.0) )*(1.0+s)/8.0
        # Mid-edge ξ = 0
        self.hess_fc[0,1,8,:]  = +(-2.0*c)*(-1.0)*(1.0-s)/4.0
        self.hess_fc[0,1,10,:] = +(-2.0*c)*(+1.0)*(1.0-s)/4.0
        self.hess_fc[0,1,12,:] = +(-2.0*c)*(-1.0)*(1.0+s)/4.0
        self.hess_fc[0,1,14,:] = +(-2.0*c)*(+1.0)*(1.0+s)/4.0
        # Mid-edge η = 0
        self.hess_fc[0,1,9,:] = +(+1.0)*(-2.0*e)*(1.0-s)/4.0
        self.hess_fc[0,1,11,:] = +(-1.0)*(-2.0*e)*(1.0-s)/4.0
        self.hess_fc[0,1,13,:] = +(+1.0)*(-2.0*e)*(1.0+s)/4.0
        self.hess_fc[0,1,15,:] = +(-1.0)*(-2.0*e)*(1.0+s)/4.0
        # Mid-edge ζ = 0
        self.hess_fc[0,1,16,:] = +(-1.0)*(-1.0)*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[0,1,17,:] = +(+1.0)*(-1.0)*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[0,1,18,:] = +(+1.0)*(+1.0)*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[0,1,19,:] = +(-1.0)*(+1.0)*(1.0+s)*(1.0-s)/4.0
        self.hess_fc[1,0,:,:] = self.hess_fc[0,1,:,:]

        # d2fdcs
        # Corners
        self.hess_fc[0,2,0,:] = -( ((-1.0)*(+1.0)+0.0)*(1.0-s)+((-1.0)*(2.0+c+e+s)+(1.0-c)*(+1.0))*(-1.0) )*(1.0-e)/8.0
        self.hess_fc[0,2,1,:] = -( ((+1.0)*(+1.0)+0.0)*(1.0-s)+((+1.0)*(2.0-c+e+s)+(1.0+c)*(-1.0))*(-1.0) )*(1.0-e)/8.0
        self.hess_fc[0,2,2,:] = -( ((+1.0)*(+1.0)+0.0)*(1.0-s)+((+1.0)*(2.0-c-e+s)+(1.0+c)*(-1.0))*(-1.0) )*(1.0+e)/8.0
        self.hess_fc[0,2,3,:] = -( ((-1.0)*(+1.0)+0.0)*(1.0-s)+((-1.0)*(2.0+c-e+s)+(1.0-c)*(+1.0))*(-1.0) )*(1.0+e)/8.0
        self.hess_fc[0,2,4,:] = -( ((-1.0)*(-1.0)+0.0)*(1.0+s)+((-1.0)*(2.0+c+e-s)+(1.0-c)*(+1.0))*(+1.0) )*(1.0-e)/8.0
        self.hess_fc[0,2,5,:] = -( ((+1.0)*(-1.0)+0.0)*(1.0+s)+((+1.0)*(2.0-c+e-s)+(1.0+c)*(-1.0))*(+1.0) )*(1.0-e)/8.0
        self.hess_fc[0,2,6,:] = -( ((+1.0)*(-1.0)+0.0)*(1.0+s)+((+1.0)*(2.0-c-e-s)+(1.0+c)*(-1.0))*(+1.0) )*(1.0+e)/8.0
        self.hess_fc[0,2,7,:] = -( ((-1.0)*(-1.0)+0.0)*(1.0+s)+((-1.0)*(2.0+c-e-s)+(1.0-c)*(+1.0))*(+1.0) )*(1.0+e)/8.0
        # Mid-edge ξ = 0
        self.hess_fc[0,2,8,:]  = +(-2.0*c)*(1.0-e)*(-1.0)/4.0
        self.hess_fc[0,2,10,:] = +(-2.0*c)*(1.0+e)*(-1.0)/4.0
        self.hess_fc[0,2,12,:] = +(-2.0*c)*(1.0-e)*(+1.0)/4.0
        self.hess_fc[0,2,14,:] = +(-2.0*c)*(1.0+e)*(+1.0)/4.0
        # Mid-edge η = 0
        self.hess_fc[0,2,9,:] = +(+1.0)*(1.0+e)*(1.0-e)*(-1.0)/4.0
        self.hess_fc[0,2,11,:] = +(-1.0)*(1.0+e)*(1.0-e)*(-1.0)/4.0
        self.hess_fc[0,2,13,:] = +(+1.0)*(1.0+e)*(1.0-e)*(+1.0)/4.0
        self.hess_fc[0,2,15,:] = +(-1.0)*(1.0+e)*(1.0-e)*(+1.0)/4.0
        # Mid-edge ζ = 0
        self.hess_fc[0,2,16,:] = +(-1.0)*(1.0-e)*(-2.0*s)/4.0
        self.hess_fc[0,2,17,:] = +(+1.0)*(1.0-e)*(-2.0*s)/4.0
        self.hess_fc[0,2,18,:] = +(+1.0)*(1.0+e)*(-2.0*s)/4.0
        self.hess_fc[0,2,19,:] = +(-1.0)*(1.0+e)*(-2.0*s)/4.0
        self.hess_fc[2,0,:,:] = self.hess_fc[0,2,:,:]

        # dfde
        # Corners
        self.grad_fc[1,0,:]  = -(1.0-c)*( (-1.0)*(2.0+c+e+s)+(1.0-e)*(+1.0) )*(1.0-s)/8.0
        self.grad_fc[1,1,:]  = -(1.0+c)*( (-1.0)*(2.0-c+e+s)+(1.0-e)*(+1.0) )*(1.0-s)/8.0
        self.grad_fc[1,2,:]  = -(1.0+c)*( (+1.0)*(2.0-c-e+s)+(1.0+e)*(-1.0) )*(1.0-s)/8.0
        self.grad_fc[1,3,:]  = -(1.0-c)*( (+1.0)*(2.0+c-e+s)+(1.0+e)*(-1.0) )*(1.0-s)/8.0
        self.grad_fc[1,4,:]  = -(1.0-c)*( (-1.0)*(2.0+c+e-s)+(1.0-e)*(+1.0) )*(1.0+s)/8.0
        self.grad_fc[1,5,:]  = -(1.0+c)*( (-1.0)*(2.0-c+e-s)+(1.0-e)*(+1.0) )*(1.0+s)/8.0
        self.grad_fc[1,6,:]  = -(1.0+c)*( (+1.0)*(2.0-c-e-s)+(1.0+e)*(-1.0) )*(1.0+s)/8.0
        self.grad_fc[1,7,:]  = -(1.0-c)*( (+1.0)*(2.0+c-e-s)+(1.0+e)*(-1.0) )*(1.0+s)/8.0
        # Mid-edge ξ = 0
        self.grad_fc[1,8,:]  = +(1.0+c)*(1.0-c)*(-1.0)*(1.0-s)/4.0
        self.grad_fc[1,10,:] = +(1.0+c)*(1.0-c)*(+1.0)*(1.0-s)/4.0
        self.grad_fc[1,12,:] = +(1.0+c)*(1.0-c)*(-1.0)*(1.0+s)/4.0
        self.grad_fc[1,14,:] = +(1.0+c)*(1.0-c)*(+1.0)*(1.0+s)/4.0
        # Mid-edge η = 0
        self.grad_fc[1,9,:] = +(1.0+c)*(-2.0*e)*(1.0-s)/4.0
        self.grad_fc[1,11,:] = +(1.0-c)*(-2.0*e)*(1.0-s)/4.0
        self.grad_fc[1,13,:] = +(1.0+c)*(-2.0*e)*(1.0+s)/4.0
        self.grad_fc[1,15,:] = +(1.0-c)*(-2.0*e)*(1.0+s)/4.0
        # Mid-edge ζ = 0
        self.grad_fc[1,16,:] = +(1.0-c)*(-1.0)*(1.0+s)*(1.0-s)/4.0
        self.grad_fc[1,17,:] = +(1.0+c)*(-1.0)*(1.0+s)*(1.0-s)/4.0
        self.grad_fc[1,18,:] = +(1.0+c)*(+1.0)*(1.0+s)*(1.0-s)/4.0
        self.grad_fc[1,19,:] = +(1.0-c)*(+1.0)*(1.0+s)*(1.0-s)/4.0

        # d2fde2
        # Corners
        self.hess_fc[1,1,0,:]  = -(1.0-c)*( (-1.0)*(+1.0)+(-1.0)*(+1.0) )*(1.0-s)/8.0
        self.hess_fc[1,1,1,:]  = -(1.0+c)*( (-1.0)*(+1.0)+(-1.0)*(+1.0) )*(1.0-s)/8.0
        self.hess_fc[1,1,2,:]  = -(1.0+c)*( (+1.0)*(-1.0)+(+1.0)*(-1.0) )*(1.0-s)/8.0
        self.hess_fc[1,1,3,:]  = -(1.0-c)*( (+1.0)*(-1.0)+(+1.0)*(-1.0) )*(1.0-s)/8.0
        self.hess_fc[1,1,4,:]  = -(1.0-c)*( (-1.0)*(+1.0)+(-1.0)*(+1.0) )*(1.0+s)/8.0
        self.hess_fc[1,1,5,:]  = -(1.0+c)*( (-1.0)*(+1.0)+(-1.0)*(+1.0) )*(1.0+s)/8.0
        self.hess_fc[1,1,6,:]  = -(1.0+c)*( (+1.0)*(-1.0)+(+1.0)*(-1.0) )*(1.0+s)/8.0
        self.hess_fc[1,1,7,:]  = -(1.0-c)*( (+1.0)*(-1.0)+(+1.0)*(-1.0) )*(1.0+s)/8.0
        # Mid-edge ξ = 0
        self.hess_fc[1,1,8,:]  = +0.0
        self.hess_fc[1,1,10,:] = +0.0
        self.hess_fc[1,1,12,:] = +0.0
        self.hess_fc[1,1,14,:] = +0.0
        # Mid-edge η = 0
        self.hess_fc[1,1,9,:] = +(1.0+c)*(-2.0)*(1.0-s)/4.0
        self.hess_fc[1,1,11,:] = +(1.0-c)*(-2.0)*(1.0-s)/4.0
        self.hess_fc[1,1,13,:] = +(1.0+c)*(-2.0)*(1.0+s)/4.0
        self.hess_fc[1,1,15,:] = +(1.0-c)*(-2.0)*(1.0+s)/4.0
        # Mid-edge ζ = 0
        self.hess_fc[1,1,16,:] = +0.0
        self.hess_fc[1,1,17,:] = +0.0
        self.hess_fc[1,1,18,:] = +0.0
        self.hess_fc[1,1,19,:] = +0.0

        # d2fdes
        # Corners
        self.hess_fc[1,2,0,:] = -(1.0-c)*( ((-1.0)*(+1.0)+0.0)*(1.0-s) + ((-1.0)*(2.0+c+e+s)+(1.0-e)*(+1.0))*(-1.0) )/8.0
        self.hess_fc[1,2,1,:] = -(1.0+c)*( ((-1.0)*(+1.0)+0.0)*(1.0-s) + ((-1.0)*(2.0-c+e+s)+(1.0-e)*(+1.0))*(-1.0) )/8.0
        self.hess_fc[1,2,2,:] = -(1.0+c)*( ((+1.0)*(+1.0)+0.0)*(1.0-s) + ((+1.0)*(2.0-c-e+s)+(1.0+e)*(-1.0))*(-1.0) )/8.0
        self.hess_fc[1,2,3,:] = -(1.0-c)*( ((+1.0)*(+1.0)+0.0)*(1.0-s) + ((+1.0)*(2.0+c-e+s)+(1.0+e)*(-1.0))*(-1.0) )/8.0
        self.hess_fc[1,2,4,:] = -(1.0-c)*( ((-1.0)*(-1.0)+0.0)*(1.0+s) + ((-1.0)*(2.0+c+e-s)+(1.0-e)*(+1.0))*(+1.0) )/8.0
        self.hess_fc[1,2,5,:] = -(1.0+c)*( ((-1.0)*(-1.0)+0.0)*(1.0+s) + ((-1.0)*(2.0-c+e-s)+(1.0-e)*(+1.0))*(+1.0) )/8.0
        self.hess_fc[1,2,6,:] = -(1.0+c)*( ((+1.0)*(-1.0)+0.0)*(1.0+s) + ((+1.0)*(2.0-c-e-s)+(1.0+e)*(-1.0))*(+1.0) )/8.0
        self.hess_fc[1,2,7,:] = -(1.0-c)*( ((+1.0)*(-1.0)+0.0)*(1.0+s) + ((+1.0)*(2.0+c-e-s)+(1.0+e)*(-1.0))*(+1.0) )/8.0

        # Mid-edge ξ = 0
        self.hess_fc[1,2,8,:]  = +(1.0+c)*(1.0-c)*(-1.0)*(-1.0)/4.0
        self.hess_fc[1,2,10,:] = +(1.0+c)*(1.0-c)*(+1.0)*(-1.0)/4.0
        self.hess_fc[1,2,12,:] = +(1.0+c)*(1.0-c)*(-1.0)*(+1.0)/4.0
        self.hess_fc[1,2,14,:] = +(1.0+c)*(1.0-c)*(+1.0)*(+1.0)/4.0
        # Mid-edge η = 0
        self.hess_fc[1,2,9,:] = +(1.0+c)*(-2.0*e)*(-1.0)/4.0
        self.hess_fc[1,2,11,:] = +(1.0-c)*(-2.0*e)*(-1.0)/4.0
        self.hess_fc[1,2,13,:] = +(1.0+c)*(-2.0*e)*(+1.0)/4.0
        self.hess_fc[1,2,15,:] = +(1.0-c)*(-2.0*e)*(+1.0)/4.0
        # Mid-edge ζ = 0
        self.hess_fc[1,2,16,:] = +(1.0-c)*(-1.0)*(-2.0*s)/4.0
        self.hess_fc[1,2,17,:] = +(1.0+c)*(-1.0)*(-2.0*s)/4.0
        self.hess_fc[1,2,18,:] = +(1.0+c)*(+1.0)*(-2.0*s)/4.0
        self.hess_fc[1,2,19,:] = +(1.0-c)*(+1.0)*(-2.0*s)/4.0
        self.hess_fc[2,1,:,:] = self.hess_fc[1,2,:,:]

        # dfds
        # Corners
        self.grad_fc[2,0,:]  = -(1.0-c)*(1.0-e)*( (-1.0)*(2.0+c+e+s)+(1.0-s)*(+1.0) )/8.0
        self.grad_fc[2,1,:]  = -(1.0+c)*(1.0-e)*( (-1.0)*(2.0-c+e+s)+(1.0-s)*(+1.0) )/8.0
        self.grad_fc[2,2,:]  = -(1.0+c)*(1.0+e)*( (-1.0)*(2.0-c-e+s)+(1.0-s)*(+1.0) )/8.0
        self.grad_fc[2,3,:]  = -(1.0-c)*(1.0+e)*( (-1.0)*(2.0+c-e+s)+(1.0-s)*(+1.0) )/8.0
        self.grad_fc[2,4,:]  = -(1.0-c)*(1.0-e)*( (+1.0)*(2.0+c+e-s)+(1.0+s)*(-1.0) )/8.0
        self.grad_fc[2,5,:]  = -(1.0+c)*(1.0-e)*( (+1.0)*(2.0-c+e-s)+(1.0+s)*(-1.0) )/8.0
        self.grad_fc[2,6,:]  = -(1.0+c)*(1.0+e)*( (+1.0)*(2.0-c-e-s)+(1.0+s)*(-1.0) )/8.0
        self.grad_fc[2,7,:]  = -(1.0-c)*(1.0+e)*( (+1.0)*(2.0+c-e-s)+(1.0+s)*(-1.0) )/8.0
        # Mid-edge ξ = 0
        self.grad_fc[2,8,:]  = +(1.0+c)*(1.0-c)*(1.0-e)*(-1.0)/4.0
        self.grad_fc[2,10,:] = +(1.0+c)*(1.0-c)*(1.0+e)*(-1.0)/4.0
        self.grad_fc[2,12,:] = +(1.0+c)*(1.0-c)*(1.0-e)*(+1.0)/4.0
        self.grad_fc[2,14,:] = +(1.0+c)*(1.0-c)*(1.0+e)*(+1.0)/4.0
        # Mid-edge η = 0
        self.grad_fc[2,9,:] = +(1.0+c)*(1.0+e)*(1.0-e)*(-1.0)/4.0
        self.grad_fc[2,11,:] = +(1.0-c)*(1.0+e)*(1.0-e)*(-1.0)/4.0
        self.grad_fc[2,13,:] = +(1.0+c)*(1.0+e)*(1.0-e)*(+1.0)/4.0
        self.grad_fc[2,15,:] = +(1.0-c)*(1.0+e)*(1.0-e)*(+1.0)/4.0
        # Mid-edge ζ = 0
        self.grad_fc[2,16,:] = +(1.0-c)*(1.0-e)*(-2.0*s)/4.0
        self.grad_fc[2,17,:] = +(1.0+c)*(1.0-e)*(-2.0*s)/4.0
        self.grad_fc[2,18,:] = +(1.0+c)*(1.0+e)*(-2.0*s)/4.0
        self.grad_fc[2,19,:] = +(1.0-c)*(1.0+e)*(-2.0*s)/4.0

        # d2fds2
        # Corners
        self.hess_fc[2,2,0,:]  = -(1.0-c)*(1.0-e)*( (-1.0)*(+1.0)+(-1.0)*(+1.0) )/8.0
        self.hess_fc[2,2,1,:]  = -(1.0+c)*(1.0-e)*( (-1.0)*(+1.0)+(-1.0)*(+1.0) )/8.0
        self.hess_fc[2,2,2,:]  = -(1.0+c)*(1.0+e)*( (-1.0)*(+1.0)+(-1.0)*(+1.0) )/8.0
        self.hess_fc[2,2,3,:]  = -(1.0-c)*(1.0+e)*( (-1.0)*(+1.0)+(-1.0)*(+1.0) )/8.0
        self.hess_fc[2,2,4,:]  = -(1.0-c)*(1.0-e)*( (+1.0)*(-1.0)+(+1.0)*(-1.0) )/8.0
        self.hess_fc[2,2,5,:]  = -(1.0+c)*(1.0-e)*( (+1.0)*(-1.0)+(+1.0)*(-1.0) )/8.0
        self.hess_fc[2,2,6,:]  = -(1.0+c)*(1.0+e)*( (+1.0)*(-1.0)+(+1.0)*(-1.0) )/8.0
        self.hess_fc[2,2,7,:]  = -(1.0-c)*(1.0+e)*( (+1.0)*(-1.0)+(+1.0)*(-1.0) )/8.0
        # Mid-edge ξ = 0
        self.hess_fc[2,2,8,:]  = +0.0
        self.hess_fc[2,2,10,:] = +0.0
        self.hess_fc[2,2,12,:] = +0.0
        self.hess_fc[2,2,14,:] = +0.0
        # Mid-edge η = 0
        self.hess_fc[2,2,9,:] = +0.0
        self.hess_fc[2,2,11,:] = +0.0
        self.hess_fc[2,2,13,:] = +0.0
        self.hess_fc[2,2,15,:] = +0.0
        # Mid-edge ζ = 0
        self.hess_fc[2,2,16,:] = +(1.0-c)*(1.0-e)*(-2.0)/4.0
        self.hess_fc[2,2,17,:] = +(1.0+c)*(1.0-e)*(-2.0)/4.0
        self.hess_fc[2,2,18,:] = +(1.0+c)*(1.0+e)*(-2.0)/4.0
        self.hess_fc[2,2,19,:] = +(1.0-c)*(1.0+e)*(-2.0)/4.0


class FEIntegrationTriangle(FEGaussIntegration):

    def __init__(self, ng: int, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.ng = ng

        self.c, self.w = self.set_gauss_coefficients_triangle(ng)

        ng = self.ng
        nbf = self.nbf
        ndim = 2

        self.bfn = np.zeros((nbf, ng))
        self.grad_fc = np.zeros((ndim, nbf, ng))
        self.hess_fc = np.zeros((ndim, ndim, nbf, ng))

        if nbf == 3:
            self.set_linear_basis_functions()
        elif nbf == 6:
            self.set_quadratic_basis_functions()
        else:
            raise NotImplementedError(f"Unsupported nbf = {nbf} in FEIntegrationTriangle!")

    def set_gauss_coefficients_triangle(self, ng: int) -> NDArray[np.float64]:

        c_2d = np.zeros((ng,2))
        w_2d = np.zeros(ng)

        if ng == 1:
            # Precision 1
            c_2d[0,0] = 1.0/3.0

            c_2d[0,1] = 1.0/3.0

            w_2d[0] = 1.0

        elif ng == 3:
            # Precision 2
            c_2d[0,0] = 0.5
            c_2d[1,0] = 0.5
            c_2d[2,0] = 0.0 + 1.0e-8

            c_2d[0,1] = 0.0 + 1.0e-8
            c_2d[1,1] = 0.5
            c_2d[2,1] = 0.5

            w_2d[0] = 1.0/3.0
            w_2d[1] = 1.0/3.0
            w_2d[2] = 1.0/3.0

        elif ng == 4:
            # Precision 3
            c_2d[0,0] = 1.0/3.0
            c_2d[1,0] = 0.2
            c_2d[2,0] = 0.2
            c_2d[3,0] = 0.6

            c_2d[0,1] = 1.0/3.0
            c_2d[1,1] = 0.2
            c_2d[2,1] = 0.6
            c_2d[3,1] = 0.2

            w_2d[0] = -27.0/48.0
            w_2d[1] =  25.0/48.0
            w_2d[2] =  25.0/48.0
            w_2d[3] =  25.0/48.0

        elif ng == 6:
            # Precision 4
            c_2d[0,0] = 0.816847572980459
            c_2d[1,0] = 0.091576213509771
            c_2d[2,0] = c_2d[1,0]
            c_2d[3,0] = 0.108103018168070
            c_2d[4,0] = 0.445948490915965
            c_2d[5,0] = c_2d[4,0]

            c_2d[0,1] = 0.091576213509771
            c_2d[1,1] = 0.816847572980459
            c_2d[2,1] = c_2d[0,1]
            c_2d[3,1] = 0.445948490915965
            c_2d[4,1] = 0.108103018168070
            c_2d[5,1] = c_2d[3,1]

            w_2d[0:3] = 0.109951743655322
            w_2d[3:6] = 0.223381589678011

        elif ng == 7:
            # Precision 5
            c_2d[0,0] = 1.0/3.0
            c_2d[1,0] = 0.10128650732345633
            c_2d[2,0] = 0.79742698535308720
            c_2d[3,0] = c_2d[1,0]
            c_2d[4,0] = 0.47014206410511505
            c_2d[5,0] = 0.05971587178976981
            c_2d[6,0] = c_2d[4,0]

            c_2d[0,1] = 1.0/3.0
            c_2d[1,1] = c_2d[1,0]
            c_2d[2,1] = c_2d[1,0]
            c_2d[3,1] = c_2d[2,0]
            c_2d[4,1] = c_2d[4,0]
            c_2d[5,1] = c_2d[4,0]
            c_2d[6,1] = c_2d[5,0]

            w_2d[0] = 0.225
            w_2d[1:4] = 0.12593918054482717
            w_2d[4:7] = 0.13239415278850616

        else:
            raise NotImplementedError(f'Unsupported ng = {ng} in set_gauss_coefficients_triangle!')

        w_2d /= 2.0

        return c_2d, w_2d

    def set_linear_basis_functions(self):

        c = self.c[:,0]
        e = self.c[:,1]

        self.bfn[0,:] = (1.0-c-e)
        self.bfn[1,:] = c
        self.bfn[2,:] = e

        # dfdc
        self.grad_fc[0,0,:] = -1.0
        self.grad_fc[0,1,:] = +1.0
        self.grad_fc[0,2,:] =  0.0

        # dfde
        self.grad_fc[1,0,:] = -1.0
        self.grad_fc[1,1,:] =  0.0
        self.grad_fc[1,2,:] = +1.0

    def set_quadratic_basis_functions(self):

        c = self.c[:,0]
        e = self.c[:,1]

        self.bfn[0,:] = (1.0-c-e)*(1.0-2.0*c-2.0*e)
        self.bfn[1,:] = c*(2.0*c-1.0)
        self.bfn[2,:] = e*(2.0*e-1.0)
        self.bfn[3,:] = 4.0*c*(1.0-c-e)
        self.bfn[4,:] = 4.0*c*e
        self.bfn[5,:] = 4.0*e*(1.0-c-e)

        # dφ/dξ
        self.grad_fc[0,0,:] = -3.0+4.0*c+4.0*e
        self.grad_fc[0,1,:] = +4.0*c-1.0
        self.grad_fc[0,2,:] = +0.0
        self.grad_fc[0,3,:] = +4.0-8.0*c-4.0*e
        self.grad_fc[0,4,:] = +4.0*e
        self.grad_fc[0,5,:] = -4.0*e

        # d2φ/dξ2
        self.hess_fc[0,0,0,:] = +4.0
        self.hess_fc[0,0,1,:] = +4.0
        self.hess_fc[0,0,2,:] = +0.0
        self.hess_fc[0,0,3,:] = -8.0
        self.hess_fc[0,0,4,:] = +0.0
        self.hess_fc[0,0,5,:] = +0.0

        # d2φ/dξdη
        self.hess_fc[0,1,0,:] = +4.0
        self.hess_fc[0,1,1,:] = +0.0
        self.hess_fc[0,1,2,:] = +0.0
        self.hess_fc[0,1,3,:] = -4.0
        self.hess_fc[0,1,4,:] = +4.0
        self.hess_fc[0,1,5,:] = -4.0
        self.hess_fc[1,0,:,:] = self.hess_fc[0,1,:,:]

        # dφ/dη
        self.grad_fc[1,0,:] = -3.0+4.0*c+4.0*e
        self.grad_fc[1,1,:] = +0.0
        self.grad_fc[1,2,:] = +4.0*e-1.0
        self.grad_fc[1,3,:] = -4.0*c
        self.grad_fc[1,4,:] = +4.0*c
        self.grad_fc[1,5,:] = +4.0-4.0*c-8.0*e

        # d2φ/dη2
        self.hess_fc[1,1,0,:] = +4.0
        self.hess_fc[1,1,1,:] = +0.0
        self.hess_fc[1,1,2,:] = +4.0
        self.hess_fc[1,1,3,:] = +0.0
        self.hess_fc[1,1,4,:] = +0.0
        self.hess_fc[1,1,5,:] = -8.0


class FEIntegrationTetrahedron(FEGaussIntegration):

    def __init__(self, ng: int, *args, **kwargs):
        super().__init__(*args, **kwargs)

        self.ng = ng

        self.c, self.w = self.set_gauss_coefficients_tetrahedron(ng)

        ng = self.ng
        nbf = self.nbf
        ndim = 3

        self.bfn = np.zeros((nbf, ng))
        self.grad_fc = np.zeros((ndim, nbf, ng))
        self.hess_fc = np.zeros((ndim, ndim, nbf, ng))

        if nbf == 4:
            self.set_linear_basis_functions()
        elif nbf == 10:
            self.set_quadratic_basis_functions()
        else:
            raise NotImplementedError(f"Unsupported nbf = {nbf} in FEIntegrationTetrahedron!")

    def set_gauss_coefficients_tetrahedron(self, ng: int) -> NDArray[np.float64]:

        c_3d = np.zeros((ng,3))
        w_3d = np.zeros(ng)

        if ng == 4:

            # Precision 2
            c_3d[:,0] = [0.5854101966249685, 0.1381966011250105,
                0.1381966011250105, 0.1381966011250105]

            c_3d[:,1] = [0.1381966011250105, 0.5854101966249685,
                0.1381966011250105, 0.1381966011250105]

            c_3d[:,2] = [0.1381966011250105, 0.1381966011250105,
                0.5854101966249685, 0.1381966011250105]

            w_3d[:] = 0.25

        elif ng == 5:
            # Precision 3
            c_3d[:,0] = [0.25, 1.0/6.0, 1.0/6.0, 1.0/6.0, 0.50]

            c_3d[:,1] = [0.25, 1.0/6.0, 1.0/6.0, 0.5, 1.0/6.0]

            c_3d[:,2] = [0.25, 1.0/6.0, 0.50, 1.0/6.0, 1.0/6.0]

            w_3d[:] = [-0.80, 0.45, 0.45, 0.45, 0.45]

        elif ng == 10:
            # Precision 3
            c_3d[:,0] = [0.5684305841968444, 0.1438564719343852,
                0.1438564719343852, 0.1438564719343852, 0.0, 0.5, 0.5,
                0.5, 0.0, 0.0]

            c_3d[:,1] = [0.1438564719343852, 0.1438564719343852,
                0.1438564719343852, 0.5684305841968444, 0.5, 0.0, 0.5,
                0.0, 0.5, 0.0]

            c_3d[:,2] = [0.1438564719343852, 0.1438564719343852,
                0.5684305841968444, 0.1438564719343852, 0.5, 0.5, 0.0,
                0.0, 0.0, 0.5]

            w_3d[:] = [0.2177650698804054, 0.2177650698804054,
                0.2177650698804054, 0.2177650698804054, 0.0214899534130631,
                0.0214899534130631, 0.0214899534130631, 0.0214899534130631,
                0.0214899534130631, 0.0214899534130631]

        elif ng == 11:
            # Precision 4
            c_3d[:,0] = [0.25, 0.7857142857142857, 0.0714285714285714,
                0.0714285714285714, 0.0714285714285714, 0.1005964238332008,
                0.3994035761667992, 0.3994035761667992, 0.3994035761667992,
                0.1005964238332008, 0.1005964238332008]

            c_3d[:,1] = [0.25, 0.0714285714285714, 0.0714285714285714,
                0.0714285714285714, 0.7857142857142857, 0.3994035761667992,
                0.1005964238332008, 0.3994035761667992, 0.1005964238332008,
                0.3994035761667992, 0.1005964238332008]

            c_3d[:,2] = [0.25, 0.0714285714285714, 0.0714285714285714,
                0.7857142857142857, 0.0714285714285714, 0.3994035761667992,
                0.3994035761667992, 0.1005964238332008, 0.1005964238332008,
                0.1005964238332008, 0.3994035761667992]

            w_3d[:] = [-0.0789333333333333, 0.0457333333333333,
                0.0457333333333333, 0.0457333333333333, 0.0457333333333333,
                0.1493333333333333, 0.1493333333333333, 0.1493333333333333,
                0.1493333333333333, 0.1493333333333333, 0.1493333333333333]

        elif ng == 15:
            # Precision 5
            c_3d[:,0] = [0.25, 0.00, 1.0/3.0, 1.0/3.0, 1.0/3.0,
                0.7272727272727273, 1.0/11.0, 1.0/11.0, 1.0/11.0,
                0.4334498464263357, 0.0665501535736643, 0.0665501535736643,
                0.0665501535736643, 0.4334498464263357, 0.4334498464263357]

            c_3d[:,1] = [0.25, 1.0/3.0, 1.0/3.0, 1.0/3.0, 0.0,
                1.0/11.0, 1.0/11.0, 1.0/11.0, 0.7272727272727273,
                0.0665501535736643, 0.4334498464263357, 0.0665501535736643,
                0.4334498464263357, 0.0665501535736643, 0.4334498464263357]

            c_3d[:,2] = [0.25, 1.0/3.0, 1.0/3.0, 0.0, 1.0/3.0,
                1.0/11.0, 1.0/11.0, 0.7272727272727273, 1.0/11.0,
                0.0665501535736643, 0.0665501535736643, 0.4334498464263357,
                0.4334498464263357, 0.4334498464263357, 0.0665501535736643]

            w_3d[:] = [0.1817020685825351, 0.0361607142857143,
                0.0361607142857143, 0.0361607142857143, 0.0361607142857143,
                0.0698714945161738, 0.0698714945161738, 0.0698714945161738,
                0.0698714945161738, 0.0656948493683187, 0.0656948493683187,
                0.0656948493683187, 0.0656948493683187, 0.0656948493683187,
                0.0656948493683187]

        else:
            raise NotImplementedError(f'Unsupported ng = {ng} in set_gauss_coefficients_tetrahedron!')

        w_3d /= 6.0

        return c_3d, w_3d

    def set_linear_basis_functions(self):

        c = self.c[:,0]
        e = self.c[:,1]
        s = self.c[:,2]

        self.bfn[0,:]  = 1.0-c-e-s
        self.bfn[1,:]  = c
        self.bfn[2,:]  = e
        self.bfn[3,:]  = s

        # dφ/dξ
        self.grad_fc[0,0,:] = -1.0
        self.grad_fc[0,1,:] = +1.0
        self.grad_fc[0,2,:] =  0.0
        self.grad_fc[0,3,:] =  0.0

        # dφ/dη
        self.grad_fc[1,0,:] = -1.0
        self.grad_fc[1,1,:] =  0.0
        self.grad_fc[1,2,:] = +1.0
        self.grad_fc[1,3,:] =  0.0

        # dφ/dζ
        self.grad_fc[2,0,:] = -1.0
        self.grad_fc[2,1,:] =  0.0
        self.grad_fc[2,2,:] =  0.0
        self.grad_fc[2,3,:] = +1.0

    def set_quadratic_basis_functions(self):

        c = self.c[:,0]
        e = self.c[:,1]
        s = self.c[:,2]

        self.bfn[0,:] = +(1.0-c-e-s)*(1.0-2.0*c-2.0*e-2.0*s)
        self.bfn[1,:] = +c*(2.0*c-1.0)
        self.bfn[2,:] = +e*(2.0*e-1.0)
        self.bfn[3,:] = +s*(2.0*s-1.0)
        self.bfn[4,:] = +4.0*(1.0-c-e-s)*c
        self.bfn[5,:] = +4.0*c*e
        self.bfn[6,:] = +4.0*(1.0-c-e-s)*e
        self.bfn[7,:] = +4.0*(1.0-c-e-s)*s
        self.bfn[8,:] = +4.0*c*s
        self.bfn[9,:] = +4.0*e*s

        # dφ/dξ
        self.grad_fc[0,0,:] = -3.0+4.0*(c+e+s)
        self.grad_fc[0,1,:] = +4.0*c-1.0
        self.grad_fc[0,2,:] = +0.0
        self.grad_fc[0,3,:] = +0.0
        self.grad_fc[0,4,:] = +4.0*(1.0-2.0*c-e-s)
        self.grad_fc[0,5,:] = +4.0*e
        self.grad_fc[0,6,:] = -4.0*e
        self.grad_fc[0,7,:] = -4.0*s
        self.grad_fc[0,8,:] = +4.0*s
        self.grad_fc[0,9,:] = +0.0

        # d2φ/dξ2
        self.hess_fc[0,0,0,:] = +4.0
        self.hess_fc[0,0,1,:] = +4.0
        self.hess_fc[0,0,2,:] = +0.0
        self.hess_fc[0,0,3,:] = +0.0
        self.hess_fc[0,0,4,:] = -8.0
        self.hess_fc[0,0,5,:] = +0.0
        self.hess_fc[0,0,6,:] = +0.0
        self.hess_fc[0,0,7,:] = +0.0
        self.hess_fc[0,0,8,:] = +0.0
        self.hess_fc[0,0,9,:] = +0.0

        # d2φ/dξdη
        self.hess_fc[0,1,0,:] = +4.0
        self.hess_fc[0,1,1,:] =  0.0
        self.hess_fc[0,1,2,:] = +0.0
        self.hess_fc[0,1,3,:] = +0.0
        self.hess_fc[0,1,4,:] = -4.0
        self.hess_fc[0,1,5,:] = +4.0
        self.hess_fc[0,1,6,:] = -4.0
        self.hess_fc[0,1,7,:] =  0.0
        self.hess_fc[0,1,8,:] =  0.0
        self.hess_fc[0,1,9,:] = +0.0
        self.hess_fc[1,0,:,:] = self.hess_fc[0,1,:,:]

        # d2φ/dξdζ
        self.hess_fc[0,2,0,:] = +4.0
        self.hess_fc[0,2,1,:] =  0.0
        self.hess_fc[0,2,2,:] = +0.0
        self.hess_fc[0,2,3,:] = +0.0
        self.hess_fc[0,2,4,:] = -4.0
        self.hess_fc[0,2,5,:] =  0.0
        self.hess_fc[0,2,6,:] =  0.0
        self.hess_fc[0,2,7,:] = -4.0
        self.hess_fc[0,2,8,:] = +4.0
        self.hess_fc[0,2,9,:] = +0.0
        self.hess_fc[2,0,:,:] = self.hess_fc[0,2,:,:]

        # dφ/dη
        self.grad_fc[1,0,:] = -3.0+4.0*(c+e+s)
        self.grad_fc[1,1,:] = +0.0
        self.grad_fc[1,2,:] = +4.0*e-1.0
        self.grad_fc[1,3,:] = +0.0
        self.grad_fc[1,4,:] = -4.0*c
        self.grad_fc[1,5,:] = +4.0*c
        self.grad_fc[1,6,:] = +4.0*(1.0-c-2.0*e-s)
        self.grad_fc[1,7,:] = -4.0*s
        self.grad_fc[1,8,:] = +0.0
        self.grad_fc[1,9,:] = +4.0*s

        # d2φ/dη2
        self.hess_fc[1,1,0,:] = +4.0
        self.hess_fc[1,1,1,:] = +0.0
        self.hess_fc[1,1,2,:] = +4.0
        self.hess_fc[1,1,3,:] = +0.0
        self.hess_fc[1,1,4,:] = +0.0
        self.hess_fc[1,1,5,:] = +0.0
        self.hess_fc[1,1,6,:] = -8.0
        self.hess_fc[1,1,7,:] = +0.0
        self.hess_fc[1,1,8,:] = +0.0
        self.hess_fc[1,1,9,:] = +0.0

        # d2φ/dηdζ
        self.hess_fc[1,2,0,:] = +4.0
        self.hess_fc[1,2,1,:] = +0.0
        self.hess_fc[1,2,2,:] =  0.0
        self.hess_fc[1,2,3,:] = +0.0
        self.hess_fc[1,2,4,:] =  0.0
        self.hess_fc[1,2,5,:] =  0.0
        self.hess_fc[1,2,6,:] = -4.0
        self.hess_fc[1,2,7,:] = -4.0
        self.hess_fc[1,2,8,:] = +0.0
        self.hess_fc[1,2,9,:] = +4.0
        self.hess_fc[2,1,:,:] = self.hess_fc[1,2,:,:]

        # dφ/dζ
        self.grad_fc[2,0,:] = -3.0+4.0*(c+e+s)
        self.grad_fc[2,1,:] = +0.0
        self.grad_fc[2,2,:] = +0.0
        self.grad_fc[2,3,:] = +4.0*s-1.0
        self.grad_fc[2,4,:] = -4.0*c
        self.grad_fc[2,5,:] = +0.0
        self.grad_fc[2,6,:] = -4.0*e
        self.grad_fc[2,7,:] = +4.0*(1.0-c-e-2.0*s)
        self.grad_fc[2,8,:] = +4.0*c
        self.grad_fc[2,9,:] = +4.0*e

        # d2φ/dζ2
        self.hess_fc[2,2,0,:] = +4.0
        self.hess_fc[2,2,1,:] = +0.0
        self.hess_fc[2,2,2,:] = +0.0
        self.hess_fc[2,2,3,:] = +4.0
        self.hess_fc[2,2,4,:] = +0.0
        self.hess_fc[2,2,5,:] = +0.0
        self.hess_fc[2,2,6,:] = +0.0
        self.hess_fc[2,2,7,:] = -8.0
        self.hess_fc[2,2,8,:] = +0.0
        self.hess_fc[2,2,9,:] = +0.0
