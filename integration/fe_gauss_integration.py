
import numpy as np

class FEGaussIntegration():

    def __init__(self, ng: int, nbf: int):

        self.ng = ng

        self.c, self.w = self.gauss_legendre_coefficients(ng)

        self.nbf = nbf

        if nbf == 2:
            self.bfn, self.dfdc, self.d2fdc2 = self.set_linear_basis_functions()
        elif nbf == 3:
            self.bfn, self.dfdc, self.d2fdc2 = self.set_quadratic_basis_functions()
        else:
            raise ValueError("Unsupported element order in FEGaussIntegration!")

    def gauss_legendre_coefficients(self, ng: int = 1
        ) -> tuple[np.ndarray[tuple[int]], np.ndarray[tuple[int]]]:
        """Solves for Gauss-Legendre points and weights of degree n."""

        points = np.zeros(ng)
        weights = np.zeros(ng)

        # We only need to find half the roots due to symmetry
        m = (ng + 1)//2

        for i in range(m):
            # Initial guess for the i-th root (Chebyshev-like spacing)
            x = np.cos(np.pi*(i + 0.75)/(ng + 0.5))

            for _ in range(100):
                p, dp = self.legendre_polynomial(x, ng)
                dx = p/dp
                x -= dx
                if abs(dx) < 1e-12:
                    break

            points[i] = -x
            points[ng-1-i] = x

            weight = 2.0/((1.0 - x**2)*(dp**2))
            weights[i] = weight
            weights[ng-1-i] = weight

        return points, weights
    
    @staticmethod
    def legendre_polynomial(x: float, ng: int) -> tuple[float, float]:
        """Evaluates the n-th Legendre polynomial and its derivative at x
        using the recursive relation."""

        if ng == 0:
            return 1.0, 0.0

        p_prev = 0.0
        p = 1.0

        for i in range(ng):
            p_next = ((2.0*i + 1.0)*x*p - i*p_prev)/(i + 1.0)
            p_prev = p
            p = p_next

        dp = ng*(x*p - p_prev)/(x**2 - 1)

        return p, dp
    
    def set_linear_basis_functions(self):
    
        c = self.c.copy()
        nbf = self.nbf
        ng = self.ng

        bfn = np.zeros((nbf, ng))
        dfdc = np.zeros((nbf, ng))
        d2fdc2 = np.zeros((nbf, ng))

        bfn[0,:] = (1.0-c)/2.0
        bfn[1,:] = (1.0+c)/2.0
        
        dfdc[0,:] = -1.0/2.0
        dfdc[1,:] = +1.0/2.0

        return bfn, dfdc, d2fdc2
    
    def set_quadratic_basis_functions(self):
        
        c = self.c.copy()
        nbf = self.nbf
        ng = self.ng

        bfn = np.zeros((nbf, ng))
        dfdc = np.zeros((nbf, ng))
        d2fdc2 = np.zeros((nbf, ng))

        bfn[0,:] = -(1.0-c)*c/2.0
        bfn[1,:] = +(1.0-c)*(1.0+c) # mid node
        bfn[2,:] = +(1.0+c)*c/2.0
        
        dfdc[0,:] = -(1.0-2.0*c)/2.0
        dfdc[1,:] = -2.0*c # mid node
        dfdc[2,:] = +(1.0+2.0*c)/2.0
        
        d2fdc2[0,:] = +1.0
        d2fdc2[1,:] = -2.0 # mid node
        d2fdc2[2,:] = +1.0
        
        return bfn, dfdc, d2fdc2
