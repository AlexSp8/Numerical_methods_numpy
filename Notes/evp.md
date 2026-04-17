# Numerical Methods: Eigenvalue Problems (EVP)

An eigenvalue problem is a homogeneous linear algebraic system of the form:

$$ (\underline{\underline{A}} - \lambda \underline{\underline{I}}) \underline{x} = \underline{0} $$

$\lambda$ are the eigenvalues and $x$ is a solution (eigenvector).
The problem lies in determining the $\lambda _i$.
For such homogeneous systems this occurs when the determinant is zero:

$$ det(\underline{\underline{A}} - \lambda \underline{\underline{I}}) = 0 $$

## Polynomial method
Expansion of the determinant leads to a polynomial whose roots can be determined with the polynomial method.

## Power method
Starting with an initial guess for an eigenvector and the form:

$$ \underline{\underline{A}} \; \underline{x} = \lambda \underline{x} $$

we update $\lambda$ and multiply the matrix with the normalized resulting eigenvector. The process is repeated until convergence at the largest eigenvalue.
The next highest eigenvalues can be obtained by repeating the process after deflation (replacing the matrix with only the remaining eigenvalues).
