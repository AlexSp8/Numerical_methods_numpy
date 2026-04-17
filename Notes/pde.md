# Numerical Methods: Partial Differential Equations (PDEs)

We want to solve a general $n$-dimensional PDE of order $m$:

$$ F(x_1, x_2, ..., x_n, u, \frac{\partial u}{\partial x_1}, \frac{\partial ^2u}{\partial x_1^2}, ...,
\frac{\partial ^mu}{\partial x_1^m}, ..., \frac{\partial u}{\partial x_n}, \frac{\partial ^2u}{\partial x_n^2}, ...,
\frac{\partial ^mu}{\partial x_n^m})= 0 $$

which requires $m$ boundary conditions in each direction.\
The general second order PDE:

$$\sum_{i=1}^{n} \sum_{j=1}^{n} a_{ij}(\underline{x}) \frac{\partial^2 u}{\partial x_i \partial x_j} + \sum_{i=1}^{n} b_i(\underline{x}) \frac{\partial u}{\partial x_i} + c(\underline{x})u = f(\underline{x})$$

The general form for physical problems is:

$$\frac{\partial u}{\partial t} + \underbrace{\nabla \cdot \mathbf{f}(u)}_{\text{Flux/Convection}} = \underbrace{\nabla \cdot (D \nabla u)}_{\text{Diffusion}} + \underbrace{S(u)}_{\text{Source}}$$

For $D=0$ the equation is hyperbolic.
At steady state $\frac{\partial u}{\partial t}=0$, it is elliptic.
