# Numerical Methods: Boundary-Value Problems (BVP)

Initial-Value Problems are ODEs in which there is only a single (initial) condition for the independent variable.
These are $1^{st}$ order (systems of) ODEs. \
Boundary-Value Problems are ODEs in which there are (boundary) conditions at different values of the independent variable.
These are higher order (systems of) ODEs. Common BVPs are 1D Mass Diffusion, 1D Heat Conduction etc..
These are second order ODEs of the general form:

$$ \frac{d^2y}{dx^2} + f(x,y,\frac{dy}{dx}) = 0 $$

At the boundaries we can have combinations of Dirichlet, Neumann or Robin conditions.\
**Dirichlet**:
$$ y(x_0)=v_0, \; y(x_n)=v_f $$
**Neumann**:
$$ \frac{dy}{dx}(x_0)=v_0, \; \frac{dy}{dx}(x_n)=v_f $$
**Robin**:
$$ \frac{dy}{dx}(x_0)+a_0y_0=v_0, \; \frac{dy}{dx}(x_n)+a_ny(x_n)=v_n $$

The BVP can also be a system of equations.

## Finite Differences
A common approach is to approximate the derivative with finite differences. For example, using central FD $O(h^2)$ for the second derivative and forward FD $O(h)$ for the first derivative (to avoid oscillations) at each internal point:

$$ \frac{y_{i+1}-2y_i+y_{i-1}}{h^2} + f(x_i,y_i,\frac{y_{i+1}-y_i}{h}) = 0 $$

Or for non-constant $h$:

$$ \frac{2}{h_{i+1}+h_i}(\frac{y_{i+1}-y_i}{h_{i+1}}-\frac{y_i-y_{i-1}}{h_i})
+ f(x_i,y_i,\frac{y_{i+1}-y_i}{h_{i+1}}) = 0 $$

where $h_{i+1} = x_{i+1}-x_i, h_i = x_i-x_{i-1}$.
At the boundaries: \
**Dirichlet**:
$$ y_0=v_0, \; y_n=v_f $$
**Neumann**:
$$ \frac{y_1-y_0}{h_0}=v_0, \; \frac{y_n-y_{n-1}}{h_n}=v_f $$
**Robin**:
$$ \frac{y_1-y_0}{h_0}+a_0y_0=v_0, \; \frac{y_n-y_{n-1}}{h_n}+a_ny_n=v_f $$
