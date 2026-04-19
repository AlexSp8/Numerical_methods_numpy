# Numerical Methods: Optimization

We can classify optimization problems based on whether they are unconstrained or constrained.
Further classification involves the number of equations (single or system), and the dimensionality (single or multi-variable equations).
We can use bracketing methods to search in a particular interval or we can use open methods around an initial guess.
With direct methods we only evaluate the function without its derivatives.
With gradient methods we also evaluate the derivatives of the function. Thus, we refer to first-order (gradient) or second-order (hessian) methods.

## Unconstrained optimization

### 1D unconstrained optimization
We have a single function that depends on a single variable, x.
We want to find the value of x where the function exhibits a local extreme.

#### Bracketing
These methods are similar to the bracketing methods for finding roots of equations.
We search for an extreme of the function at a specified interval.
These methods will find the optimum value in the interval if the function is unimodal.

##### Golden-section search
It is the analog of the bisection method for finding the roots of an equation.
We start with two points, $x_l$ and $x_u$, that bracket a local extremum of the function.
We choose two more points:

$$x_1 = x_l + d, \; x_2 = x_u - d$$

where $d = \frac{\sqrt{5}-1}{2}(x_u - x_l)$ is the distance between the interval points multiplied by the golden ratio.
Thus, $x_l < x_2 < x_1 < x_u$. Next, we evaluate $f(x_1), f(x_2)$.
When searching for a maximum:
*   *If $f(x_1) > f(x_2)$:*\
The maximum is in the $[x_2, x_1, x_u]$ interval.\
We disregard the domain $[x_l, x_2]$.\
$x_l$ is replaced by $x_2$ ($x_l^{(k)} \to x_2^{(k-1)}$).\
$x_2$ is replaced by $x_1$ ($x_2^{(k)} \to x_1^{(k-1)}$).\
The new $x_1^{(k)} = x_l^{(k)} + d^{(k)}$ is the new optimum.\
The domain reduces by a factor of $\frac{\sqrt{5}-1}{2}$:

$$x_{int}^{(k)} = x_u^{(k)} - x_l^{(k)} = x_u^{(k)} - x_2^{(k-1)} = $$
$$x_u^{(k)} - x_u^{(k-1)} + d^{(k-1)} \to $$
$$x_{int}^{(k)} = d^{(k-1)}$$

*   *If $f(x_2) > f(x_1)$:*\
The maximum is in the $[x_l, x_2, x_1]$ interval.\
We remove the domain from $x_1$ to $x_u$.\
$x_u$ is replaced by $x_1$ ($x_u^{(k)} \to x_1^{(k-1)}$)\
$x_1$ is replaced by $x_2$ ($x_1^{(k)} \to x_2^{(k-1)}$)\
The new $x_2^{(k)} = x_u^{(k)} - d^{(k)}$ is the new optimum.\
The domain again reduces by a factor of $\frac{\sqrt{5}-1}{2}$:

$$x_{int}^{(k)} = x_u^{(k)} - x_l^{(k)} = x_1^{(k-1)} - x_l^{(k)} = $$
$$x_l^{(k-1)} + d^{(k-1)} - x_l^{(k)} \to $$
$$x_{int}^{(k)} = d^{(k-1)}$$

#### Open methods

##### Parabolic interpolation
A $2^{\text{nd}}$ order polynomial is often a good approximation of a function near an optimum value.
If we fit the function with a parabola using 3 points ($x_0 < x_1 < x_2$) around an interval, we can differentiate the parabola and set it to 0.
This will give an estimation of the optimal $x_{opt}$:

$$x_{opt} = x_2+\frac{1}{2} \frac{(x_2^2 - x_0^2)[f(x_2)-f(x_1)] - (x_2^2 - x_1^2)[f(x_2)-f(x_0)]}
{(x_2 - x_0)[f(x_2)-f(x_1)] - (x_2 - x_1)[f(x_2)-f(x_0)]}$$

Depending on the values of $x_{opt}$ and $f(x_{opt})$ we update the interval and repeat the process similar to the golden-section search method.

##### Secant method
We can use the Secant method for finding roots of equations to find the root of the derivative $f'(x)$.

$$x^{(k)} = x^{(k-1)} - f'(x^{(k)})\frac{x^{(k)}-x^{(k-1)}}{f'(x^{(k)})-f'(x^{(k-1)})}$$

##### Newton’s method
We can use Newton’s method for finding roots of equations to find the root of the derivative $f'(x)$.
If we define $g(x) = f'(x)$, we can iteratively find $g(x) = 0$:

$$x^{(k)} = x^{(k-1)} - \frac{g(x^{(k)})}{g'(x^{(k)})} \to x^{(k)} = x^{(k-1)} - \frac{f'(x^{(k)})}{f''(x^{(k)})}$$

This can be derived from $2^{\text{nd}}$ order Taylor series expansion of $f(x)$ and setting the derivative $f'(x) = 0$.

#### Hybrid Methods

##### Brent’s method
Analogous to Brent’s method for root finding, it uses parabolic interpolation whenever possible to converge fast, but reverts to golden-section search when necessary.

### Multi-dimensional unconstrained optimization
#### Direct methods
In direct methods, we find the optimum point $\underline{x}$ without evaluating derivatives of the function $f(\underline{x})$.

##### Powell’s method
We start with an initial guess $\underline{x}_0$ and a set of direction vectors, $\underline{d}_i$ along each dimension $i = 1, \dots, n$, where $n$ is the dimensionality of the variables (how many variables the function $f$ depends on).
At each direction, $i$, we find the optimal step size, $\lambda$, such that $f(\underline{x} + \lambda \underline{d}_i)$ is minimized or maximized.
To this end, we solve a 1D optimization problem. We update each coordinate $j$ of the optimum point:

$$x_j^{(k)} = x_j^{(k-1)} + \lambda_{i,opt}d_{ij}$$

We calculate a new direction vector $\underline{d}_n^{(k)} = \underline{x}^{(k)} - \underline{x}^{(k-1)}$ which is appended in direction vectors.
We also remove from the direction vectors, the vector, $\underline{d}_i$, that caused the biggest change in $f$, i.e. $|f(\underline{x}^{(k)}) - f(\underline{x}^{(k-1)})| = |f(\underline{x}^{(k-1)} + \lambda_{i,opt}\underline{d}_i) - f(\underline{x}^{(k-1)})| \to max$. The process is repeated until $|f(\underline{x}^{(k)}) - f(\underline{x}^{(k-1)})| < \varepsilon$.

#### Gradient methods
For gradient methods, we need to evaluate the function’s derivatives. We need the gradient:

$$\underline{\nabla} f = \left[ \frac{\partial f}{\partial x_1} \quad \frac{\partial f}{\partial x_2} \quad \dots \quad \frac{\partial f}{\partial x_n} \right]^T$$

We also need the Hessian:
$$\underline{\underline{H}} f =
\begin{bmatrix}
\frac{\partial^2 f}{\partial x_1^2} & \frac{\partial^2 f}{\partial x_1 \partial x_2} & \dots & \frac{\partial^2 f}{\partial x_1 \partial x_n} \\
\frac{\partial^2 f}{\partial x_2 \partial x_1} & \frac{\partial^2 f}{\partial x_2^2} & \dots & \frac{\partial^2 f}{\partial x_2 \partial x_n} \\
\vdots & \vdots & \ddots & \vdots \\
\frac{\partial^2 f}{\partial x_n \partial x_1} & \frac{\partial^2 f}{\partial x_n \partial x_2} & \dots & \frac{\partial^2 f}{\partial x_n^2}
\end{bmatrix}$$

##### Steepest ascent/descent
The method resembles Powell’s method, but instead of a set of direction vectors which we update in every iteration, we use $\underline{\nabla}f$ as the direction vector in every iteration.
We start at an initial point $\underline{x}_0$, where we evaluate the gradient $\underline{\nabla}f$.
We move along this direction until $f$ reaches a local extreme at point $\underline{x}^{(k)}$ where we re-evaluate the gradient.
To find the local extreme we need to transform $f(\underline{x})$ to $g(h)$, where $h$ is the gradient direction.

$$g(h) = f(\underline{x} + h\underline{\nabla}f)$$

Then, we find the extreme, i.e. $g'(h_{opt}^{(k)}) = 0$, using a 1D optimization algorithm.
Finally, we update:

$$\underline{x}^{(k)} = \underline{x}^{(k-1)} + h_{opt}^{(k)}\underline{\nabla}f$$

The method is linearly convergent. It slows down significantly close to the optimum value where it performs many small “criss-cross” steps.
Each new search direction is orthogonal to the previous one.

##### Conjugate Gradient
The search direction of steepest descent can be improved by making use of conjugate gradients, leading to quadratic convergence.
The algorithm by Fletcher and Reeves imposes that successive gradient search directions are mutually conjugate.
Thus, the new search direction, $\underline{p}$, is orthogonal to all previous directions ($\underline{p}_i^T \cdot \underline{\underline{A}} \cdot \underline{p}_j = 0, \forall \ i \neq j$), instead of always being the gradient without any use of the search history.
The new solution is now:

$$\underline{x}^{(k)} = \underline{x}^{(k-1)} + h_{opt}^{(k)}\underline{p}^{(k-1)}$$

where we find $h_{opt}^{(k)}$ from a 1D optimization algorithm similar to steepest descent.
The new gradient is: $\underline{\nabla}f^{(k)} = \underline{\nabla}f^{(k)}(\underline{x}^{(k)})$. The direction vector is updated as:

$$\underline{p}^{(k)} = \underline{\nabla}f^{(k)} + \beta \underline{p}^{(k-1)}$$

For Fletcher-Reeves $\beta$ is:

$$\beta = \frac{\underline{\nabla}f^{(k)^T} \cdot \underline{\nabla}f^{(k)}}{\underline{\nabla}f^{(k-1)^T} \cdot \underline{\nabla}f^{(k-1)}}$$

For Polak-Ribiere:

$$\beta = \frac{\underline{\nabla}f^{(k)^T} \cdot (\underline{\nabla}f^{(k)} - \underline{\nabla}f^{(k-1)})}{\underline{\nabla}f^{(k-1)^T} \cdot \underline{\nabla}f^{(k-1)}}$$

##### Newton
With Newton’s method we calculate the Hessian of the function as well.
At each iteration we update: $\underline{x}^{(k)} = \underline{x}^{(k-1)} + d\underline{x}^{(k)}$, where $d\underline{x}^{(k)}$ is calculated by:

$$\underline{\underline{Hf}}^{(k-1)} \cdot d\underline{x}^{(k)} = -\underline{\nabla f}^{(k-1)}$$

This is derived from a second-order Taylor series around $\underline{x}^{(k-1)}$:

$$f(\underline{x}) = f(\underline{x}^{(k-1)}) + \underline{\nabla f}^{(k-1)^T} \cdot d\underline{x}^{(k)} + \frac{1}{2} d\underline{x}^{(k)^T} \cdot \underline{\underline{Hf}}^{(k-1)} \cdot d\underline{x}^{(k)} \to $$
$$\underline{\nabla f}(\underline{x}) = \underline{\nabla f}(\underline{x}^{(k-1)}) + \underline{\underline{Hf}}^{(k-1)} \cdot d\underline{x}^{(k)}$$

At the optimum point:
$$\underline{\nabla f}(\underline{x}) = \underline{0} \to
\underline{\underline{Hf}}^{(k-1)} \cdot d\underline{x}^{(k)} = -\underline{\nabla f}(\underline{x}^{(k-1)})$$

##### Marquardt
It modifies the Newton method by using a modified Hessian:

$$\underline{\underline{\tilde{H}f}}^{(k-1)} = \underline{\underline{Hf}}^{(k-1)} + \lambda^{(k-1)} \underline{\underline{I}}$$

When $\lambda \to 0$ we get the Newton method, while when $\lambda \to \infty$ we get the steepest descent method.

##### BFGS
A quasi-Newton method that approximates the Hessian instead of calculating it.
We define a search direction vector $\underline{p}$, which is updated at every iteration:

$$\underline{p}^{(k-1)} = -\underline{\underline{\tilde{H}}}^{-1} \cdot \underline{\nabla f}(\underline{x}^{(k-1)})$$

With a 1D optimization algorithm we find the optimal step size, $a^{(k)}$ and update the position:

$$\underline{x}^{(k)} = \underline{x}^{(k-1)} + a^{(k)} \underline{p}^{(k-1)}$$

We also define the differences: $d\underline{x}^{(k)} = \underline{x}^{(k)} - \underline{x}^{(k-1)}$ and $d\underline{y}^{(k)} = \underline{\nabla f}(\underline{x}^{(k)}) - \underline{\nabla f}(\underline{x}^{(k-1)})$.
Then, we update the inverse Hessian approximation (starting from $\underline{\underline{\tilde{H}}}^{-1^{(0)}} = \underline{\underline{I}}$, which corresponds to steepest descent):

$$\underline{\underline{\tilde{H}}}^{-1^{(k)}} = \left( \underline{\underline{I}} - \frac{d\underline{x}^{(k-1)} d\underline{y}^{(k-1)^T}}{d\underline{y}^{(k-1)^T} d\underline{x}^{(k-1)}} \right) \underline{\underline{\tilde{H}}}^{-1^{(k-1)}} \left( \underline{\underline{I}} - \frac{d\underline{y}^{(k-1)} d\underline{x}^{(k-1)^T}}{d\underline{y}^{(k-1)^T} d\underline{x}^{(k-1)}} \right) + \frac{d\underline{x}^{(k-1)} d\underline{x}^{(k-1)^T}}{d\underline{y}^{(k-1)^T} d\underline{x}^{(k-1)}}$$

## Constrained Optimization
We want to find the optimum value of an objective function $f(\underline{x})$ subject to constraints.
The constraints can be equality constraints, $g_i(\underline{x}) = c_i$ for $i = 1, \dots, n$ equality constraints.
They can also be inequality constraints, $h_j(\underline{x}) \ge d_j$ for $j = 1, \dots, m$ inequality constraints.

### Linear
We want to optimize the objective function:

$$Z = c_1x_1 + c_2x_2 + \dots + c_nx_n$$

$x_j$ is the magnitude of an activity $j$, and $c_j$ is the payoff of the activity.
The constraints are generally:

$$a_{i1}x_1 + a_{i2}x_2 + \dots + a_{in}x_n \le b_i$$

$a_{ij}$ is the amount of resource $i$ that is consumed by activity $j$, and $b_i$ is the total amount of resource $i$.
Another common constraint is: $x_i \ge 0$.
The objective function and the constraint specify a linear programming problem.

### Lagrange Multipliers

For equality constraints we can use Lagrange multipliers, $\lambda_i$.
The dimensionality of the problem is $n_x + n_e$, where $n_x$ is the number of original variables, $\underline{x}$, and $n_e$ is the number of equality constraints.
We define the Lagrangian function:

$$L(\underline{x}, \underline{\lambda}) = f(\underline{x}) + \underline{\lambda} \cdot \underline{g}(\underline{x})$$

Then, we set partial derivatives of $L$ to zero (at the optimum the gradient of the objective function, $\underline{\nabla}f$ is parallel to the gradient of the constraints, $\underline{\nabla}g_i$):

$$\frac{\partial L}{\partial \underline{x}} = \frac{\partial f}{\partial \underline{x}} + \underline{\lambda} \cdot \frac{\partial g}{\partial \underline{x}} = \underline{0}$$

$$\frac{\partial L}{\partial \underline{\lambda}} = \underline{g}(\underline{x}) = \underline{0}$$

This system of non-linear equations is solved with the Newton Raphson method iteratively for $\underline{x}$ and $\underline{\lambda}$:

$$\underline{\underline{J}} \cdot d\underline{X} = -\underline{r}$$

$$\begin{bmatrix}
\frac{\partial}{\partial \underline{x}} \left( \frac{\partial L}{\partial \underline{x}} \right) & \frac{\partial}{\partial \underline{\lambda}} \left( \frac{\partial L}{\partial \underline{x}} \right) \\ \frac{\partial}{\partial \underline{x}} \left( \frac{\partial L}{\partial \underline{\lambda}} \right) & 0
\end{bmatrix} \cdot
\begin{bmatrix}
d\underline{x} \\ d\underline{\lambda}
\end{bmatrix} = -
\begin{bmatrix}
\frac{\partial L}{\partial \underline{x}} \\ \frac{\partial L}{\partial \underline{\lambda}}
\end{bmatrix} \to $$

$$\begin{bmatrix}
\frac{\partial}{\partial \underline{x}} \left( \frac{\partial L}{\partial \underline{x}} \right) & \frac{\partial}{\partial \underline{\lambda}} \left( \frac{\partial L}{\partial \underline{x}} \right) \\ \frac{\partial g}{\partial \underline{x}} & 0
\end{bmatrix} \cdot
\begin{bmatrix} d\underline{x} \\ d\underline{\lambda}
\end{bmatrix} = -
\begin{bmatrix} \frac{\partial f}{\partial \underline{x}} + \underline{\lambda} \cdot \frac{\partial g}{\partial \underline{x}} \\ \underline{g}(\underline{x})
\end{bmatrix}$$

The solution is always a saddle point of the Lagrangian function.
Aside from being a mathematical bridge, $\lambda$ has a physical meaning.
In economics and engineering, $\lambda$ is often called the shadow price.
It represents the rate of change of the optimal value of the objective function with respect to the constraint.

### Augmented Lagrangian
Combines Lagrange multipliers, $\lambda$, with penalty.
Instead of the function, $f$, with constraints, $g$, we optimize:

$$ L_A(x,\lambda,\rho) = f(x) + \lambda g(x) +
\rho g^2(x)$$

We perform unconstrained optimization on $L_A$, update $\lambda$ and increase the penalty $\rho$:

$$ \lambda ^{(k+1)} = \lambda ^{(k)} + \rho ^{(k)} g(x) $$
$$ \rho ^{(k+1)} = \rho ^{(k)} \gamma $$

with $\gamma$ a constant.
We repeat until convergence, i.e. $norm(g) < tol$.

