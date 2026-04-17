# Numerical Methods: Roots of Equations

We seek solutions to non-linear equations:

$$ f(x) = 0 $$

## Bracketing Methods
We seek for solutions in a closed interval $[a,b]$ where the function changes sign.
Bracketing methods are slow, but will converge to a root provided that $f(a)f(b)<0$.
However, they will miss other possible roots.

### Bisection
We start $(k = 1)$ with the initial interval $[x_l^{(k)},x_u^{(k)}]=[a,b]$.
We calculate $f(x_r^{(k)})$ with $x_r^{(k)} = \frac{x_l^{(k)}+x_u^{(k)}}{2}$:

* If: $f(x_l^{(k)})f(x_r^{(k)}) < 0: x_u^{(k+1)} = x_r^{(k)}$. The root is in the lower sub-interval.
* If: $f(x_l^{(k)})f(x_r^{(k)}) > 0: x_l^{(k+1)} = x_r^{(k)}$. The root is in the upper sub-interval.

Thus, we iteratively decrease the searching interval, $\delta x^{(k+1)} = x_u^{(k+1)}-x_l^{(k-1)}$.\
The absolute error is related to the interval:

$$ e_a^{(k+1)} = \left|x_r^{(k+1)}-x_r^{(k)} \right| =
\frac{\delta x^{(k+1)}}{2} =
\frac{\delta x^{(1)}}{2^k}$$

Thus, to reach a desired absolute error, $e$, we need $k = log_2(\frac{\delta x^{(1)}}{e})$ iterations.
The approximate error is greater than the true error.

### False position (Regula falsi)
We check the same conditions to identify the sub-interval of the root.
The solution is updated by taking into account the magnitudes $f(x_l), f(x_u)$.
The intersection of x-axis with the line connecting $f(x_l)$ and $f(x_u)$ is the new root estimation:

$$ x_r^{(k+1)} = x_u^{(k)} - \frac{f(x_u^{(k)})(x_l^{(k)}-x_u^{(k)})}{f(x_l^{(k)})-f(x_u^{(k)})} $$

The False position is faster than bisection except for when one of the two bounds remains stagnant (unchanged for more than $k$ iterations).
To circumvent that, the function value is divided by half at the stagnant bound.

## Open Methods
Here, we do not search in a closed interval. Instead we seek for a root close to a point (initial condition).
The choice of the initial condition will dictate whether the method will converge or not.
Upon convergence, open methods are usually faster than bracketing methods.

### Fixed-point iteration
We reformulate the function $f(x) = 0 \to g(x) = x $. We calculate the new estimation as:

$$ x^{(k+1)} = g(x^{(k)}) $$

The method converges linearly if $ |g'(x)| < 1 $. It diverges otherwise.

### Secant
The Secant method approximates the new solution with a straight line connecting the two previous solutions:

$$ x^{(k+1)} = x^{(k)} - f(x^{(k)})\frac{x^{(k)}-x^{(k-1)}}{f(x^{(k)})-f(x^{(k-1)})}$$

### IQI (Inverse Quadratic Interpolation)
The IQI method approximates the new solution with a quadratic line connecting the two previous solutions:

$$ x^{(k+1)} = x^{(k-2)}\frac{f^{(k)}f^{(k-1)}}{(f^{(k-2)}-f^{(k-1)})(f^{(k-2)}-f^{(k)})} $$
$$ + x^{(k-1)}\frac{f^{(k)}f^{(k-2)}}{(f^{(k-1)}-f^{(k-2)})(f^{(k-1)}-f^{(k)})} $$
$$ + x^{(k)}\frac{f^{(k-2)}f^{(k-1)}}{(f^{(k)}-f^{(k-2)})(f^{(k)}-f^{(k-1)})} $$

### Newton-Raphson
At each iteration we calculate the derivative

$$ f'(x^{(k)}) = \frac{f(x^{(k)})-f(x^{(k+1)})}{x^{(k)}-x^{(k+1)}} $$

and we set $f(x^{(k+1)}) = 0$. Thus, the new approximation is:

$$ x^{(k+1)} = x^{(k)} - \frac{f(x^{(k)})}{f'(x^{(k)})}$$

The method converges quadratically, but it performs poorly around inflection points, i.e. $ f''(x) = 0 $ and when $ f'(x) \to 0 $.
When the derivative is calculated using finite differences, we refer to the modified Secant method.

### Ralston-Rabinowitz (Modified Newton-Raphson)
This method is preferable when there are multiple roots around the initial guess.
We seek the solution to $u(x)$:

$$ u(x) = \frac{f(x)}{f'(x)} \to
u'(x) = \frac{[f'(x)]^2 - f(x)f''(x)}{[f'(x)]^2} $$

Thus:

$$ x^{(k+1)} = x^{(k)} - \frac{u(x^{(k)})}{u'(x^{(k)})}\to x^{(k+1)} = x^{(k)} - \frac{f(x)f'(x)}{[f'(x)]^2-f(x)f''(x)} $$

## Hybrid Methods
Hybrid methods combine a reliable but slow bracketing method away from the root, and a fast open method close to the root.

### Brent's Method
It performs an IQI (if possible) or Secant step unless some conditions are met which force bisection.\
The new root estimation, $s$, replaces the best, $b$, if $f(a)f(s) < 0$. Otherwise, $s$ replaces the other end of the bracket, $a$.\
It is similar to MATLAB's fzero function.

### Ridders' Method
It evaluates the function at the midpoint:

$$ x_m^{(k+1)} = \frac{x_1^{(k)}+x_2^{(k)}}{2} $$

The false position method is applied at $x_0, x_1$ for the function:

$$ h(x) = f(x)e^{(ax)} $$

For which:

$$ h(x_m) = \frac{h(x_1)+h(x_2)}{2} $$

This results in the new root estimation:

$$ x_{r}^{(k+1)} = x_{m}^{(k+1)}+(x_{m}^{(k+1)}-x_{1}^{(k)})
\frac {sign[f_1^{(k)}]f_m^{(k+1)}}
{\sqrt {{f_m^{(k+1)}}^{2}-f_1^{(k)}f_2^{(k)}}} $$

If $f_m^{(k+1)}f_r^{(k+1)} < 0$, $x_1, x_2$ are replaced by $x_m, x_r$.\
If $f_1^{(k+1)}f_r^{(k+1)} < 0$, $x_r$ replaces $x_2$.\
Otherwise, $f_2^{(k+1)}f_r^{(k+1)} < 0$ and $x_r$ replaces $x_1$.

### Chandrupatla's Method
Performs IQI if:

$$ 0 < \phi < 1, \; \xi < 1-\frac{1}{1-\phi} $$
$$ \phi = \frac{f_a - f_b}{f_c - f_b} $$
$$ \xi = \frac{a-b}{c-b} $$

Otherwise, it performs bisection.
The new estimate, $s$, replaces $b$ if $f(a)f(s) < 0$.
Otherwise, it replaces $a$.
