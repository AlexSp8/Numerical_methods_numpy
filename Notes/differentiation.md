# Numerical Methods: Differentiation

## Finite Differences
We can approximate derivatives around a value x with finite differences (FD). We can use forward points (forward FD), backward points (backward FD), or points in both directions (central FD). For example, using two points:

**Forward $O(h)$:**
$$f'(x_i) = \frac{f(x_{i+1}) - f(x_i)}{h}$$

**Backward $O(h)$:**
$$f'(x_i) = \frac{f(x_i) - f(x_{i-1})}{h}$$

**Central $O(h^2)$:**
$$f'(x_i) = \frac{f(x_{i+1}) - f(x_{i-1})}{2h}$$

We can also approximate higher order derivatives. For example, the second derivative using 3 points:

**Forward $O(h)$:**
$$f''(x_i) = \frac{f(x_{i+2}) - 2f(x_{i+1}) + f(x_i)}{h^2}$$

**Backward $O(h)$:**
$$f''(x_i) = \frac{f(x_i) - 2f(x_{i-1}) + f(x_{i-2})}{h^2}$$

**Central $O(h^2)$:**
$$f''(x_i) = \frac{f(x_{i+1}) - 2f(x_i) + f(x_{i-1})}{h^2}$$

We can define more accurate relations by utilizing more points. For two points:

**Forward $O(h^2)$:**
$$f'(x_i) = \frac{-f(x_{i+2}) + 4f(x_{i+1}) - 3f(x_i)}{2h}$$

**Backward $O(h^2)$:**
$$f'(x_i) = \frac{3f(x_i) - 4f(x_{i-1}) + f(x_{i-2})}{2h}$$

**Central $O(h^4)$:**
$$f'(x_i) = \frac{-f(x_{i+2}) + 8f(x_{i+1}) - 8f(x_{i-1}) + f(x_{i-2})}{12h}$$

### Richardson Extrapolation
Similar to Romberg's method for integrals, we can iteratively make better derivative predictions by utilizing two previous approximations of different step size $(h_2 < h_1)$. For example, we can combine two $O(h^2)$ approximations (eg. 2-point central FD) to cancel out their errors and obtain an $O(h^4)$ approximation. For step halving $(h_2 = \frac{h_1}{2})$:

$$D \approx \frac{4}{3}D(h_2) - \frac{1}{3}D(h_1) $$

For $j$ consecutive halving steps:

$$D_{i,j} \approx \frac{4^{j-1}D_{i+1,j-1} - D_{i,j-1}}{4^{j-1} - 1}$$
