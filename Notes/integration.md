# Numerical Methods: Integration

## Integration
The goal is to integrate:

$$I = \int_{}f(\underline{x}) d\underline{x}$$

$f(\underline{x})$ can be an analytical function or a set of data points. $\underline{x}$ is the vector of dependent variables. For 1D integration we have:

$$I = \int_{a}^{b} f(x) dx$$

### Newton-Cotes
Newton-Cotes formulas replace a complicated function with a simpler (polynomial) approximation:

$$I = \int_{a}^{b} f(x) dx \approx \int_{a}^{b} f_n(x) dx$$

${n}$ is the order of the polynomial approximation.

**Trapezoidal Rule**: Uses a linear polynomial as approximation. It integrates exactly a linear function.

$$f_1(x) = f(a) + \frac{f(b) - f(a)}{b - a}(x - a)$$

*   *Integral Approximation:*

$$I \approx (b - a) \frac{f(a) + f(b)}{2}$$

*   *Error:*

$$E_a = -\frac{(b-a)^3}{12}f''(\xi) = -\frac{1}{12}f''(\xi) h^3$$

where $h = b - a$ and $\xi \in [a, b]$.

**Composite Trapezoidal Rule**: We divide the domain into $n$ equally spaced segments and apply the Trapezoidal rule. The segment length is:

$$h = \frac{b - a}{n}$$

*   *Integral approximation:*

$$I = \frac{h}{2} \left[ f(x_0) + 2 \sum_{i=1}^{n-1} f(x_i) + f(x_n) \right]$$

*   *Error:*

$$E_a = -\frac{(b - a)^3}{12n^2} f''=-\frac{(b - a)}{12}h^2f''$$

**Simpson 1/3**: Uses a second order polynomial as approximation. It integrates exactly a cubic function.

*   *Integral Approximation:*

$$I \approx \frac{b - a}{6} [f(x_0) + 4f(x_1) + f(x_2)]$$

*   *Error:*

$$E_a = -\frac{(b-a)^5}{2880} f^{(4)}(\xi) = -\frac{1}{90} f^{(4)}(\xi)h^5$$

where $h = \frac{b-a}{2}$ and $\xi\in[a, b]$. The midpoint is $x_1 = \frac{b+a}{2}$.

**Composite Simpson 1/3**: We divide the domain into $n$ even and equally spaced segments. The segment length is:

$$h = \frac{b - a}{n}$$

*   *Integral approximation:*

$$I = \frac{h}{3} \left[ f(x_0) + 4 \sum_{i=1, odd}^{n-1} f(x_i) + 2 \sum_{j=2, even}^{n-2} f(x_j) + f(x_n) \right]$$

*   *Error:*

$$E_a = -\frac{(b - a)^5}{180n^4}\overline{f}^{(4)}=-\frac{(b - a)}{180}\overline{f}^{(4)}h^4$$

**Simpson 3/8**: Uses a third order polynomial as approximation. It integrates exactly a cubic function. It practically has the same accuracy as Simpson 1/3. For this reason it is only used at the last segment of an odd number of segments, while for the rest of the segments we use Simpson 1/3.

*   *Integral Approximation:*
$$I \approx \frac{b - a}{8} [f(x_0) + 3f(x_1) + 3f(x_2) + f(x_3)]$$

*   *Error:*

$$E_a = -\frac{(b-a)^5}{6480} f^{(4)}(\xi) = -\frac{3}{80} f^{(4)}(\xi)h^5$$

where $h = \frac{b-a}{3}$ and $\xi \in [a, b]$. The midpoints are $x_1 = a+h$, $x_2 = a+2h$.

#### Uneven segments
When we want to integrate datasets segments are usually uneven. We can modify the Newton Cotes formulas to make the integration with:

$$ h_i = x_{i+1} - x_i $$

**Trapezoidal rule**:

$$I \approx \sum_{i=1}^{n-1} h_i [f(x_i)+f(x_{i+1})]$$

**Simpson 1/3**:

$$I \approx \frac{1}{6}\sum_{i=1}^{n-1} (h_i+h_{i+1})
[(2-\frac{h_{i+1}}{h_i})f(x_i) + \frac{(h_{i+1}+h_i)^2}{h_{i+1}h_i}f(x_{i+1}) + (2-\frac{h_i}{h_{i+1}})f(x_{i+2})]$$

### Romberg Integration
We estimate the integral with different step sizes ($h_2 < h_1$). We use these estimations to cancel out errors and get a better prediction. For example, we can combine two $O(h^2)$ to cancel out their errors and obtain an $O(h^4)$ approximation. The process can be repeated by progressively decreasing $h$.

$$ I = I(h_1)+E(h_1) = I(h_2)+E(h_2)$$

**Trapezoidal rule**: The errors are $E(h)\sim O(h^2)$. Thus:

$$ I(h_1)+E(h_2)\frac{h_1^2}{h_2^2} = I(h_2)+E(h_2) \to E(h_2) \approx \frac{I(h_2)-I(h_1)}{\frac{h_1^2}{h_2^2}-1}$$

The improved integral estimation is:

$$ I = I(h_2)+E(h_2) \approx I(h_2)+\frac{I(h_2)-I(h_1)}{\frac{h_1^2}{h_2^2}-1}$$

For step halving $(h_2 = \frac{h_1}{2})$:

$$ I \approx I(h_2)+\frac{I(h_2)-I(h_1)}{4-1} \to
I \approx \frac{4}{3}I(h_2)-\frac{1}{3}I(h_1)$$

For $j$ consecutive halving steps:

$$I_{i,j} \approx \frac{4^{j-1}I_{i+1,j-1} - I_{i,j-1}}{4^{j-1} - 1}$$

### Gauss Quadrature
We find the optimal values for the points and weights of integration to minimize the error of integrating certain basis functions.

#### Gauss-Legendre
For Gauss-Legendre quadrature, the basis functions are polynomials, i.e. we seek to accurately integrate polynomials.

**General Formula**

$$I \approx \sum_{i=1}^{n} w_i f(x_i)$$

For $n$ points of integration we evaluate the $2n$ constants by requiring that the formula accurately integrates polynomials up to order $2n-1$.

$$c_11+c_21 + ... + c_n1 = \int_{-1}^{1}1dx = 2$$
$$c_1\xi_1+c_2\xi_2 + ... + c_n\xi_n = \int_{-1}^{1}xdx = 0$$
$$c_1\xi_1^2+c_2\xi_2^2 + ... + c_n\xi_n^2 = \int_{-1}^{1}x^2dx = \frac{2}{3}$$
$$...$$
$$c_1\xi_1^{2n-1}+c_2\xi_2^{2n-1} + ... + c_n\xi_n^{2n-1} = \int_{-1}^{1}x^{2n-1}dx $$

The limits of integrations are $[-1, 1]$ for ease of integration. We can rewrite any integral with limits $[a, b]$ by a change of variables:

$$ x = \frac{b-a}{2}\xi + \frac{b+a}{2} \to x = m\xi + d$$
$$ dx = \frac{b-a}{2}d\xi \to dx = md\xi$$

Thus:

$$ I =\int_{a}^{b}f(x)dx = \int_{-1}^{1}f(\xi)md\xi \approx m\sum_{i = 1}^{n} w_if(\xi_i)$$

The $n\times n$ non-linear system can be solved for $(c_i, \xi_i)$. However, for large $n$ it becomes singular. A better way to find the coefficients is to find the roots of the Legendre polynomial, $P_n(x)$ Using Bonnet’s recursion:

$$(n+1)P_{n+1}(x) = (2n+1)xP_n(x) - nP_{n-1}(x)$$
$$P'_n(x) = \frac{n}{x^2-1}[xP_n(x)-P_{n-1}(x)]$$

We can find the roots of $P_n(x)$ with Newton's method. A good initial guess is:

$$x_i = \cos(\pi\frac{i-0.25}{n-0.25})$$

After finding $x_i$ we can calculate the weights analytically:

$$w_i = \frac{2}{(1-x_i^2)[P'_n(x_i)]^2}$$

In multiple dimensions we follow the same procedure for each dimension. For example, in 2D:

$$ I =\int_{a_1}^{b_1} \int_{a_2}^{b_2} f(x)dx =
\int_{-1}^{1} \int_{-1}^{1} f(\xi,\eta)m_1m_2d\xi d\eta \approx m_1m_2\sum_{i = 1}^{n_1} \sum_{j = 1}^{n_2} w_iw_jf(\xi_i,\eta_i)$$
