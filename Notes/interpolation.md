# Numerical Methods: Interpolation

In interpolation, we are interested in the accurate description of precise data and their intermediate values.
An $n^{th}$ order polynomial can accurately fit $n + 1$ data points.
In all cases, it is best to sort the data before proceeding with the interpolation.

## Vandermonde
We can calculate the $n + 1$ coefficients of a polynomial $f(x) = a_0 + a_1x + \dots + a_nx^n$ that passes through all $n + 1$ data points:

$$f(x_0) = a_0 + a_1x_0 + \dots + a_nx_0^n$$
$$f(x_1) = a_0 + a_1x_1 + \dots + a_nx_1^n$$
$$\dots$$
$$f(x_n) = a_0 + a_1x_n + \dots + a_nx_n^n$$

This is a linear system of equations:

$$\begin{bmatrix}
1 & x_0 & \dots & x_0^n \\
1 & x_1 & \dots & x_1^n \\
\dots & \dots & \dots & \dots \\
1 & x_n & \dots & x_n^n
\end{bmatrix} \cdot
\begin{bmatrix}
a_0 \\
a_1 \\
\dots \\
a_n
\end{bmatrix} =
\begin{bmatrix} f(x_0) \\
f(x_1) \\
\dots \\
f(x_n)
\end{bmatrix}$$

The method is generally unstable and slow.

## Gregory-Newton
Newton’s divided-difference interpolating polynomials are usually used when the order of interpolation is unknown, and we want to examine multiple orders to find the best approximation.
Linear polynomial:
To accurately describe a linear function given two points $(x_0, f(x_0))$ and $(x_1, f(x_1))$:

$$f_1(x) = f(x_0) + \frac{f(x_1) - f(x_0)}{x_1 - x_0}(x - x_0) = b_0 + b_1(x - x_0)$$

For non-linear functions, the prediction at $x$ is more accurate the closer the interval $x_1 - x_0$.
Quadratic polynomial:

$$f_2(x) = b_0 + b_1(x - x_0) + b_2(x - x_0)(x - x_1)$$

where

$$b_0 = f(x_0)$$
$$b_1 = \frac{f(x_1) - f(x_0)}{x_1 - x_0} = fd(x_1, x_0)$$
$$b_2 = \frac{\frac{f(x_2) - f(x_1)}{x_2 - x_1} - \frac{f(x_1) - f(x_0)}{x_1 - x_0}}{x_2 - x_0} = fd(x_2, x_1, x_0)$$

Here, $fd$ denotes the finite difference approximation.General $n^{th}$ order polynomial

$$f_n(x) = f(x_0) + fd(x_1, x_0)(x - x_0) + fd(x_2, x_1, x_0)(x - x_0)(x - x_1) + \dots + fd(x_n, \dots, x_0)(x - x_0) \dots (x - x_{n-1})$$

The method requires storing the finite differences at every data point.\
Error Estimation: The error estimate of the prediction can be calculated if we have an extra data point: $e = f_{n+1}(x) - f_n(x)$.
This estimate is converging fast and can be less than the true error, which is not desirable. Moreover, high-order polynomials are ill-conditioned and prone to error.
Consequently, as the order of interpolation increases, there is a point of diminishing returns.

## Lagrange Polynomials
Lagrange polynomial interpolation is a reformulation of Newton-Gregory polynomials to avoid calculating finite differences.
It is easier to implement and it is used when the order of interpolation is known a priori.

$$f_n(x) = \sum_{i=0}^{n} L_i(x)f(x_i)$$

where the Lagrange coefficients are:

$$L_i(x) = \prod_{j=0, j \neq i}^{n} \frac{x - x_j}{x_i - x_j}$$

The Lagrange coefficients $L_i(x)$ are zero at every data point except $x = x_i$ where they are equal to $1$.
These formulas are for fitting an $n^{th}$ order polynomial to $n + 1$ data points.
If we want to fit with a given polynomial order, $m$, we must use $m + 1$ data points.
We choose the points that are closest to the point of interest $x$.

## Splines
An alternative to high order polynomials that interpolate many data points, is to use lower order polynomials in subsets of the data points.

### Linear Splines
We can use linear splines where a straight line approximates each interval.
The slope of the function is then discontinuous between data points.

$$f(x) = f(x_i) + m_i(x - x_i)$$

where $m_i = \frac{f(x_{i+1}) - f(x_i)}{x_{i+1} - x_i}$ is the slope at interval $i$.
This is the only unknown at each interval.

### Quadratic Splines
To preserve $1^{\text{st}}$ derivative continuity, we can use quadratic splines.
Now we have $n + 1$ data points and 3 unknowns ($a_i, b_i, c_i$) per interval.
So, the total unknowns for $n$ intervals are $3n$.
The spline and its derivatives are:

$$f_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2, x_i \leq x \leq x_{i+1}$$
$$f_i'(x) = b_i + 2c_i(x - x_i)$$
$$f_i''(x) = 2c_i$$

For the spline of interval $i$ at point $x_i$, we have $f_i(x_i) = a_i = y_i$ for $i = 1, \dots, n$.
Thus, in this formulation, the coefficients $a_i$ are known and equal to the dependent variable.
This reduces the total unknowns to $2n$.
We have $n - 1$ conditions of function continuity at interior points, $i = 1, \dots, n - 1$:

$$f_i(x_{i+1}) = f_{i+1}(x_{i+1}) \to $$
$$y_i + b_i(x_{i+1} - x_i) + c_i(x_{i+1} - x_i)^2 = y_{i+1}$$

Also, at the last point $x_{n+1}$:

$$f_n(x_{n+1}) = y_n + b_n(x_{n+1} - x_n) + c_n(x_{n+1} - x_n)^2 = y_{n+1}$$

We also have $n - 1$ conditions of derivative continuity at interior points, $i = 1, \dots, n - 1$:

$$f_i'(x_{i+1}) = f_{i+1}'(x_{i+1}) \to $$
$$b_i + 2c_i(x_{i+1} - x_i) = b_{i+1} \to c_i = \frac{b_{i+1} - b_i}{2(x_{i+1} - x_i)}$$

Inserting this into the function continuity equation:

$$y_i + b_i(x_{i+1} - x_i) + \frac{b_{i+1} - b_i}{2}(x_{i+1} - x_i) = y_{i+1} \to $$
$$\frac{b_i + b_{i+1}}{2}(x_{i+1} - x_i) = y_{i+1} - y_i \to b_i + b_{i+1} = 2\frac{y_{i+1} - y_i}{x_{i+1} - x_i}$$

Now we have totally $2n - 1$ equations.
We usually set a linear start, i.e., the $2^{\text{nd}}$ derivative to 0 at the first point:

$$f_1''(x_1) = 2c_1 = 0 \to b_1 - b_2 = 0$$

Other options are:
*   *Clamped start*:

$$f_1'(x_1) = b_1 = v$$

*   *Average slope*: to follow the trend of the first two points.
$$f_1'(x_1) = b_1 = \frac{y_2 - y_1}{x_2 - x_1}$$

*   *Equal $2^{\text{nd}}$ derivatives*:

$$c_1 = c_2 \to \frac{b_2 - b_1}{(x_2 - x_1)} = \frac{b_3 - b_2}{(x_3 - x_2)}$$

though this breaks the tri-diagonality of the resulting matrix.\
Having constructed a system of $n \times n$ equations, we solve for the coefficients $b_i$ using the Thomas algorithm and then calculate:

$$c_i = \frac{b_{i+1} - b_i}{2(x_{i+1} - x_i)}$$

For a linear start, the system is:

$$\begin{bmatrix}
1 & -1 & \dots & \dots & \dots \\
1 & 1 & \dots & \dots & \dots \\
\dots & 1 & 1 & \dots & \dots \\
\dots & \dots & \dots & \dots & \dots \\
\dots & \dots & \dots & 1 & 1
\end{bmatrix} \cdot
\begin{bmatrix} b_1 \\
b_2 \\
b_3 \\
\dots \\
b_n
\end{bmatrix} =
\begin{bmatrix} 0 \\
2\frac{y_2 - y_1}{x_2 - x_1} \\
2\frac{y_3 - y_2}{x_3 - x_2} \\
\dots \\
2\frac{y_{n+1} - y_n}{x_{n+1} - x_n}
\end{bmatrix}$$

### Cubic Splines
To ensure continuity of the $2^{\text{nd}}$ derivative, we can use cubic splines.
With $n + 1$ data points, there are 4 unknowns ($a_i, b_i, c_i, d_i$) per interval, totaling $4n$ unknowns for $n$ intervals.
The cubic spline and its first two derivatives are defined as:

$$f_i(x) = a_i + b_i(x - x_i) + c_i(x - x_i)^2 + d_i(x - x_i)^3, \quad x_i \leq x \leq x_{i+1}$$
$$f_i'(x) = b_i + 2c_i(x - x_i) + 3d_i(x - x_i)^2$$
$$f_i''(x) = 2c_i + 6d_i(x - x_i)$$
At each point $x_i$, $f_i(x_i) = a_i = y_i$ for $i = 1, \dots, n$.
Because these coefficients $a_i$ are known and equal to the dependent variable, the total number of unknowns is reduced to $3n$.
To solve for the remaining unknowns, several continuity conditions are applied at interior points ($i = 1, \dots, n - 1$):
The spline must match the data point at the end of each interval:

$$f_i(x_{i+1}) = f_{i+1}(x_{i+1}) \to $$
$$y_i + b_i(x_{i+1} - x_i) + c_i(x_{i+1} - x_i)^2 + d_i(x_{i+1} - x_i)^3 = y_{i+1}$$

At the last point: $x_{n+1}$:

$$f_n(x_{n+1}) = y_n + b_n(x_{n+1} - x_n) + c_n(x_{n+1} - x_n)^2 + d_n(x_{n+1} - x_n)^3 = y_{n+1}$$

The slopes must match at interior points:

$$f_i'(x_{i+1}) = f_{i+1}'(x_{i+1}) \to b_i + 2c_i(x_{i+1} - x_i) + 3d_i(x_{i+1} - x_i)^2 = b_{i+1}$$

The curvature must match at interior points:

$$f_i''(x_{i+1}) = f_{i+1}''(x_{i+1}) \to 2c_i + 6d_i(x_{i+1} - x_i) = 2c_{i+1}$$

From this, we can solve for $d_i$:

$$d_i = \frac{c_{i+1} - c_i}{3(x_{i+1} - x_i)}$$

By substituting the expression for $d_i$ into the continuity conditions, the function and derivative continuity equations become:

$$\begin{cases}
y_i + b_i(x_{i+1} - x_i) + c_i(x_{i+1} - x_i)^2 + \frac{c_{i+1} - c_i}{3}(x_{i+1} - x_i)^2 = y_{i+1} \\
b_i + 2c_i(x_{i+1} - x_i) + (c_{i+1} - c_i)(x_{i+1} - x_i) = b_{i+1}
\end{cases}$$

Which simplifies to:

$$\begin{cases}
y_i + b_i(x_{i+1} - x_i) + \frac{c_{i+1} + 2c_i}{3}(x_{i+1} - x_i)^2 = y_{i+1} \\
b_i + (c_{i+1} + c_i)(x_{i+1} - x_i) = b_{i+1}
\end{cases}$$

Solving the first equation for $b_i$ gives:

$$\begin{cases}
b_i = \frac{y_{i+1} - y_i}{x_{i+1} - x_i} - \frac{c_{i+1} + 2c_i}{3}(x_{i+1} - x_i) \\
b_i + (c_{i+1} + c_i)(x_{i+1} - x_i) = b_{i+1}
\end{cases}$$

We can change the indices $i \to i - 1$ to yield:

$$b_{i-1} = \frac{y_i - y_{i-1}}{x_i - x_{i-1}} - \frac{c_i + 2c_{i-1}}{3}(x_i - x_{i-1})$$
$$b_{i-1} + (c_i + c_{i-1})(x_i - x_{i-1}) = b_i$$

To find the final equation, we combine the formulas for $b_i$ and $b_{i-1}$ by substituting the second index-shifted equation into the first original equation:

$$\frac{y_i - y_{i-1}}{x_i - x_{i-1}} - \frac{c_i + 2c_{i-1}}{3}(x_i - x_{i-1}) + (c_i + c_{i-1})(x_i - x_{i-1}) = \frac{y_{i+1} - y_i}{x_{i+1} - x_i} - \frac{c_{i+1} + 2c_i}{3}(x_{i+1} - x_i)$$

We finally obtain a tri-diagonal system to solve for the $c_i$ coefficients:

$$\frac{(x_i - x_{i-1})}{3}c_{i-1} + \frac{2(x_i - x_{i-1})}{3}c_i + \frac{2(x_{i+1} - x_i)}{3}c_i + \frac{(x_{i+1} - x_i)}{3}c_{i+1} = \frac{y_{i+1} - y_i}{x_{i+1} - x_i} - \frac{y_i - y_{i-1}}{x_i - x_{i-1}}$$

By multiplying both sides of the previous equation by 3, we simplify the tri-diagonal system for the cubic spline coefficients to:

$$(x_i - x_{i-1})c_{i-1} + 2[(x_i - x_{i-1}) + (x_{i+1} - x_i)]c_i + (x_{i+1} - x_i)c_{i+1} = 3 \left( \frac{y_{i+1} - y_i}{x_{i+1} - x_i} - \frac{y_i - y_{i-1}}{x_i - x_{i-1}} \right)$$

To solve this system, we must define boundary conditions at the endpoints $x_1$ and $x_{n+1}$:
*   *Natural Splines*: In a natural spline, the second derivatives at the endpoints are set to zero, which makes the spline a straight line at its ends:

$$f_1''(x_1) = 0 \to c_1 = 0$$
$$f_n''(x_{n+1}) = 0 \to c_n = 0$$

*   *Clamped Splines*: In a clamped spline, the first derivatives at the boundaries are set to specific values $v_1$ and $v_{n+1}$:
At the first point ($x_1$):

$$f_1'(x_1) = v_1 \to b_1 = v_1$$

Substituting the expression for $b_1$, we get:

$$\frac{y_2 - y_1}{x_2 - x_1} - \frac{c_2 + 2c_1}{3}(x_2 - x_1) = v_1 \to $$
$$2c_1 + c_2 = \frac{3}{x_2 - x_1} \left( \frac{y_2 - y_1}{x_2 - x_1} - v_1 \right)$$

At the last point ($x_{n+1}$):

$$f_n'(x_{n+1}) = v_{n+1} \to $$
$$b_n + 2c_n(x_{n+1} - x_n) + 3d_n(x_{n+1} - x_n)^2 = v_{n+1}$$

Substituting the known expressions for $b_n$ and $d_n$:

$$\frac{y_{n+1} - y_n}{x_{n+1} - x_n} - \frac{c_{n+1} + 2c_n}{3}(x_{n+1} - x_n) + 2c_n(x_{n+1} - x_n) + 3 \frac{c_{n+1} - c_n}{3(x_{n+1} - x_n)}(x_{n+1} - x_n)^2 = v_{n+1}$$

Simplifying leads to:

$$\left( \frac{2}{3}c_{n+1} + \frac{1}{3}c_n \right)(x_{n+1} - x_n) = v_{n+1} - \frac{y_{n+1} - y_n}{x_{n+1} - x_n} \to $$
$$2c_{n+1} + c_n = \frac{3}{x_{n+1} - x_n} \left( v_{n+1} - \frac{y_{n+1} - y_n}{x_{n+1} - x_n} \right)$$

## Multi-dimensional Interpolation
We can perform interpolation in multiple dimensions by performing 1D interpolation in one direction while keeping the rest of the coordinates fixed.
Then, we use the new points for interpolating in the remaining directions.

## Fourier Transform
Fourier transform is used when we want to interpolate trigonometric functions, for example for vibrating or oscillating engineering systems.
In Fourier analysis, we deal with both the time and the frequency domains. \
The frequency is: $f = \frac{1}{T} = \frac{\omega_0}{2\pi}$.
A periodic function is one such that: $f(t) = f(t + T)$.
Sinusoids are common periodic functions: $f(t) = A + c_1 \cos(\omega_0t + \theta)$.
They consist of 4 parameters: the mean value, $A_0$, the amplitude, $c_1$, the fundamental frequency, $\omega_0$, and the phase shift, $\theta$.
The presence of the phase shift is problematic for curve fitting. We can reformulate the sinusoid function as:

$$f(t) = A_0 + A_1 \cos(\omega_0t) + B_1 \sin(\omega_0t)$$

where $A_1 = c_1 \cos(\theta)$, $B_1 = -c_1 \sin(\theta)$.

### Linear least-squares
The sinusoid function can be used as a linear least-squares model.
In a regression problem, we minimize:

$$S_r = \sum_{i=1}^{n} \{y_i - [A_0 + A_1 \cos(\omega_0t_i) + B_1 \sin(\omega_0t_i)]\}^2$$

The normal equations that result from $\frac{\partial S_r}{\partial a} = 0$ give the solution:

$$A_0 = \frac{\sum y_i}{n}, A_1 = \frac{2}{n} \sum y_i \cos(\omega_0t_i), B_1 = \frac{2}{n} \sum y_i \sin(\omega_0t_i)$$

For a general model:

$$f(t) = A_0 + \sum_{k=1}^{m} [A_k \cos(k\omega_0t) + B_k \sin(k\omega_0t)]$$

$$A_0 = \frac{\sum y_i}{n}, A_j = \frac{2}{n} \sum y_i \cos(j\omega_0t_i), B_j = \frac{2}{n} \sum y_i \sin(j\omega_0t_i)$$

### Continuous Fourier Series
Fourier showed that an arbitrary periodic function can be represented by an infinite series of sinusoids:

$$f(t) = a_0 + \sum_{k=1}^{\infty} [a_k \cos(k\omega_0t) + b_k \sin(k\omega_0t)] = $$
$$\sum_{k=-\infty}^{\infty} \tilde{c}_k e^{ik\omega_0t}$$

The coefficients are:

$$a_0 = \frac{1}{T} \int_{0}^{T} f(t) dt$$
$$a_k = \frac{2}{T} \int_{0}^{T} f(t) \cos(k\omega_0t) dt$$
$$b_k = \frac{2}{T} \int_{0}^{T} f(t) \sin(k\omega_0t) dt$$
$$\tilde{c}_k = \frac{1}{T} \int_{-T/2}^{T/2} f(t) e^{-ik\omega_0t} dt$$

### Fourier transform
For nonperiodic functions, the Fourier series reduces to (as the period $T \to \infty$):

$$f(t) = \frac{1}{2\pi} \int_{-\infty}^{+\infty} F(i\omega_0) e^{i\omega_0t} d\omega_0$$

This is the inverse Fourier transform of $F(i\omega_0)$ to transform the frequency domain to the time domain.
Conversely, the Fourier transform of $f(t)$ transforms the time domain to the frequency domain:

$$F(i\omega_0) = \frac{1}{2\pi} \int_{-\infty}^{+\infty} f(t) e^{-i\omega_0t} dt$$

The Fourier transform converts a continuous time domain function to a continuous frequency domain function.
In contrast, the Fourier series converts a continuous and periodic time domain function to a discrete frequency domain.

### Discrete Fourier Transform (DFT)
The data from a signal are discrete measurements, rather than a continuous function.
The Fourier Transform for $N$ discrete values of a continuous function is the discrete Fourier Transform:

$$F_k = \sum_{n=0}^{N-1} f_n e^{-ik\omega_0n} = \sum_{n=0}^{N-1} f_n \cos(k\omega_0n) - if_n \sin(k\omega_0n), k = 0, \dots, N-1$$

The discrete inverse Fourier Transform is:

$$f_n = \frac{1}{N} \sum_{k=0}^{N-1} F_k e^{ik\omega_0n} = \frac{1}{N} \sum_{k=0}^{N-1} F_k \cos(k\omega_0n) + iF_k \sin(k\omega_0n), k = 0, \dots, N-1$$

The fundamental frequency is $\omega_0 = \frac{2\pi}{N}$.
DFT requires $N^2$ operations.

### Fast Fourier Transform (FFT)
Alternative DFT algorithms reduce the operations to $O(N\log_2 N)$.
A DFT of length $N$ is decomposed (“decimated”) into successively smaller DFTs.
There are decimation-in-time techniques and decimation-in-frequency techniques.
A well-known decimation-in-time technique is the Sande-Tukey algorithm.
We need the number of data to be an integer power of 2, i.e., $N = 2^M$.
$M$ is an integer that is essentially the number of stages of decimation.
First, the sample is divided in half:

$$F_k = \sum_{n=0}^{\frac{N}{2}-1} f_n e^{-i\frac{2\pi}{N}nk} + \sum_{n=N/2}^{N-1} f_n e^{-i\frac{2\pi}{N}nk}, \quad k = 0, \dots, N-1$$

We define $m = n - \frac{N}{2}$:

$$F_k = \sum_{n=0}^{\frac{N}{2}-1} f_n e^{-i\frac{2\pi}{N}nk} + \sum_{m=0}^{\frac{N}{2}-1} f_{m+\frac{N}{2}} e^{-i\frac{2\pi}{N}k(m+\frac{N}{2})} = $$
$$\sum_{n=0}^{\frac{N}{2}-1} \left(f_n + f_{n+\frac{N}{2}} e^{-i\pi k}\right) e^{-i\frac{2\pi}{N}nk}$$

Using $e^{-i\pi k} = (-1)^k$, we separate the summation into odd and even values of $k$.
For even $k$, $e^{-i\pi k} = 1$:

$$F_{2k} = \sum_{n=0}^{\frac{N}{2}-1} \left(f_n + f_{n+\frac{N}{2}}\right) e^{-i\frac{2\pi}{N}n2k} = $$
$$ \sum_{n=0}^{\frac{N}{2}-1} \left(f_n + f_{n+\frac{N}{2}}\right) e^{-\frac{i2\pi kn}{N/2}}$$

For odd $k$, $e^{-i\pi k} = -1$:

$$F_{2k+1} = \sum_{n=0}^{\frac{N}{2}-1} \left(f_n - f_{n+\frac{N}{2}}\right) e^{-i\frac{2\pi}{N}n(2k+1)} = $$
$$ \sum_{n=0}^{\frac{N}{2}-1} \left(f_n - f_{n+\frac{N}{2}}\right) e^{-\frac{i2\pi kn}{N/2}} e^{-i\frac{2\pi}{N}n}$$

We define $W = e^{-\frac{2\pi}{N}i} = \cos\left(\frac{2\pi}{N}\right) - i \sin\left(\frac{2\pi}{N}\right)$.
Then, for even $k$:

$$F_{2k} = \sum_{n=0}^{\frac{N}{2}-1} \left(f_n + f_{n+\frac{N}{2}}\right) W^{2kn}$$

For odd $k$:

$$F_{2k+1} = \sum_{n=0}^{\frac{N}{2}-1} \left(f_n - f_{n+\frac{N}{2}}\right) W^{2kn}W^n$$

Now a single $N$-point DFT has been replaced by two $\frac{N}{2}$-point DFTs.
Thus, the computations are reduced from $N^2$ to $\frac{N^2}{2}$.
The process is repeated until $\frac{N}{2}$ 2-point DFTs are required.
This reduces the computations to $N \log_2 N$.
The Cooley-Tukey algorithm is the reverse of the Sande-Tukey algorithm, performing decimation in time.

### Power spectrum
The power of a periodic signal in the time domain is:

$$P = \frac{1}{T} \int_{-T/2}^{+T/2} f^2(t) dt$$

Using the Fourier series for $f(t) = \sum_{k=-\infty}^{\infty} F_k e^{ik\omega_0t}$:

$$P = \frac{1}{T} \int_{-T/2}^{+T/2} f^2(t) dt = \sum_{k=-\infty}^{\infty} |F_k^2|$$

The power in the $k^{th}$ real harmonic is:

$$p_k = 2|F_k^2|$$

The power spectrum is the plot $p_k - k\omega_0$.
