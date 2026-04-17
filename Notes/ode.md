# Numerical Methods: Ordinary Differential Equations (ODEs)

To solve:

$$\frac{dy}{dx} = f(x, y)$$

with initial condition $y[x_0]=y_0$, we approximate the solution in general as:

$$y_{i+1} = y_i + \phi(x_i, y_i, h)h$$

We call $\phi(x_i, y_i, h)$ the increment function.\
We can generalize to a system of $n$ ODEs:

$$ \frac{dy_1}{dx} = f_1(x, y_1, y_2, ..., y_n) $$
$$ \frac{dy_2}{dx} = f_2(x, y_1, y_2, ..., y_n) $$
$$...$$
$$ \frac{dy_n}{dx} = f_n(x, y_1, y_2, ..., y_n) $$

Now we need $n$ initial conditions:

$$y_1[x_0]=y_{1,0}$$
$$y_2[x_0]=y_{2,0}$$
$$...$$
$$y_n[x_0]=y_{n,0}$$

## Explicit Methods
Explicit methods solve for $y_{i+1}$ using information from previous steps only.

### Runge-Kutta (RK) Methods
The increment function is for $s$ stages per step:

$$ \phi(x_i, y_i, h) = b_1k_1+b_2k_2+...+b_sk_s$$

$k_j$ and $b_j$ are constants:

$$ k_1 = f(x_i,\; y_i) $$
$$ k_2 = f(x_i+p_2h,\; y_i+q_{21}k_1h) $$
$$ k_3 = f(x_i+p_3h,\; y_i+q_{31}k_1h+q_{32}k_2h) $$
$$ ... $$
$$ k_s = f(x_i+p_{s}h,\; y_i+q_{s,1}k_1h+q_{s,2}k_2h+...+q_{s,s-1}k_{s-1}h) $$

#### RK1 - Explicit Euler
The increment function is $\phi(x_i, y_i, h) = f(x_i, y_i) = k_1$. Thus:

$$y_{i+1} = y_i + f(x_i, y_i)h$$

The method requires one function evaluation per step. It is exact for linear $f(x,y)$.
The local error approximation is $O(h^2)$ and the global $O(h)$.\
The method can be improved by utilizing for terms from the Taylor expansion:

$$y_{i+1} = y_i + f(x_i, y_i)h + \frac{f'(x_i,y_i)}{2}h^2 + O(h^3)$$

#### RK2 methods
The increment function is $\phi(x_i, y_i, h) = b_1k_1 + b_2k_2$.

$$ k_1 = f(x_i,\; y_i) $$
$$ k_2 = f(x_i+p_2h,\; y_i+q_{21}k_2h) \to
k_2 = f(x_i, y_i) + p_2h\frac{\partial f}{\partial x} + q_{21}k_1h\frac{\partial f}{\partial y} $$

These lead to:

$$ y_{i+1} = y_i + (b_1+b_2)f(x_i,y_i)h +
[b_2p_2\frac{\partial f}{\partial x} + b_2q_{21}k_1\frac{\partial f}{\partial y}]h^2 $$

To find the constants $b_j, k_j$ we compare these formulas with the second order Taylor expansion:

$$ y_{i+1} = y_i + f(x_i,y_i)h + \frac{f'(x_i,y_i)}{2}h^2 \to y_{i+1} = y_i + f(x_i,y_i)h + (\frac{\partial f}{\partial x} + \frac{\partial f}{\partial y}\frac{\partial y}{\partial x})_{(x_i,y_i)} \frac{h^2}{2} $$

Also, $\frac{\partial y}{\partial x}(x_i,y_i) = f(x_i,y_i) = k_1$. Thus we get the following equations:

$$ b_1 + b_2 = 1 $$
$$ b_2p_2 = \frac{1}{2} $$
$$ b_2q_{21} = \frac{1}{2} $$

We have 3 equations and 4 unknowns. We set one unknown and calculate the rest. This leads to a specific RK2 method.\
**Heun's method**: For $b_2 = \frac{1}{2}$:

$$ b_1 = \frac{1}{2}, \quad p_2 = q_{21} = 1 $$
$$ y_{i+1} = y_i + \frac{1}{2}(k_1 + k_2)h $$
$$ k_1 = f(x_i, y_i) $$
$$ k_2 = f(x_i+h, y_i+k_1h) $$

**Midpoint method**: For $b_2 = 1$:

$$ b_1 = 0, \quad p_2 = q_{21} = \frac{1}{2} $$
$$ y_{i+1} = y_i + k_2h $$
$$ k_1 = f(x_i, y_i) $$
$$ k_2 = f(x_i+\frac{h}{2}, y_i+\frac{k_1}{2}h) $$

**Ralston's method**: For $b_2 = \frac{2}{3}$:

$$ b_1 = \frac{1}{3}, \quad p_2 = q_{21} = \frac{3}{4} $$
$$ y_{i+1} = y_i + \frac{1}{3} (k_1+2k_2)h $$
$$ k_1 = f(x_i, y_i) $$
$$ k_2 = f(x_i+\frac{3h}{4}, y_i+\frac{3k_1}{4}h) $$

The method requires two function evaluations per step. It is exact for quadratic $f(x,y)$.
The local error approximation is $O(h^3)$ and the global $O(h^2)$.

#### RK3 methods
The increment function is $\phi(x_i, y_i, h) = b_1k_1 + b_2k_2+b_3k_3$.

$$ k_1 = f(x_i,\; y_i) $$
$$ k_2 = f(x_i+p_2h,\; y_i+q_{21}k_1h) $$
$$ k_3 = f(x_i+p_3h,\; y_i+q_{31}k_1h+q_{32}k_2h) $$

To find the constants $b_j, k_j$ we compare these formulas with the third order Taylor expansion:

$$ y_{i+1} = y_i + f(x_i,y_i)h + \frac{f'(x_i,y_i)}{2!}h^2 + \frac{f''(x_i,y_i)}{3!}h^3 $$

Following the same procedure as with RK2 we get 6 equations for 8 unknowns:

$$ b_1 + b_2 + b_3 = 1$$
$$ b_2p_2 + b_3p_3 = \frac{1}{2}$$
$$ b_2q_{21} + b_3(q_{31}+q_{32}) = \frac{1}{2}$$
$$ b_2p_2^2 + b_3p_3^2 = \frac{1}{3}$$
$$ b_3q_{32}p_2 = \frac{1}{6}$$
$$ b_2p_2q_{21} + b_3p_3(q_{31}+q_{32}) = \frac{1}{3}$$

We set two unknowns and calculate the rest. This leads to a specific RK3 method. Setting $p_2, p_3$ we get a $2 \times 2$ system for $b_2, b_3$:

$$ \begin{bmatrix}
p_2 & p_3 \\
p_2^2 & p_3^2
\end{bmatrix}
\begin{bmatrix}
b_2 \\
b_3
\end{bmatrix} =
\begin{bmatrix}
1/2 \\
1/3
\end{bmatrix}$$

Solving for $b_2, b_3$ we, then, calculate:

$$ b_1 = 1 - b_2 - b_3, \quad q_{32} = \frac{1}{6p_2b_3} $$

Then, we have another $2 \times 2$ system for $q_{21}, q_{31}$. Instead of solving this, we usually set:

$$ p_2 = q_{21},\quad p_3 = q_{31} + q_{32} $$

Thus, setting the values of $p_3, p_3$ we obtain different RK3 versions. For the classic RK3 method:

$$ p_2 = \frac{1}{2}, p_3 = 1 $$
$$ b_1 = \frac{1}{6}, b_2 = \frac{4}{6}, b_3 = \frac{1}{6} $$
$$ q_{21} = \frac{1}{2}, q_{31} = -1, q_{32} = 2 $$

$$ y_{i+1} = y_i + \frac{1}{6}(k_1 + 4k_2+k_3)h $$
$$ k_1 = f(x_i,\; y_i) $$
$$ k_2 = f(x_i+\frac{1}{2}h,\; y_i+\frac{k_1}{2}h) $$
$$ k_3 = f(x_i+h,\; y_i-k_1h+2k_2h) $$

The method requires three function evaluations per step. It is exact for cubic $f(x,y)$.
The local error approximation is $O(h^4)$ and the global $O(h^3)$.

#### RK4 methods
The increment function is $\phi(x_i, y_i, h) = b_1k_1 + b_2k_2+b_3k_3+b_4k_4$.

$$ k_1 = f(x_i,\; y_i) $$
$$ k_2 = f(x_i+p_2h,\; y_i+q_{21}k_1h) $$
$$ k_3 = f(x_i+p_3h,\; y_i+q_{31}k_1h+q_{32}k_2h) $$
$$ k_4 = f(x_i+p_4h,\; y_i+q_{41}k_1h+q_{42}k_2h+q_{43}k_3h) $$

To find the constants $b_j, k_j$ we compare these formulas with the fourth order Taylor expansion:

$$ y_{i+1} = y_i + f(x_i,y_i)h + \frac{f'(x_i,y_i)}{2!}h^2 + \frac{f''(x_i,y_i)}{3!}h^3 + \frac{f^{(3)}(x_i,y_i)}{4!}h^4 $$

Following the same procedure we get 11 equations for 13 unknowns. We need to set two unknowns and calculate the rest. However, the system is much more coupled now.
The classic RK4 is:

$$ y_{i+1} = y_i + (k_1 + 2k_2 + 2k_3 + k_4)\frac{h}{6} $$
$$ k_1 = f(x_i, y_i) $$
$$ k_2 = f(x_i+\frac{h}{2},\; y_i+\frac{k_1}{2}h) $$
$$ k_3 = f(x_i+\frac{h}{2},\; y_i+\frac{k_2}{2}h) $$
$$ k_4 = f(x_i+h,\; y_i+k_3h) $$

The method requires four function evaluations per step. It is exact up to fourth order $f(x,y)$.
The local error approximation is $O(h^5)$ and the global $O(h^4)$.

#### Higher order RK methods
Butcher's 5th order RK (actually an RK6 method):

$$ y_{i+1} = y_i + (7k_1 + 32k_3 + 12k_4 + 32k_5 + 7k_6)\frac{h}{90} $$
$$ k_1 = f(x_i, y_i) $$
$$ k_2 = f(x_i+\frac{h}{4},\; y_i+\frac{k_1}{4}h) $$
$$ k_3 = f(x_i+\frac{h}{4},\; y_i+\frac{k_1+k_2}{8}h) $$
$$ k_4 = f(x_i+\frac{h}{2},\; y_i+(-\frac{k_2}{2}+k_3)h) $$
$$ k_5 = f(x_i+\frac{3}{4}h,\; y_i+\frac{3k_1+9k_4}{16}h) $$
$$ k_6 = f(x_i+h,\; y_i+\frac{-3k_1+2k_2+12k_3-12k_4+8k_5}{7}h) $$

The method requires six function evaluations per step. The local error approximation is $O(h^6)$ and the global $O(h^5)$.

#### General Butcher Tableau
The generalization of the RK methods for $s$ stages:

$$ y_{i+1} = y_i + h\sum_{i=1}^{s} b_ik_i $$
$$ k_1 = f(x_i,\; y_i) $$
$$ k_2 = f(x_i+p_2h,\; y_i+q_{21}k_1h) $$
$$ k_3 = f(x_i+p_3h,\; y_i+q_{31}k_1h+q_{32}k_2h) $$
$$ ... $$
$$ k_s = f(x_i+p_{s}h,\; y_i+q_{s,1}k_1h+q_{s,2}k_2h+...+q_{s,s-1}k_{s-1}h) $$

can be written in a table of $p_i$ (nodes), $q_{ij}$ (coefficients), $b_i$ (weights):

$$ \begin{array}{c|ccccc}
0      &        &        &        &        & \\
p_2    & q_{21} &        &        &        & \\
p_3    & q_{31} & q_{32} &        &        & \\
\vdots & \vdots & \vdots & \ddots &        & \\
p_s    & q_{s1} & q_{s2} & \dots  & q_{s,s-1} & \\
\hline
       & b_1    & b_2    & \dots  & b_{s-1} & b_s \\
\end{array} $$

For each stage $i$, the sum of the weights must match the time offset:

$$\sum_{j=1}^{i-1} q_{ij} = p_i$$

The weights must sum to 1 to ensure the method is consistent:

$$\sum_{i=1}^{s} b_i = 1$$

The Butcher tableaus for common methods are:

**RK1 - Explicit Euler**:

$$ \begin{array}{c|ccccc}
0 \\
\hline
       & 1
\end{array} $$

**RK2**: \
Heun's:
$$ \begin{array}{c|ccccc}
0 &  \\
1 & 1 \\
\hline
& 1/2    & 1/2 \\
\end{array} $$

Midpoint:
$$ \begin{array}{c|ccccc}
0   &    \\
1/2 & 1/2 \\
\hline
& 0 & 1 \\
\end{array} $$

Ralston's:
$$ \begin{array}{c|ccccc}
0   &    \\
3/4 & 3/4 \\
\hline
& 1/3 & 2/3 \\
\end{array} $$

**Classic RK3**:
$$ \begin{array}{c|ccccc}
0    &    \\
1/2  & 1/2 \\
1    & -1 & 2 \\
\hline
    & 1/6 & 4/6 & 1/6 \\
\end{array} $$

**Classic RK4**:
$$ \begin{array}{c|ccccc}
0   &     &     &   & \\
1/2 & 1/2 &     &   & \\
1/2 & 0   & 1/2 &   & \\
1   & 0   & 0   & 1 & \\
\hline
    & 1/6 & 2/6 & 2/6 & 1/6 \\
\end{array} $$

#### Adaptive RK methods
When the solution of the ODE is abrupt in some region we need to adaptively control the step size $h$.
To do this, we calculate the error between two consecutive predictions.
These predictions are either with the same RK method, but different step $h$, or with the same step $h$ but different method.\
**RK45 Fehlberg**: We estimate $y_{i+1}$ with a 4th and a 5th order method. The algorithm requires 6 function evaluations per step ($k_1 - k_6$). The Butcher tableau is:

$$ \begin{array}{c|ccccc}
0 \\
1/4 & 1/4 \\
3/8 & 3/32 & 9/32 \\
12/13 &	1932/2197 &	−7200/2197 & 7296/2197 \\
1 &	439/216 & −8 & 3680/513 & -845/4104 \\
1/2 & −8/27 & 2 & −3544/2565 & 1859/4104 & −11/40 \\
\hline
& 16/135 & 0 & 6656/12825 & 28561/56430 &	−9/50 &	2/55 \\
& 25/216 & 0 & 1408/2565 & 2197/4104 & −1/5 &	0 \\
\end{array} $$

**RK45 Cash-Karp**: We estimate $y_{i+1}$ with a 4th and a 5th order method. The algorithm requires 6 function evaluations per step ($k_1 - k_6$). The Butcher tableau is:

$$\begin{array}{c|cccccc}
0 & & & & & & \\
1/5 & 1/5 & & & & & \\
3/10 & 3/40 & 9/40 & & & & \\
3/5 & 3/10 & -9/10 & 6/5 & & & \\
1 & -11/54 & 5/2 & -70/27 & 35/27 & & \\
7/8 & 1631/55296 & 175/512 & 575/13824 & 44275/110592 & 253/4096 & \\
\hline
& 37/378 & 0 & 250/621 & 125/594 & 0 & 512/1771 \\
& 2825/27648 & 0 & 18575/48384 & 13525/55296 & 277/14336 & 1/4
\end{array}$$

**RK45 Dormand-Prince**: We estimate $y_{i+1}$ with a 4th and a 5th order method. The algorithm requires 6 function evaluations per step ($k_1 - k_7$, but $k_1$ on the next iteration is $k_7$ of the previous). The Butcher tableau is:

$$\begin{array}{c|ccccccc}
0 & & & & & & & \\
1/5 & 1/5 & & & & & & \\
3/10 & 3/40 & 9/40 & & & & & \\
4/5 & 44/45 & -56/15 & 32/9 & & & & \\
8/9 & 19372/6561 & -25360/2187 & 64448/6561 & -212/729 & & & \\
1 & 9017/3168 & -355/33 & 46732/5247 & 49/176 & -5103/18656 & & \\
1 & 35/384 & 0 & 500/1113 & 125/192 & -2187/6784 & 11/84 & \\
\hline
& 35/384 & 0 & 500/1113 & 125/192 & -2187/6784 & 11/84 & 0 \\
& 5179/57600 & 0 & 7571/16695 & 393/640 & -92097/339200 & 187/2100 & 1/40
\end{array}$$

In all cases, we calculate the error (difference) between the two predictions and update the step:

$$err = max(\; abs(\; y(O(h^5)) - y(O(h^4)) \;) \;)$$
$$ h_{new} = h (\frac{eps}{err})^a $$

where $eps$ is a pre-defined tolerance and $a \approx 0.2$ a parameter.

### Implicit Methods
For stiff (systems of) ODEs explicit methods are unstable unless very small step $h$ is used.
Implicit methods remedy this problem by evaluating the derivative at the current step instead of the previous one.
Now, we have to solve a system of (non-linear) equations for the slopes $k_j$ or intermediate values ($Y_j$) the at each step:

$$ y_{i+1} = y_i + h\sum_{j=1}^{n_s} b_jk_j $$
$$ k_j = f(t_i + p_j h, Y_j) $$
$$ Y_j = y_i + h \sum_{m=1}^{n_s} q_{jm} f(t_i + p_m h, Y_m) = y_i + h \sum_{m=1}^{n_s} q_{jm} k_m $$

**Implicit (Backward) Euler**:
The simplest single-stage ($n_s=1$) method. The Butcher tableau is:

$$ \begin{array}{c|ccccc}
1 & 1 \\
\hline
       & 1
\end{array} $$

So:

$$y_{i+1} = y_i + f(x_{i+1}, y_{i+1})h$$

The method is unconditionally stable.
The local error approximation is $O(h^2)$ and the global $O(h)$.

**Crank-Nicolson (Trapezoidal)**:
This is a hybrid between an implicit and an explicit method.
It combines the stability of implicit Euler and it is also second order accurate (The local error approximation is $O(h^3)$ and the global $O(h^2)$).
The Butcher tableau is:

$$ \begin{array}{c|ccccc}
0 & 0 & 0 \\
1 & 1/2 & 1/2 \\
\hline
    & 1/2 & 1/2
\end{array} $$

And:

$$y_{i+1} = y_i + \frac{h}{2}[f(x_{i}, y_{i})+f(x_{i+1}, y_{i+1})]$$

**Implicit Midpoint**:
The Butcher tableau is:

$$ \begin{array}{c|ccccc}
1/2 & 1/2 \\
\hline
    & 1
\end{array} $$

And:

$$y_{i+1} = y_i + f(x_{i}+\frac{h}{2}, y_{i}+\frac{h}{2}y_{i+1})h$$

**Gauss-Legendre-2**:
Two-stage method ($n_s=2$) that is $O(h^4)$ accurate. The Butcher tableau is:

$$ \begin{array}{c|ccccc}
\frac{1}{2}-\gamma & 1/4 & 1/4-\gamma \\
\frac{1}{2}+\gamma & 1/4+\gamma & 1/4 \\
\hline
    & 1/2 & 1/2
\end{array} $$

where $\gamma = \frac{\sqrt{3}}{6}$.

**Gauss-Legendre-3**:
Three-stage method ($n_s=3$) that is $O(h^6)$ accurate. The Butcher tableau is:

$$ \begin{array}{c|ccccc}
\frac{1}{2}-\frac{\sqrt{15}}{10} & 5/36 & 2/9-\frac{\sqrt{15}}{10} & 5/36-\frac{\sqrt{15}}{30} \\
\frac{1}{2} & 5/36+\frac{\sqrt{15}}{24} & 2/9 & 5/36-\frac{\sqrt{15}}{24} \\
\frac{1}{2}+\sqrt{15} & 5/36+\frac{\sqrt{15}}{30} & 2/9+\frac{\sqrt{15}}{15} & 5/36 \\
\hline
    & 5/18 & 8/18 & 5/18
\end{array} $$

**Randau-IIA-2**:
Two-stage method ($n_s=2$) that is $O(h^3)$ accurate.
It is L-stable for stiffer equations.
The Butcher tableau is:

$$ \begin{array}{c|cc}
1/3 & 5/12 & -1/12 \\
1 & 3/4 & 1/4 \\
\hline
    & 3/4 & 1/4
\end{array} $$

**Randau-IIA-3**:
Three-stage method ($n_s=3$) that is $O(h^5)$ accurate.
It is L-stable for stiffer equations.
The Butcher tableau is:
$$\begin{array}{c|ccc}
\frac{4-\sqrt{6}}{10} & \frac{88-7\sqrt{6}}{360} & \frac{296-169\sqrt{6}}{1800} & \frac{-2+3\sqrt{6}}{225} \\
\frac{4+\sqrt{6}}{10} & \frac{296+169\sqrt{6}}{1800} & \frac{88+7\sqrt{6}}{360} & \frac{-2-3\sqrt{6}}{225} \\
1 & \frac{16-\sqrt{6}}{36} & \frac{16+\sqrt{6}}{36} & \frac{1}{9} \\
\hline
& \frac{16-\sqrt{6}}{36} & \frac{16+\sqrt{6}}{36} & \frac{1}{9}
\end{array}$$

**Lobatto-IIIC-2**:
Two-stage method ($n_s=2$) that is $O(h^2)$ accurate.
It is L-stable for stiffer equations.
The Butcher tableau is:

$$\begin{array}{c|cc}
0 & 1/2 & -1/2 \\
1 & 1/2 & 1/2 \\
\hline
& 1/2 & 1/2
\end{array}$$

**Lobatto-IIIA-3**:
Three-stage method ($n_s=3$) that is $O(h^4)$ accurate.
It is L-stable for stiffer equations.
The Butcher tableau is:

$$\begin{array}{c|ccc}
0 & 0 & 0 & 0 \\
1/2 & 5/24 & 1/3 & -1/24 \\
1 & 1/6 & 2/3 & 1/6 \\
\hline
& 1/6 & 2/3 & 1/6
\end{array}$$

### Multi-step Methods
We use information from previous points apart from $i$.

**Non-self starting Heun**:
Heun's predictor:

$$ y_{i+1}^0 = y_i + f(x_i, y_i)h $$

is of order $O(h^2)$. The corrector:

$$ y_{i+1} = y_i + \frac{f(x_i, y_i)+f(x_{i+1}, y_{i+1}^0)}{2}h $$

is of order $O(h^3)$. The lower order predictor can be improved to $O(h^3)$ by using a previous point $i-1$:

$$ y_{i+1}^0 = y_{i-1} + f(x_i, y_i)2h $$

The method is non-self starting because we need to specify the value $y_{-1}$ as well as $y_0$.\
The corrector step (as in the regular Heun's method) can be repeated for a better prediction.\
The predictor is derived from the open midpoint integration method which has local error $O(h^3)$.
The step error can be expressed in terms of the predictor, $y_{i+1}^0$ and the corrector, $y_{i+1}^m$ (using $m$ correction steps):

$$ e_{i+1} = -\frac{y_{i+1}^0-y_{i+1}^m}{5} $$

The error can be used as a corrector modifier:

$$  y_{i+1}^m \leftarrow  y_{i+1}^m-\frac{y_{i+1}^0-y_{i+1}^m}{5} $$

The same can be done for the predictor:

$$  y_{i+1}^0 \leftarrow  y_{i+1}^0+\frac{4}{5}(y_i^m-y_i^0) $$

The method can be improved by using a higher order open Newton-Cotes (usually for the predictor):

$$ n=1: y_{i+1} = y_{i-1} + 2hf_i $$
$$ n=2: y_{i+1} = y_{i-2} + \frac{3h}{2}(f_i+f_{i-1}) $$
*Milne's method*:
$$ n=3: y_{i+1} = y_{i-3} + \frac{4h}{3}(2f_i-f_{i-1}+2f_{i-2}) $$

And the corresponding close Newton-Cotes for the corrector.

**Adams-Bashforth**: Explicit methods can be derived by Taylor expansion:

$$ y_{i+1} = y_i + f_ih + \frac{f'_i}{2}h^2
+ \frac{f''_i}{6}h^3 + ... = y_i + h(f_i + \frac{f'_i}{2}h + \frac{f''_i}{6}h^2 + ...) $$

Approximating $n$ derivatives with forward FD:

$$ f'_i = \frac{f_i - f_{i-1}}{h} $$
$$ ... $$

we obtain the $n^{th}$ order Adams-Bashforth formulas:

$$ y_{i+1} = y_i + h\sum_{k=0}^{n-1}\beta _k f_{i-k} + O(h^{n+1})$$

The coefficients $\beta _k$:

$$\begin{array}{cccccccc}
\hline
\text{Order} & \beta_0 & \beta_1 & \beta_2 & \beta_3 & \beta_4 & \beta_5 & \text{Local Truncation Error} \\
\hline
1 & 1 & & & & & & \frac{1}{2} h^2 f'(\xi) \\
2 & 3/2 & -1/2 & & & & & \frac{5}{12} h^3 f''(\xi) \\
3 & 23/12 & -16/12 & 5/12 & & & & \frac{9}{24} h^4 f^{(3)}(\xi) \\
4 & 55/24 & -59/24 & 37/24 & -9/24 & & & \frac{251}{720} h^5 f^{(4)}(\xi) \\
5 & 1901/720 & -2774/720 & 2616/720 & -1274/720 & 251/720 & & \frac{475}{1440} h^6 f^{(5)}(\xi) \\
6 & 4277/720 & -7923/720 & 9982/720 & -7298/720 & 2877/720 & -475/720 & \frac{19,087}{60,480} h^7 f^{(6)}(\xi) \\
\hline
\end{array}$$

**Adams-Moulton**: Implicit methods can be derived by backward Taylor expansion around $x_{i+1}$:

$$ y_i = y_{i+1} + f_{i+1}h + \frac{f'_{i+1}}{2}h^2
+ \frac{f''_{i+1}}{6}h^3 + ... \to $$
$$ y_i =  y_{i+1} + h(f_{i+1} + \frac{f'_{i+1}}{2}h + \frac{f''_{i+1}}{6}h^2 + ...) $$

Approximating $n$ derivatives with forward FD:

$$ f'_i = \frac{f_i - f_{i-1}}{h} $$
$$ ... $$

we obtain the $n^{th}$ order Adams-Moulton formulas:

$$ y_{i+1} = y_i + h\sum_{k=-1}^{n-1}\beta _k f_{i-k} + O(h^{n+1}) \to $$
$$ y_{i+1} = y_i + h\beta _{-1}f_{i+1} + h\sum_{k=0}^{n-1}\beta _k f_{i-k} + O(h^{n+1}) $$

The summation starts from $k=-1$ which corresponds to the implicit term $f_{i+1}$.
The coefficients $\beta _k$:

$$\begin{array}{cccccccc}
\hline
\text{Order} & \beta_{-1} & \beta_0 & \beta_1 & \beta_2 & \beta_3 & \beta_4 & \text{Local Truncation Error} \\
\hline
2 & 1/2 & 1/2 & & & & & -\frac{1}{12} h^3 f''(\xi) \\
3 & 5/12 & 8/12 & -1/12 & & & & -\frac{1}{24} h^4 f^{(3)}(\xi) \\
4 & 9/24 & 19/24 & -5/24 & 1/24 & & & -\frac{19}{720} h^5 f^{(4)}(\xi) \\
5 & 251/720 & 646/720 & -264/720 & 106/720 & -19/720 & & -\frac{27}{1440} h^6 f^{(5)}(\xi) \\
6 & 475/1440 & 1427/1440 & -798/1440 & 482/1440 & -173/1440 & 27/1440 & -\frac{863}{60,480} h^7 f^{(6)}(\xi) \\
\hline
\end{array}$$
