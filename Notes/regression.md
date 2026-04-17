# Numerical Methods: Regression

Experimental data is often subject to error, and we cannot use interpolation to pass a function exactly through all points.
In such cases, we perform regression by deriving an approximate function that captures the trend of the data.
We need to establish a criterion to minimize the discrepancy between the data points and the function predictions.
This is done through least-squares regression.
In least-squares regression we minimize the sum (for all n data points) of squares of the residuals between experiment and calculation.

$$S_r = \sum_{i=1}^{n} (y_i - y_{i,model})^2$$

To minimize the sum of squares we set $\frac{\partial S_r}{\partial a_i} = 0$ for all $a_i$ parameters of the model.

## Linear regression
For linear regression the model function in 1D is given by: $y_{model}(x) = a_0 + a_1x$. Thus:

$$S_r = \sum_{i=1}^{n} (y_i - a_0 - a_1x_i)^2$$

$$\frac{\partial S_r}{\partial a_0} = -2 \sum_{i=1}^{n} y_i - a_0 - a_1x_i$$
$$\frac{\partial S_r}{\partial a_1} = -2 \sum_{i=1}^{n} (y_i - a_0 - a_1x_i)x_i$$

Setting the derivatives to zero we get the normal equations:

$$na_0 + a_1 \sum_{i=1}^{n} x_i = \sum_{i=1}^{n} y_i$$
$$a_0 \sum_{i=1}^{n} x_i + a_1 \sum_{i=1}^{n} x_i^2 = \sum_{i=1}^{n} x_iy_i$$

Or:

$$\begin{bmatrix}
n & \sum x_i \\ \sum x_i & \sum x_i^2
\end{bmatrix} \cdot
\begin{bmatrix}
a_0 \\ a_1
\end{bmatrix} =
\begin{bmatrix}
\sum y_i \\ \sum x_iy_i
\end{bmatrix}$$

Solving for $a_0, a_1$ we get:

$$a_0 = \frac{\sum y_i}{n} - a_1 \frac{\sum x_i}{n} = \bar{y} - a_1\bar{x}$$
$$a_1 = \frac{n \sum x_iy_i - \sum x_i \sum y_i}{n \sum x_i^2 - (\sum x_i)^2}$$

We quantify the goodness of the fit by calculating a “standard deviation” called standard error of the estimate

$$s_{y/x} = \sqrt{\frac{S_r}{n-2}}$$

which quantifies the spread of data around the regression line (in contrast to the standard deviation which quantifies the spread of data around the mean value).
The smaller $S_r$ (residual error from regression) is, the smaller the deviation.
We can also calculate the sum of squares around the mean: $S_t = \sum(y_i - \bar{y})^2$ and use it to calculate the coefficient of determination, $r^2 = \frac{S_t - S_r}{S_t}$, i.e. how much the error has reduced by the regression compared to describing the data with the mean value.
$r^2 = 1$ means that the model line can explain 100% of the data variability, while $r^2 = 0$ means that there was no improvement from the regression.

## Polynomial regression
Linear regression can be generalized to fit the data using an $m^{\text{th}}$ order polynomial:

$$y_{model}(x) = a_0 + a_1x + \dots + a_mx^m$$

Thus:

$$S_r = \sum_{i=1}^{n} (y_i - a_0 - a_1x_i - \dots - a_mx_i^m)^2$$
$$\frac{\partial S_r}{\partial a_0} = 2 \sum_{i=1}^{n} (y_i - a_0 - a_1x_i - \dots - a_mx_i^m)(-1)$$
$$\frac{\partial S_r}{\partial a_1} = 2 \sum_{i=1}^{n} (y_i - a_0 - a_1x_i - \dots - a_mx_i^m)(-x_i)$$
$$\dots$$
$$\frac{\partial S_r}{\partial a_m} = 2 \sum_{i=1}^{n} (y_i - a_0 - a_1x_i - \dots - a_mx_i^m)(-x_i^m)$$

Or:

$$\begin{bmatrix}
n & \sum x_i & \dots & \sum x_i^m \\
\sum x_i & \sum x_i^2 & \dots & \sum x_i^{m+1} \\
\dots & \dots & \dots & \dots \\
\sum x_i^m & \sum x_i^{m+1} & \dots & \sum x_i^{2m}
\end{bmatrix} \cdot
\begin{bmatrix}
a_0 \\
a_1 \\
\dots \\
a_m
\end{bmatrix} =
\begin{bmatrix} \sum y_i \\ \sum x_iy_i \\ \dots \\ \sum x_i^m y_i
\end{bmatrix}$$

The standard error of the estimate is now, $s_{y/x} = \sqrt{\frac{S_r}{n-(m+1)}}$.

## Multiple linear regression
When there are two independent variables, the multi-linear model for regression is: $y = a_0 + a_1x_1 + a_2x_2$. The sum of squares of the residuals is:

$$S_r = \sum_{i=1}^{n} (y_i - a_0 - a_1x_{1i} - a_2x_{2i})^2$$

We differentiate:

$$\frac{\partial S_r}{\partial a_0} = -2 \sum_{i=1}^{n} (y_i - a_0 - a_1x_{1i} - a_2x_{2i})$$
$$\frac{\partial S_r}{\partial a_1} = -2 \sum_{i=1}^{n} (y_i - a_0 - a_1x_{1i} - a_2x_{2i})x_{1i}$$
$$\frac{\partial S_r}{\partial a_2} = -2 \sum_{i=1}^{n} (y_i - a_0 - a_1x_{1i} - a_2x_{2i})x_{2i}$$

Setting the derivatives to zero, the matrix-vector form is:
$$\begin{bmatrix}
n & \sum x_{1i} & \sum x_{2i} \\
\sum x_{1i} & \sum x_{1i}^2 & \sum x_{1i}x_{2i} \\
\sum x_{2i} & \sum x_{1i}x_{2i} & \sum x_{2i}^2
\end{bmatrix} \cdot
\begin{bmatrix}
a_0 \\ a_1 \\ a_2
\end{bmatrix} =
\begin{bmatrix}
\sum y_i \\ \sum x_{1i}y_i \\ \sum x_{2i}y_i
\end{bmatrix}$$

For $m = ndim$ independent variables: $y = a_0 + a_1x_1 + a_2x_2 + \dots + a_m x^m$. The sum of squares of the residuals is:

$$S_r = \sum_{i=1}^n (y_i - a_0 - a_1x_{1i} - a_2x_{2i} + \dots + a_mx_{mi})^2$$

We differentiate:

$$\frac{\partial S_r}{\partial a_0} = -2 \sum_{i=1}^n y_i - a_0 - a_1x_{1i} - a_2x_{2i} - \dots - a_mx_{mi}$$
$$\frac{\partial S_r}{\partial a_1} = -2 \sum_{i=1}^n (y_i - a_0 - a_1x_{1i} - a_2x_{2i} \dots - a_mx_{mi})x_{1i}$$
$$\dots$$
$$\frac{\partial S_r}{\partial a_m} = -2 \sum_{i=1}^n (y_i - a_0 - a_1x_{1i} - a_2x_{2i} \dots - a_mx_{mi})x_{mi}$$

Setting the derivatives to zero, the matrix-vector form is:

$$\begin{bmatrix}
n & \sum x_{1i} & \dots & \sum x_{mi} \\
\sum x_{1i} & \sum x_{1i}^2 & \dots & \sum x_{1i}x_{mi} \\
\dots & \dots & \dots & \dots \\
\sum x_{mi} & \sum x_{1i}x_{mi} & \dots & \sum x_{mi}^2
\end{bmatrix} \cdot
\begin{bmatrix}
a_0 \\ a_1 \\ \dots \\ a_m
\end{bmatrix} =
\begin{bmatrix}
\sum y_i \\ \sum x_{1i}y_i \\ \dots \\ \sum x_{mi}y_i
\end{bmatrix}$$

## General linear regression
Linear regression refers to linearity in the coefficients $a_i$.
The model can still be non-linear in terms of the independent variables.
Multi-linear regression can be used for linear regression of polynomial or generally non-linear models.
We need to convert the list of $x_i$ data to multiple columns of data, one for each term of the model.
For example, when working with a second order polynomial: $y = a_0 + a_1x + a_2x^2$ we can treat it as $y = a_0 + a_1x_1 + a_2x_2$, where $x_1 = x$ and $x_2 = x^2$.
Thus, we created two columns of independent variable data: one for $x$ as it was and one for $x^2$.
We then proceed with multiple linear regression.
In general, for a function
:
$$f(z) = a_0 + a_1z_1 + \dots + a_mz_m$$

the sum of squares is (for $n$ data and $m+1$ coefficients):

$$S_r = \sum_{i=1}^{n} \left( y_i - \sum_{j=0}^{m} a_j z_{ji} \right)^2$$

This is minimized by making the derivatives $\frac{\partial S_r}{\partial a_i} = 0$ which results in:

$$(\underline{\underline{Z}}^T \cdot \underline{\underline{Z}}) \cdot \underline{a} = \underline{\underline{Z}}^T \cdot \underline{y}$$

where

$$\underline{\underline{Z}} = \begin{bmatrix} z_{01} & z_{11} & \dots & z_{m1} \\
z_{02} & z_{12} & \dots & z_{m2} \\
\dots & \dots & \dots & \dots \\ z_{0n} & z_{1n} & \dots & z_{mn} \end{bmatrix}$$
$$\underline{y}^T = [y_1 \quad \dots \quad y_n]$$
$$\underline{a}^T = [a_1 \quad \dots \quad a_m]$$

## Non-linear Regression
When the model is non-linear in terms of its parameters we cannot perform the generalized multi-linear regression.
Instead, we need to solve iterative for the optimal parameters following the Gauss-Newton approach.
We linearize the non-linear expression around the previous value:

$$f(x_i^{(k)}) = f(x_i^{(k-1)}) + \frac{\partial f(x_i^{(k-1)})}{\partial a_0} \delta a_0 + \dots + \frac{\partial f(x_i^{(k-1)})}{\partial a_m} \delta a_m$$

where $\delta a_i = a_i^{(k)} - a_i^{(k-1)}$. The sum of squares is:

$$S_r = \sum_{i=1}^{n} \left( y_i - f(x_i^{(k-1)}) - \sum_{j=0}^{m} \frac{\partial f(x_i^{(k-1)})}{\partial a_j} \delta a_j \right)^2$$

This is minimized by making the derivatives $\frac{\partial S_r}{\partial a_i} = 0$ which results in:

$$(\underline{\underline{J}}^T \cdot \underline{\underline{J}}) \cdot \underline{\delta a}^{(k)} = \underline{\underline{J}}^T \cdot \underline{b}$$

where
$$\underline{\underline{J}} =
\begin{bmatrix}
\frac{\partial f(x_1^{(k-1)})}{\partial a_0} & \frac{\partial f(x_1^{(k-1)})}{\partial a_1} & \dots & \frac{\partial f(x_1^{(k-1)})}{\partial a_m} \\
\frac{\partial f(x_2^{(k-1)})}{\partial a_0} & \frac{\partial f(x_2^{(k-1)})}{\partial a_1} & \dots & \frac{\partial f(x_2^{(k-1)})}{\partial a_m} \\
\dots & \dots & \dots & \dots \\ \frac{\partial f(x_n^{(k-1)})}{\partial a_0} & \frac{\partial f(x_n^{(k-1)})}{\partial a_1} & \dots & \frac{\partial f(x_n^{(k-1)})}{\partial a_m}
\end{bmatrix}$$
$$\underline{b}^T = [y_1 - f(x_1) \quad \dots \quad y_n - f(x_n)]$$
$$\underline{\delta a}^{(k)^T} = [\delta a_1^{(k)} \quad \dots \quad \delta a_m^{(k)}]$$

We approximate the Jacobian matrix with finite differences and solve for the correction $\underline{\delta a}^{(k)}$. The Gauss-Newton method converges fast close to the solution but might diverge away from it.
Optimization algorithms can also be used (for example, steepest descent which can be slow around the minimum).
The standard alternative is the Levenberg-Marquardt algorithm, which combines steepest descent steps away from the solution and Gauss-Newton precision close to the minimum.
It is the standard in many optimization libraries.
It uses a damping factor for the diagonal of the matrix:

$$(\underline{\underline{J}}^T \cdot \underline{\underline{J}} + \lambda \underline{\underline{I}}) \cdot \underline{\delta a}^{(k)} = \underline{\underline{J}}^T \cdot \underline{b}$$

When $\lambda \to 0$ we obtain the Gauss-Newton method, while for $\lambda \to \infty$ it acts like steepest descent.
Unlike linear regression, where there is a global minimum, in non-linear regression we might obtain a local minimum.
