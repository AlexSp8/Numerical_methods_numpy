# Numerical Methods: Non-Linear Systems

We want to solve the $n \times n$ non-linear system:

$$ \underbar{F}(\underbar{x}) = \underbar{0} $$

## Newton-Raphson Method
We linearize the equations around the previous solution $\underline{x}^{(k-1)}$ using a Taylor series expansion for each equation:

$$\underline{F}(\underline{x}^{(k)}) = \underline{F}(\underline{x}^{(k-1)}) + dx_1^{(k)} \frac{\partial \underline{F}(\underline{x}^{(k-1)})}{\partial x_1} + \dots + dx_n^{(k)} \frac{\partial \underline{F}(\underline{x}^{(k-1)})}{\partial x_n}$$

The correction of the current iteration is: $dx_i^{(k)} = x_i^{(k)} - x_i^{(k-1)}$.
We set $\underline{F}(\underline{x}^{(k)}) = \underline{0}$:

$$dx_1^{(k)} \frac{\partial \underline{F}(\underline{x}^{(k-1)})}{\partial x_1} + \dots + dx_n^{(k)} \frac{\partial \underline{F}(\underline{x}^{(k-1)})}{\partial x_n} = -\underline{F}(\underline{x}^{(k-1)})$$

In matrix form:

$$\underline{\underline{J}}^{(k-1)} \cdot d\underline{x}^{(k)} = -\underline{F}(\underline{x}^{(k-1)})$$

$\underline{\underline{J}} = \frac{\partial \underline{F}}{\partial \underline{x}}$ is the Jacobian of the system. We solve for the correction $d\underline{x}^{(k)}$.

