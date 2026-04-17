# Numerical Methods: Linear Systems

We want to solve the $n \times n$ linear system:

$$ \underbar{\underbar{A}} \underbar{x} = \underbar{b} $$

The system is physically interpreted as [Interactions][response]=[stimuli].

## Direct methods
### Cramer's rule

$$ x_j = \frac{detA_j}{detA} $$

where $A_j$ is the matrix for which the column $j$ has been replaced by the right-hand side vector $\underbar{b}$.

It is used for small systems. We need to calculate $n+1$ determinants and the flops are $(n+1)!$.

### Gauss Elimination
The method involves two steps:\
**Forward elimination**:\
$x_i$ is eliminated from equations $i+1$ through $n$. An upper triangular matrix is formed.\
**Back substitution**:\
Solving for $x_i$.\
**Partial pivoting**:\
Performed to avoid division by zero and minimize round-off errors.
At each row, $k$, we compare the diagonal element $A_{kk}$ with $A_{ik} (i > k)$, i.e. the elements below at the same column, $k$.
We swap row $k$ with the row that has the largest absolute value.\
**Complexity**:\
The flops are $O(n^3)$: Elimination $O(\frac{2}{3} n^3)$, Back substitution $O(2n^2)$.

#### Gauss-Jordan
A variation of Gauss elimination where unknowns are eliminated from all equations rather than subsequent ones only.
Also, all rows are normalized by diving by their pivot element.
Thus, instead of an upper triangular matrix, the identity matrix is formed after forward elimination.
The method involves more operations compared to Gauss elimination and it is, thus, not preferred.

### LU Decomposition
The matrix is decomposed into a lower, $\underbar{\underbar{L}}$, and an upper, $\underbar{\underbar{U}}$, triangular matrix such that: $\underbar{\underbar{L}} \underbar{\underbar{U}} = \underbar{\underbar{A}}$.
The right-hand side vector is given by:

$$ \underbar{b} = \underbar{\underbar{L}} \underbar{d} $$

where $\underbar{d}$ is an intermediate vector. The solution is:

$$ \underbar{\underbar{U}} \underbar{x} = \underbar{d} $$

The process of LU decomposition is performed during Gauss Elimination. The number of flops is the same.
The Doolittle variation results in $\underbar{\underbar{U}}$ containing the diagonal elements and $\underbar{\underbar{L}}$ having $1$ on the diagonal.
The Crout decomposition is the opposite.\
Once the LU decomposition of the matrix has been performed, it can be used to find solutions with different right-hand side vectors $\underbar{b}$.
This can be practical in various instances.
For example, the inverse $\underbar{\underbar{A}}^{-1}$ can be calculated as:

$$ \underbar{\underbar{A}} \underbar{x}_j = \underbar{b}_j $$

$\underbar{b}_j$ is the zero vector with $1$ in the $j^{th}$ entry.
The solution $\underbar{x}_j$ is the $j^{th}$ column of $\underbar{\underbar{A}}^{-1}$.

#### Cholesky Decomposition
A special case of LU decomposition for symmetric and positive-definite matrices.
The decomposition results in:

$$ \underbar{\underbar{L}} \underbar{\underbar{L}}^T  = \underbar{\underbar{A}} $$

The lower diagonal $\underbar{\underbar{L}}$ elements are:

$$ L_{ki} = \frac{A_{ki}-\sum_{j=1}^{k-1}L_{ij}L_{kj}}{L_{ii}}, i = 1, ..., k-1 $$

And the diagonal ones are:
$$ L_{kk} = \sqrt{A_{kk}-\sum_{j=1}^{k-1}L_{kj}^2} $$

The number of flops is reduced to half, $O(\frac{1}{3}n^3)$.

### Thomas Algorithm
It is used for tri-diagonal matrices. The total number of flops for decomposition and solution is $O(8n)$.

## Iterative Methods
We start with an initial guess for the solution, $\underbar{x}$ and improve until convergence.
They are utilized for large sparse systems where LU decomposition is costly.

### Jacobi
The solution is updated as:

$$x_i^{(k)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j=1, j \neq i}^n a_{ij} x_j^{(k-1)} \right)$$

$$x_i^{(k)} = x_i^{(k-1)} + \frac{1}{a_{ii}} r_i^{(k-1)}$$

where the residual of row $i$ is:

$$ r_i = b_i - \sum_{j=1}^n a_{ij} x_j^{(k-1)} =
\left( b_i - \sum_{j=1, j \neq i}^n a_{ij} x_j^{(k-1)} \right) + a_{ii} x_i^{(k-1)} $$


In the Jacobi method, newly calculated values of $x_i^{(k)}$ are not used immediately, but only at the next iteration.
The method converges for diagonally dominant matrices:

$$ \left| {a_{ii}} \right| > \sum_{j=1}^n \left| {a_{ij}} \right|$$

### Gauss-Seidel
The solution is updated as:

$$x_i^{(k)} = \frac{1}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k)} - \sum_{j=i+1}^{n} a_{ij} x_j^{(k-1)} \right)$$

$$x_i^{(k)} = x_i^{(k-1)} + \frac{1}{a_{ii}} r_i$$

where the residual of row $i$ is:

$$ r_i = b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k)} - \sum_{j=i}^{n} a_{ij} x_j^{(k-1)} $$
$$ r_i = b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k)} - \sum_{j=i+1}^{n} a_{ij} x_j^{(k-1)} + a_{ii} x_i^{(k-1)} $$

In the Gauss-Seidel method, newly calculated values of $x_i^{(k)}$ are used immediately.
The method converges for diagonally dominant matrices.

### Successive Over-Relaxation (SOR)
The solution is updated as:

$$x_i^{(k)} = \omega x_i^{(k)} + (1-\omega)x_i^{(k-1)}$$

$$x_i^{(k)} = \frac{\omega}{a_{ii}} \left( b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k)} - \sum_{j=i+1}^{n} a_{ij} x_j^{(k-1)} \right) + (1-\omega)x_i^{(k-1)}$$

$$x_i^{(k)} = x_i^{(k-1)} + \frac{\omega}{a_{ii}} r_i$$

where the residual of row $i$ is:

$$ r_i = b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k)} - \sum_{j=i}^{n} a_{ij} x_j^{(k-1)} $$
$$ r_i = b_i - \sum_{j=1}^{i-1} a_{ij} x_j^{(k)} - \sum_{j=i+1}^{n} a_{ij} x_j^{(k-1)} + a_{ii} x_i^{(k-1)} $$

For $\omega=1$ we obtain the Gauss-Seidel method.
To dampen oscillations we set $0 < \omega < 1$.
For faster convergence we use $\omega > 1$, provided that the solution moves in the correct direction.

### Steepest Descent
The optimization algorithm of steepest descent can be used to solve a linear system.
The method works for symmetric and positive-definite matrix.
The cost function:

$$f(\underline{x}) = \frac{1}{2}\underline{x}^T \cdot \underline{\underline{A}} \cdot \underline{x} - \underline{x}^T \cdot \underline{b}$$

has gradient:

$$\nabla f = \underline{\underline{A}} \cdot \underline{x} - \underline{b} = -\underline{r}$$

Thus, the linear system is solved when the function is minimized. The optimal step size can be calculated by setting:

$$\frac{\partial f(\underline{x}^{(k)})}{\partial a} = \frac{\partial f(\underline{x}^{(k-1)} + a^{(k)}\underline{r}^{(k-1)})}{\partial a} = 0$$

$$a^{(k)} = \frac{\underline{r}^{(k-1)T} \cdot \underline{r}^{(k-1)}}{\underline{r}^{(k-1)T} \cdot \underline{\underline{A}} \cdot \underline{r}^{(k-1)}}$$

The new solution is:

$$\underline{x}^{(k)} = \underline{x}^{(k-1)} + a^{(k)}\underline{r}^{(k-1)}$$

and the new residual:

$$\underline{r}^{(k)} = \underline{b} - \underline{\underline{A}} \cdot \underline{x}^{(k)} = \underline{r}^{(k-1)} - a^{(k)}\underline{\underline{A}} \cdot \underline{r}^{(k-1)}$$

### Conjugate Gradient
The method also works for symmetric, positive-definite matrices.
In steepest descent, the new search direction is orthogonal to the previous one and it does not consider all previous ones.
With conjugate gradient, the new search direction is orthogonal to all previous directions:

$$\underline{p}_i^T \cdot \underline{\underline{A}} \cdot \underline{p}_i = 0, \forall i \neq j$$

Instead of using the residual $\underline{r}$ as the new direction, we use a vector, $\underline{p}$, which is a combination of the residual and the previous direction.
Starting with the first search direction $\underline{p}_0 = \underline{r}_0$, the optimal step size is calculated with vector $\underline{p}$:

$$\frac{\partial f(\underline{x}^{(k)})}{\partial a} = \frac{\partial f(\underline{x}^{(k-1)} + a^{(k)}\underline{p}^{(k-1)})}{\partial a} = 0$$

$$a^{(k)} = \frac{\underline{r}^{(k-1)T} \cdot \underline{r}^{(k-1)}}{\underline{p}^{(k-1)T} \cdot \underline{\underline{A}} \cdot \underline{p}^{(k-1)}}$$

The new solution is now updated with vector $\underline{p}$ towards the new best direction (instead of $\underline{r}$):

$$\underline{x}^{(k)} = \underline{x}^{(k-1)} + a^{(k)}\underline{p}^{(k-1)}$$

and the new residual is:

$$\underline{r}^{(k)} = \underline{b} - \underline{\underline{A}} \cdot \underline{x}^{(k)} = \underline{r}^{(k-1)} - a^{(k)}\underline{\underline{A}} \cdot \underline{p}^{(k-1)}$$

The direction vector is updated:
$$\underline{p}^{(k)} = \underline{r}^{(k)} + \beta\underline{p}^{(k-1)}$$

$\beta$ is given by enforcing conjugacy:

$$\underline{p}_i^{T(k)} \cdot \underline{\underline{A}} \cdot \underline{p}_i^{(k-1)} = \left( \underline{r}^{(k)} + \beta\underline{p}^{(k-1)} \right) \cdot \underline{\underline{A}} \cdot \underline{p}_i^{(k-1)} = 0$$

$$\beta = \frac{\underline{r}^{(k)T} \cdot \underline{r}^{(k)}}{\underline{r}^{(k-1)T} \cdot \underline{r}^{(k-1)}}$$

