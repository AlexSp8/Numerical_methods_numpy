
from typing import Any

import numpy as np
import numpy.typing as npt

def trapezoidal(x: npt.NDArray, y: npt.NDArray) -> float:
    """Returns the value of the integral of a dataset, x-y,
    using the trapezoidal rule."""

    n = len(x)

    if n != len(y):
        raise ValueError("x and y must have the same length for trapezoidal")
    if n < 2:
        raise ValueError("Not enough data points for trapezoidal")

    s = 0.0
    for i in range(n-1):
        hi = x[i+1] - x[i]
        s += hi*(y[i+1]+y[i])/2

    return s

def simpson_13(x: npt.NDArray, y: npt.NDArray) -> float:
    """Returns the value of the integral of a dataset, x-y,
    using the Simpson 1/3 rule."""

    n = len(x)
    ns = n-1

    if n != len(y):
        raise ValueError("x and y must have the same length for Simpson 1/3")
    if ns < 2:
        raise ValueError("Not enough data points for Simpson 1/3")
    if ns%2 != 0:
        raise ValueError("Even number of segments required for Simpson 1/3")

    s = 0.0
    for i in range(0, ns, 2):

        h0 = x[i+1] - x[i]
        h1 = x[i+2] - x[i+1]

        term1 = (2 - h1/h0)*y[i]
        term2 = ((h0 + h1)**2/(h0 * h1))*y[i+1]
        term3 = (2 - h0/h1)*y[i+2]

        s += (h0 + h1)*(term1 + term2 + term3)/6

    return s

def simpson_mixed(x: npt.NDArray, y: npt.NDArray) -> float:
    """Returns the value of the integral of a dataset, x-y, using
    mixed Simpson and trapezoidal rules. For even segments Simpson's 1/3 is used.
    For odd segments trapezoidal is used for the last segment only."""

    n = len(x)
    ns = n-1
    if n != len(y):
        raise ValueError("x and y must have the same length for Simpson 1/3")
    if n < 2:
        raise ValueError("Not enough data points for Simpson mixed")

    if ns == 1:
        return trapezoidal(x, y)

    m = ns
    s = 0.0
    if ns%2 == 1:
        s += trapezoidal(x[-2:], y[-2:])
        m -= 1

    if m > 0:
        s += simpson_13(x[:m+1], y[:m+1])

    return s

def simpson_mixed_2D(x: npt.NDArray, y: npt.NDArray, z: npt.NDArray) -> float:
    """Returns the value of the integral of a dataset in two dimensions
    using the mixed Simpson method."""

    n = len(x)
    s = np.zeros(n)
    for i in range(n):
        s[i] = simpson_mixed(y, z[i])

    return simpson_mixed(x, s)

def simpson_mixed_multi(x: npt.NDArray, y: Any) -> float:
    """Returns the value of the integral of a dataset in
    multiple dimensions using the mixed Simpson method."""

    n = len(x)
    if n == 1:
        return simpson_mixed(x[0], y)

    ndim = len(x[0])
    s = []
    for row in y:
        area = simpson_mixed_multi(x[1:], row)
        s.append(area)

    return simpson_mixed(x[0], s)
