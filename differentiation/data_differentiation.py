
from typing import List

import numpy as np
import numpy.typing as npt

def divided_difference(x: npt.NDArray, y: npt.NDArray) -> List[float]:
    """Returns the values of the derivatives of a dataset, x-y,
    using the divided differences."""

    npts = len(x)

    if npts != len(y):
        raise ValueError("x and y must have the same length for divided difference")
    if npts < 2:
        raise ValueError("Not enough data points for divided difference")

    df = np.zeros(npts-1)
    for i in range(npts-1):
        hi = x[i+1] - x[i]
        df[i] = (y[i+1]-y[i])/hi

    return df
