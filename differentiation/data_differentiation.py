
from typing import List

def divided_difference(x: List[float] = None, y: List[float] = None) -> List[float]:
    """Returns the values of the derivatives of a dataset, x-y,
    using the divided differences."""

    np = len(x)

    if np != len(y):
        raise ValueError("x and y must have the same length for divided difference")
    if np < 2:
        raise ValueError("Not enough data points for divided difference")

    df = [0.0]*(np-1)
    for i in range(np-1):
        hi = x[i+1] - x[i]
        df[i] = (y[i+1]-y[i])/hi

    return df
