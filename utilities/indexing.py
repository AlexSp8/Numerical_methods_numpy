
from typing import Tuple

import numpy as np

def nearest_index(xp: float, x: np.NDArray,) -> int:
    """Returns the index of the value in sorted xi closest to xp."""

    idx = np.searchsorted(x, xp)

    if idx == 0:
        return 0

    if idx == len(x):
        return len(x) - 1

    if abs(xp - x[idx-1]) < abs(xp - x[idx]):
        return idx - 1
    else:
        return idx

def starting_index(xp: float, x: np.NDArray, m: int = 1) -> int:
    """Returns the starting index of the sub-interval of
    size m+1 of xi that is centered around xp."""

    n = x.shape[0]
    idx = nearest_index(xp, x)

    i_start = idx - (m//2)

    return max(0, min(i_start, n - (m + 1)))

def bubble_sort(x: np.NDArray, y: np.NDArray) -> Tuple[np.NDArray]:

    n = len(x)
    x_sorted, y_sorted = x.copy(), y.copy()

    for i in range(n):
        for j in range(n-i-1):
            if x_sorted[j] > x_sorted[j+1]:
                x_sorted[j], x_sorted[j+1] = x_sorted[j+1], x_sorted[j]
                y_sorted[j], y_sorted[j+1] = y_sorted[j+1], y_sorted[j]

    return x_sorted, y_sorted
