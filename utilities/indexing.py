
from typing import List, Tuple

def nearest_index(xp: float, x: List[float] = None) -> int:
    """Returns the index of the value in sorted xi closest to xp."""

    idx = 0
    min_diff = abs(xp - x[0])
    for i in range(1, len(x)):
        diff = abs(xp - x[i])
        if diff < min_diff:
            min_diff = diff
            idx = i
    return idx

def starting_index(xp: float, x: List[float] = None, m: int = 1) -> int:
    """Returns the starting index of the sub-interval of
    size m+1 of xi that is centered around xp."""

    n = len(x)
    idx = nearest_index(xp, x)

    i_start = idx - (m//2)

    if i_start < 0:
        i_start = 0
    if i_start + m >= n:
        i_start = n - (m + 1)

    return i_start

def bubble_sort(x: List[float], y: List[float]) -> Tuple[List[float]]:

    n = len(x)
    x_sorted, y_sorted = x[:], y[:]

    for i in range(n):
        for j in range(n-i-1):
            if x_sorted[j] > x_sorted[j+1]:
                temp_x = x_sorted[j]
                x_sorted[j] = x_sorted[j+1]
                x_sorted[j+1] = temp_x

                temp_y = y_sorted[j]
                y_sorted[j] = y_sorted[j+1]
                y_sorted[j+1] = temp_y

    return x_sorted, y_sorted
