
import numpy as np

def nearest_index(xp: float,
    x: np.ndarray[tuple[int], np.dtype[np.float64]]) -> int:
    """Returns the index of the value in sorted x closest to xp.

    Args:
        xp: the point of interest
        x: the sorted array to search in
    Returns:
        The index of the sorted x such that the value of x is closest to xp."""

    idx = np.searchsorted(x, xp)

    if idx == 0:
        return 0

    if idx == len(x):
        return len(x) - 1

    if abs(xp - x[idx-1]) < abs(xp - x[idx]):
        return idx - 1
    else:
        return idx

def starting_index(xp: float,
    x: np.ndarray[tuple[int], np.dtype[np.float64]], m: int = 1) -> int:
    """Returns the starting index of the sub-interval of
    size m+1 of xi that is centered around xp.

    Args:
        xp: the point of interest
        x: the list of points
        m: the number of points to consider for each sub-interval
    Returns:
        The index of the sub-interval that is centered around xp."""

    idx = nearest_index(xp, x)

    i_start = idx - (m//2)

    n = x.shape[0]

    return max(0, min(i_start, n - (m + 1)))

def bubble_sort(x: np.ndarray[tuple[int], np.dtype[np.float64]],
    y: np.ndarray[tuple[int], np.dtype[np.float64]]
    ) -> tuple[np.ndarray[tuple[int], np.dtype[np.float64]], np.ndarray[tuple[int], np.dtype[np.float64]]]:
    """Returns the sorted versions of arrays x, y by performing
    bubble sort with respect to the array x.

    Args:
        x: the array on which the sorting is based
        y: the corresponding independent array
    Returns:
        The sorted versions of arrays x, y."""

    n = len(x)
    x_sorted, y_sorted = x.copy(), y.copy()

    for i in range(n):
        for j in range(n-i-1):
            if x_sorted[j] > x_sorted[j+1]:
                x_sorted[j], x_sorted[j+1] = x_sorted[j+1], x_sorted[j]
                y_sorted[j], y_sorted[j+1] = y_sorted[j+1], y_sorted[j]

    return x_sorted, y_sorted
