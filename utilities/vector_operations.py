
from typing import List

def norm2(v: List[float]) -> float:
    """Returns the L2 (Euclidean) norm of a vector, v."""

    s = sum(x**2 for x in v)
    return s**0.5

def dot_product(v1: List[float], v2: List[float]) -> float:
    """Returns the dot product of two vectors, v1, v2."""

    if len(v1) != len(v2):
        raise ValueError('Vectors must be of the same length for dot product!')

    return sum( a*b for a, b in zip(v1, v2) )
