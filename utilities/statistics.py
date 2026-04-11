
from typing import List

def mean(xi: List[float] = None) -> float:
    """Returns the mean value of a set of data, xi."""

    n = len(xi)
    sum_x = 0.0
    for i in range(n):
        sum_x += xi[i]

    return sum_x/n

def sum_of_squares_around_mean(xi: List[float] = None) -> float:
    """Returns the sum of squares of a set of data, xi."""

    n = len(xi)

    x_mean = mean(xi)

    St = 0.0
    for i in range(n):
        St += (xi[i]-x_mean)**2

    return St

def standard_deviation(xi: List[float] = None) -> float:
    """Returns the standard deviation of a set of data, xi."""

    n = len(xi)

    St = sum_of_squares_around_mean(xi)

    return (St/(n-1))**0.5

def variance(xi: List[float] = None) -> float:
    """Returns the variance of a set of data, xi."""

    n = len(xi)
    sum_x, sum_x2 = 0.0, 0.0
    for i in range(n):
        sum_x += xi[i]
        sum_x2 += xi[i]**2

    sy = (sum_x2 - (sum_x**2)/n)/(n-1)
    # sy = standard_deviation(xi)**2
    return sy

def coefficient_of_variation(xi: List[float] = None) -> float:
    """Returns the coefficient of variation of a set of data, xi."""

    return 100*standard_deviation(xi)/mean(xi)
