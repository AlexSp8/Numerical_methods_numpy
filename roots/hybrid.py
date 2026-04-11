
from typing import Callable

def brent(f: Callable[[float], float], a: float, b: float,
    eps: float = 1e-8, k_max: int = 1000) -> float|None:
    """Returns the root, f(x) = 0, of a function f(x)
    in an interval [a,b] using Brent's method"""

    fa, fb = f(a), f(b)
    if fa*fb >= 0:
        return None

    c = a
    fc = fa
    bisection_flag = True

    for k in range(1, k_max+1):

        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa

        if fa != fc and fb != fc:
            # IQI
            x = (a*fb*fc/( (fa - fb)*(fa - fc) ) +
                 b*fa*fc/( (fb - fa)*(fb - fc) ) +
                 c*fa*fb/( (fc - fa)*(fc - fb) ))
        else:
            # Secant Method
            x = b - fb*(b - a)/(fb - fa)

        check1 = not ((3*a + b)/4 <= x <= b or b <= x <= (3*a + b)/4)
        check2 = bisection_flag and (abs(x - b) >= abs(b - c) / 2)
        check3 = not bisection_flag and (abs(x - b) >= abs(c - d) / 2)
        bisection_flag = check1 or check2 or check3
        if bisection_flag:
            x = (a + b)/2

        fr = f(x)
        d = c # d is used to track "two steps ago"
        c, fc = b, fb

        if fa*fr < 0:
            b, fb = x, fr
        else:
            a, fa = x, fr

        if fb == 0 and abs(b - a) < eps:
            break

    return b
