
from typing import Callable

def df_h(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-8) -> float:
    if fx is None:
        fx = f(x)
    return (fx - f(x-h))/h

def df_h2(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-8) -> float:
    if fx is None:
        fx = f(x)
    return (3*fx - 4*f(x-h) + f(x-2*h))/(2*h)

def d2f_h(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-6) -> float:
    if fx is None:
        fx = f(x)
    return (fx - 2*f(x-h) + f(x-2*h))/(h**2)

def d2f_h2(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-6) -> float:
    if fx is None:
        fx = f(x)
    return (2*fx - 5*f(x-h) + 4*f(x-2*h) - f(x-3*h))/(h**2)

def d3f_h(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-4) -> float:
    if fx is None:
        fx = f(x)
    return (fx - 3*f(x-h) + 3*f(x-2*h) - f(x-3*h))/(h**3)

def d3f_h2(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-4) -> float:
    if fx is None:
        fx = f(x)
    return (5*fx - 18*f(x-h) + 24*f(x-2*h) - 14*f(x-3*h) + 3*f(x-4*h))/(2*(h**3))

def d4f_h(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-3) -> float:
    if fx is None:
        fx = f(x)
    return (fx - 4*f(x-h) + 6*f(x-2*h) - 4*f(x-3*h) + f(x-4*h))/(h**4)

def d4f_h2(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-3) -> float:
    if fx is None:
        fx = f(x)
    return (3*fx - 14*f(x-h) + 26*f(x-2*h) - 24*f(x-3*h) + 11*f(x-4*h) - 2*f(x-5*h))/(h**4)
