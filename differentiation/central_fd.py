
from typing import Callable

def df_h2(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-8) -> float:
    return (f(x+h) - f(x-h))/(2*h)

def df_h4(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-8) -> float:
    return (-f(x+2*h) + 8*f(x+h) - 8*f(x-h) + f(x-2*h))/(12*h)

def d2f_h2(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-6) -> float:
    if fx is None:
        fx = f(x)
    return (f(x+h) - 2*fx + f(x-h))/(h**2)

def d2f_h4(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-6) -> float:
    if fx is None:
        fx = f(x)
    return (-f(x+2*h) + 16*f(x+h) - 30*fx + 16*f(x-h) - f(x-2*h))/(12*(h**2))

def d3f_h2(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-4) -> float:
    return (f(x+2*h) - 2*f(x+h) + 2*f(x-h) - f(x-2*h))/(2*(h**3))

def d3f_h4(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-4) -> float:
    return (-f(x+3*h) + 8*f(x+2*h) - 13*f(x+h) + 13*f(x-h) - 8*f(x-2*h) + f(x-3*h))/(8*(h**3))

def d4f_h2(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-3) -> float:
    if fx is None:
        fx = f(x)
    return (f(x+2*h) - 4*f(x+h) + 6*fx - 4*f(x-h) + f(x-2*h))/(h**4)

def d4f_h4(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-3) -> float:
    if fx is None:
        fx = f(x)
    return (-f(x+3*h) + 12*f(x+2*h) - 39*f(x+h) + 56*fx - 39*f(x-h) + 12*f(x-2*h) - f(x-3*h))/(6*(h**4))
