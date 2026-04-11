
from typing import Callable

def df_h(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-8) -> float:
    if fx is None:
        fx = f(x)
    return (f(x+h) - fx)/h

def df_h2(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-8) -> float:
    if fx is None:
        fx = f(x)
    return (-f(x+2*h) + 4*f(x+h) - 3*fx)/(2*h)

def d2f_h(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-6) -> float:
    if fx is None:
        fx = f(x)
    return ( f(x+2*h) - 2*f(x+h) + fx )/(h**2)

def d2f_h2(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-6) -> float:
    if fx is None:
        fx = f(x)
    return ( -f(x+3*h) + 4*f(x+2*h) - 5*f(x+h) + 2*fx )/(h**2)

def d3f_h(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-4) -> float:
    if fx is None:
        fx = f(x)
    return ( f(x+3*h) - 3*f(x+2*h) + 3*f(x+h) - fx )/(h**3)

def d3f_h2(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-4) -> float:
    if fx is None:
        fx = f(x)
    return ( -3*f(x+4*h) + 14*f(x+3*h) - 24*f(x+2*h) + 18*f(x+h) - 5*fx )/(2*(h**3))

def d4f_h(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-3) -> float:
    if fx is None:
        fx = f(x)
    return ( f(x+4*h) - 4*f(x+3*h) + 6*f(x+2*h) - 4*f(x+h) + fx )/(h**4)

def d4f_h2(f: Callable[[float], float], x: float,
    fx: float = None, h: float = 1e-3) -> float:
    if fx is None:
        fx = f(x)
    return ( -2*f(x+5*h) + 11*f(x+4*h) - 24*f(x+3*h) + 26*f(x+2*h) - 14*f(x+h) + 3*fx )/(h**4)
