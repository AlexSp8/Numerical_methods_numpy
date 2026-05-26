
import numpy as np

g2 = (3.0**0.5)/6.0
sqrt6 = 6.0**0.5
sqrt15 = 15.0**0.5
IMPLICIT_METHODS = {
    'euler': {
        'ns': 1,
        'order': 1,
        'p': np.array([1.0]),
        'q': np.array([[1.0]]),
        'b': np.array([1.0])
    },
    'midpoint': {
        'ns': 1,
        'order': 2,
        'p': np.array([0.5]),
        'q': np.array([[0.5]]),
        'b': np.array([1.0])
    },
    'crank-nicolson': {
        'ns': 2,
        'order': 2,
        'p': np.array([0.0, 1.0]),
        'q': np.array([[0.0, 0.0], [0.5, 0.5]]),
        'b': np.array([0.5, 0.5])
    },
    'gauss-legendre-2': {
        'ns': 2,
        'order': 4,
        'p': np.array([0.5-g2, 0.5+g2]),
        'q': np.array([[0.25, 0.25-g2], [0.25+g2, 0.25]]),
        'b': np.array([0.5, 0.5])
    },
    'gauss-legendre-3': {
        'ns': 3,
        'order': 6,
        'p': np.array([0.5-sqrt15/10.0, 0.5, 0.5+sqrt15/10.0]),
        'q': np.array([[5.0/36.0, 2.0/9.0-sqrt15/15.0, 5.0/36.0-sqrt15/30.0],
              [5.0/36.0+sqrt15/24.0, 2.0/9.0, 5.0/36.0-sqrt15/24.0],
              [5.0/36.0+sqrt15/30.0, 2.0/9.0+sqrt15/15.0, 5.0/36.0]]),
        'b': np.array([5.0/18.0, 8.0/18.0, 5.0/18.0])
    },
    'radau-iia-2': {
        'ns': 2,
        'order': 3,
        'p': np.array([1.0/3.0, 1.0]),
        'q': np.array([[5.0/12.0, -1.0/12.0], [0.75, 0.25]]),
        'b': np.array([0.75, 0.25])
    },
    'radau-iia-3': {
        'ns': 3,
        'order': 5,
        'p': np.array([(4-sqrt6)/10, (4+sqrt6)/10, 1.0]),
        'q': np.array([[(88-7*sqrt6)/360, (296-169*sqrt6)/1800, (-2+3*sqrt6)/225],
              [(296+169*sqrt6)/1800, (88+7*sqrt6)/360, (-2-3*sqrt6)/225],
              [(16-sqrt6)/36, (16+sqrt6)/36, 1/9]]),
        'b': np.array([(16+sqrt6)/36, (16-sqrt6)/36, 1/9])
    },
    'lobatto-iiic-2': {
        'ns': 2,
        'order': 2,
        'p': np.array([0.0, 1.0]),
        'q': np.array([[0.5, -0.5],
              [0.5, 0.5]]),
        'b': np.array([0.5, 0.5]),
    },
    'lobatto-iiia-3': {
        'ns': 3,
        'order': 4,
        'p': np.array([0.0, 0.5, 1.0]),
        'q': np.array([[0.0, 0.0, 0.0],
              [5.0/24.0, 1.0/3.0, -1.0/24.0],
              [1.0/6.0, 2.0/3.0, 1.0/6.0]]),
        'b': np.array([1.0/6.0, 2.0/3.0, 1.0/6.0]),
    },
}
