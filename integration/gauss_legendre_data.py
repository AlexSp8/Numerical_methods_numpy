import numpy as np

sqrt3 = np.sqrt(3.0)
sqrt5 = np.sqrt(5.0)
sqrt06 = np.sqrt(0.6)
sqrt7 = np.sqrt(7.0)
sqrt10 = np.sqrt(10.0)
sqrt30 = np.sqrt(30.0)
sqrt6_5 = np.sqrt(6.0/5.0)
sqrt10_7 = np.sqrt(10.0/7.0)

GAUSS_POINTS = {
    1: [0.0],
    2: [-1.0/sqrt3, 1.0/sqrt3],
    3: [-sqrt06, 0.0, sqrt06],
    4: [
        -np.sqrt( (3.0 + 2.0*sqrt6_5)/7.0 ),
        -np.sqrt( (3.0 - 2.0*sqrt6_5)/7.0 ),
         np.sqrt( (3.0 - 2.0*sqrt6_5)/7.0 ),
         np.sqrt( (3.0 + 2.0*sqrt6_5)/7.0 )
    ],
    5: [
        -(1.0/3.0)*np.sqrt(5.0 + 2.0*sqrt10_7),
        -(1.0/3.0)*np.sqrt(5.0 - 2.0*sqrt10_7),
         0.0,
         (1.0/3.0)*np.sqrt(5.0 - 2.0*sqrt10_7),
         (1.0/3.0)*np.sqrt(5.0 + 2.0*sqrt10_7)
    ]
}

GAUSS_WEIGHTS = {
    1: [2.0],
    2: [1.0, 1.0],
    3: [5.0/9.0, 8.0/9.0, 5.0/9.0],
    4: [
        (18.0 - sqrt30)/36.0,
        (18.0 + sqrt30)/36.0,
        (18.0 + sqrt30)/36.0,
        (18.0 - sqrt30)/36.0
    ],
    5: [
        (322.0 - 13.0*sqrt7*sqrt10)/900.0,
        (322.0 + 13.0*sqrt7*sqrt10)/900.0,
        128.0/225.0,
        (322.0 + 13.0*sqrt7*sqrt10)/900.0,
        (322.0 - 13.0*sqrt7*sqrt10)/900.0
    ]
}
