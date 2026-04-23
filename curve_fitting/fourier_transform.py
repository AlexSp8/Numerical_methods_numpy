
from typing import Tuple

import numpy as np
import numpy.typing as npt

def discrete_fourier_transform(f: npt.NDArray) -> Tuple[npt.NDArray]:
    """Returns the coefficients of discrete Fourier Transform
    of a set of N discrete data, f."""

    N = f.shape[0]

    w0 = 2*np.pi/N

    n = np.arange(N)
    k = n.reshape((N, 1))

    M = np.exp(-1j*w0*k*n)

    X = np.dot(M,f)

    re, im = X.real, X.imag

    # mag = (re**2 + im**2)
    # print(mag)

    return re, im

def fft_sande_tukey(f: npt.NDArray) -> Tuple[npt.NDArray]:
    """Returns the coefficients of Sande-Tukey FFT
    of a set of N discrete data, f."""

    N = len(f)
    if (N & (N - 1)) != 0:
        raise ValueError("N must be a power of 2")

    data = np.array(f, dtype=complex)

    n_stages = int(np.log2(N))

    for stage in range(n_stages):

        step = N // (2**stage) # = N
        half_step = step // 2 # = N/2

        twiddles = np.exp(-2j*np.pi*np.arange(half_step)/step)
        for i in range(0, N, step):

            top = data[i : i + half_step]
            bottom = data[i + half_step : i + step]

            temp_top = top + bottom
            data[i + half_step : i + step] = (top - bottom)*twiddles
            data[i : i + half_step] = temp_top

    return bit_reversal(data)

def bit_reversal(data: npt.NDArray) -> Tuple[npt.NDArray]:

    re, im = data.real, data.imag
    N = len(re)
    num_bits = int(np.log2(N))

    for i in range(N):
        reversed_j = 0
        temp_i = i
        for _ in range(num_bits):
            last_digit = temp_i%2
            reversed_j = (reversed_j*2) + last_digit
            temp_i = temp_i//2

        if i < reversed_j:
            re[i], re[reversed_j] = re[reversed_j], re[i]
            im[i], im[reversed_j] = im[reversed_j], im[i]

    return re, im
