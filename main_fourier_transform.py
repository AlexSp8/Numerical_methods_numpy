
from pathlib import Path
from typing import List
import math

from utilities import io_utils, indexing
from curve_fitting import fourier_transform

# Path to the Datasets directory
DATASETS_DIR = Path(__file__).parent / 'curve_fitting/data'


def f_sin():
    N = 16
    fs, dt = 12.5, 0.01
    f = [0.0]*N
    for n in range(N):
        t = n*dt
        f[n] = math.cos(2*math.pi*fs*t)
    return f
    # N = 64
    # f = [0.0]*N
    # for n in range(N):
    #     val = math.sin(2 * math.pi * 3 * n / N)
    #     if val >= 0:
    #         f[n] = 1.0
    #     else:
    #         f[n] = -1.0
    # return f

def main():

    print('\Discrete Fourier Transform')
    f = f_sin()
    re, im = fourier_transform.discrete_fourier_transform(f)
    for k in range(len(re)):
        print(f"k = {k}: ({re[k]:.3f}, {im[k]:.3f})")

    print('\nFFT Sande Tukey')
    re, im = fourier_transform.fft_sande_tukey(f)
    for k in range(len(re)):
        print(f"k = {k}: ({re[k]:.3f}, {im[k]:.3f})")


if __name__ == '__main__':
    main()
