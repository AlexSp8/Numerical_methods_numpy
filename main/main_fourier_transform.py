
import sys
from pathlib import Path

import numpy as np

current_dir = Path(__file__).resolve().parent
parent_dir = current_dir.parent

if str(parent_dir) not in sys.path:
    sys.path.insert(0, str(parent_dir))

from curve_fitting import fourier_transform

# Path to the Datasets directory
DATASETS_DIR = parent_dir/'curve_fitting/data'

def f_sin():
    N = 16
    fs, dt = 12.5, 0.01
    f = np.zeros(N)
    for n in range(N):
        t = n*dt
        f[n] = np.cos(2*np.pi*fs*t)
    return f
    # N = 64
    # f = np.zeros(N)
    # for n in range(N):
    #     val = np.sin(2 * np.pi * 3 * n / N)
    #     if val >= 0:
    #         f[n] = 1.0
    #     else:
    #         f[n] = -1.0
    # return f

def main():

    print('\nDiscrete Fourier Transform')
    f = f_sin()
    re, im = fourier_transform.discrete_fourier_transform(f)
    for k in range(re.shape[0]):
        print(f"k = {k}: ({re[k]:.3f}, {im[k]:.3f})")

    print('\nFFT Sande Tukey')
    re, im = fourier_transform.fft_sande_tukey(f)
    for k in range(re.shape[0]):
        print(f"k = {k}: ({re[k]:.3f}, {im[k]:.3f})")


if __name__ == '__main__':
    main()
