
from typing import List, Tuple
import math

def discrete_fourier_transform(f: List[float]) -> Tuple[List[float]]:
    """Returns the coefficients of discrete Fourier Transform
    of a set of N discrete data, f."""

    N = len(f)
    w0 = 2*math.pi/N

    re, im = [0.0]*N, [0.0]*N
    for k in range(N):
        for n in range(N):
            angle = k*w0*n
            re[k] += f[n]*math.cos(angle)
            im[k] -= f[n]*math.sin(angle)

    mag = [ (r**2 + i**2) for r, i in zip(re,im) ]
    # print(mag)

    return re, im

def fft_sande_tukey(f: List[float]) -> Tuple[List[float]]:
    """Returns the coefficients of Sande-Tukey FFT
    of a set of N discrete data, f."""

    N = len(f)
    if (N & (N - 1)) != 0:
        raise ValueError("N must be a power of 2")

    # f is the input signal; f_imag starts as zeros
    re = [float(x) for x in f]
    im = [0.0]*N

    n_stages = int(math.log2(N))

    # 1. Butterfly Stages
    for stage in range(n_stages):

        step = N // (2**stage) # = N
        half_step = step // 2 # = N/2

        angle_step = (2.0*math.pi)/step

        for i in range(0, N, step):
            for j in range(i, i + half_step):

                theta = angle_step*(j - i)

                cos_t = math.cos(theta)
                sin_t = math.sin(theta)

                # Difference components
                diff_re = re[j] - re[j + half_step]
                diff_im = im[j] - im[j + half_step]

                # Top part:
                re[j] = re[j] + re[j + half_step]
                im[j] = im[j] + im[j + half_step]

                # Bottom part: (diff_re + j*diff_im) * (cos_t - j*sin_t)
                re[j + half_step] = diff_re*cos_t + diff_im*sin_t
                im[j + half_step] = diff_im*cos_t - diff_re*sin_t

    return bit_reversal(re, im)

def bit_reversal(re: List[float], im: List[float]) -> Tuple[List[float]]:

    N = len(re)
    # How many bits do we need? (e.g., for N=8, we need 3 bits because 2^3 = 8)
    num_bits = int(math.log2(N))

    for i in range(N):

        j = 0
        temp_i = i

        for step in range(num_bits):
            # 1. Get the last 'digit' of the binary version of temp_i
            # (In math, i % 2 tells us if a number is even (0) or odd (1))
            last_digit = temp_i % 2

            # 2. Shift our new 'j' to make room, then add the digit
            # This is like building a number from right to left
            j = (j * 2) + last_digit

            # 3. Move to the next digit of temp_i
            temp_i = temp_i // 2

        # Swap the elements, but only once!
        if i < j:
            re[i], re[j] = re[j], re[i]
            im[i], im[j] = im[j], im[i]

    return re, im
