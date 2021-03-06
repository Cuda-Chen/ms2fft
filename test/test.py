#!/usr/local/python

import numpy as np
import matplotlib.pyplot as plt

def fft_and_print_result(inputfile, freq):
    s = np.loadtxt(inputfile)
    start = 20000
    end = 40000
    s = s[start:end]
    n = end - start
    t = np.linspace(0, n / freq, n);

    # plot time domain figure
    plt.ylabel("Amplitude")
    plt.xlabel("Time [s]")
    plt.plot(t, s)
    plt.show()

    fft = np.fft.fft(s)
    #print(fft[0])

    for i in range(2):
        print("Value at index {}:\t{}".format(i, fft[i + 1]), "\nValue at index {}:\t{}".format(fft.size -1 - i, fft[-1 - i]))

    T = t[1] - t[0]
    N = n

    # 1 / T = frequency
    f = np.linspace(0, 1 / T, N)

    # plot frequency domain figure
    plt.ylabel("Amplitude")
    plt.xlabel("Frequency [Hz]")
    #plt.bar(f[:N // 2], np.abs(fft)[:N // 2] * 1 / N, width=1.5)  # 1 / N is a normalization factor
    plt.plot(f[:N // 2], np.abs(fft)[:N // 2] * 1 / N)
    plt.show()

def print_result(fftInputFile, freq):
    fft = np.loadtxt(fftInputFile, usecols = (0,))
    for i in range(2):
        print("Value at index {}:\t{}".format(i, fft[i + 1]), "\nValue at index {}:\t{}".format(fft.size -1 - i, fft[-1 - i]))

    N = fft.size - 1
    t = np.linspace(0, N / freq, N)
    f = np.linspace(0, 1 / (t[1] - t[0]), N)

    # plot frequency domain figure
    plt.ylabel("Amplitude")
    plt.xlabel("Frequency [Hz]")
    plt.plot(f[:N // 2], np.abs(fft)[:N // 2] * 1 / N)
    plt.show()

if __name__ == '__main__':
    fft_and_print_result('dumpdata.txt', 100)
    #print_result('fftoutput.txt', 100)
