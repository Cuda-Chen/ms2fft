# ms2fft
Read miniSEED data then make FFT on it.

# Dependencies
- libmseed
- fftw

# How to Compile and Run
First clone this repo with this command:
```
$ git clone --recursive https://github.com/Cuda-Chen/ms2fft.git
```

Then type:
```
$ make
$ ./ms2fft
```

# What is the meaning of out[0]?
- http://www.fftw.org/fftw3_doc/Complex-One_002dDimensional-DFTs.html#Complex-One_002dDimensional-DFTs
    - The DFT results are stored in-order in the array out, with the zero-frequency (DC) component in out[0]. 
- Zero-frequency component
    - https://www.quora.com/What-is-a-zero-Hz-frequency-component
    - Short answer: the zero Hz component is the average of the signal in the time-domain. Also called its DC value.
