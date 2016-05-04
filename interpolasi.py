#!/usr/bin/python
import numpy as np
from scipy.io.wavfile import read
import scipy.signal as sig
import sys

def load(path='/home/cilsat/dev/speech/test.wav'):
    awav = read(path)
    return awav

def linear(x, n, o):
    pass

def lagrange(x, n, o):
    pass

def hermite(x, n, o):
    pass

def evaluate(func):
    pass

# wav to oversample. n oversample ratio
def upsample(wav, n):
    # zero-stuffing
    w = wav[1].reshape((wav[1].shape[0], 1))
    npad = (0,n-1)
    w = np.pad(w, npad, 'constant', constant_values=0)[:-n+1]
    ws = np.hstack(w)
    # lowpass filter
    del_f = wav[0]*(n*0.5 - 1.)
    b = sig.firwin(numtaps, 1./n)
    wr = sig.lfilter(b, [1.0], ws)
    return np.array(wr*n, dtype=np.int16)

def resample(wav, target_sr):
    w = wav[1]
    sr = wav[0]

    

if __name__ == "__main__":
    main(sys.argv)
