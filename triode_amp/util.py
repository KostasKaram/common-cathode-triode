
import numpy
import math
from scipy import signal


def gen_sine(A, f, fs, t):
    'Generate a sine signal'
    dt = numpy.arange(0, t, 1/fs)
    sine = A*numpy.sin(2.0*numpy.pi*dt*f)
    return sine, dt


def dc_block(x, fs, fc=30.0, order=4):
    'Apply DC-blocking high-pass filter'
    fcn = fc / (0.5*fs)
    b, a = signal.butter(order, fcn, btype='highpass', analog=False)
    return signal.lfilter(b, a, x)


def soft_clip(x, alpha=2.0):
    'Gentle soft-clip function'
    return x / numpy.sqrt(alpha + numpy.square(x))




