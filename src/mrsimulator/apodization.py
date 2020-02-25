# -*- coding: utf-8 -*-
import matplotlib
import matplotlib.pyplot as plt
import numpy as np
from numpy.fft import fft
from numpy.fft import fftshift
from numpy.fft import ifft
from numpy.fft import ifftshift


class Apodization:
    def __init__(self, sim, dimension=0):
        self.dimension = dimension
        self.fft = True
        self.x = sim.dimensions[dimension].coordinates
        self.y = sim.dependent_variables[0].components[0]
        self.inverse_x = reciprocal_coordinates(self.x)

    def Lorentzian(self, sigma):
        return np.exp(-sigma * np.pi * np.abs(self.inverse_x))

    def Gaussian(self, sigma):
        return np.exp(-((self.inverse_x * sigma * np.pi) ** 2) * 2)

    #     def Stretch_Exponential(self, beta, dim):
    # return np.exp(-(self.inverse_x**beta))
    # Bartlett
    # Blackman
    # Connes
    # Cosine
    # Hamming
    # Hanning
    # Uniform
    # Welch
    def apodize(self, fn, **kwargs):
        axis = -self.dimension - 1
        fapp = fn(self, **kwargs)
        TimeDomain = ifft(ifftshift(self.y, axes=axis), axis=axis)
        TimeDomain = np.roll(TimeDomain, int(self.x.size / 2), axis=axis)

        appodized = np.roll(TimeDomain * fapp, -int(self.x.size / 2), axis=axis)

        return fftshift(fft(appodized, axis=axis), axes=axis).real


def reciprocal_coordinates(x):
    delta_x = abs(x[1] - x[0]).to("Hz")
    t_indices = np.arange(x.size) - int(x.size / 2)
    t_increment = (1 / (x.size * delta_x)).to("s").value
    return t_indices * t_increment


# def Lorentzian( sigma):
#     return np.exp(-sigma * np.pi * np.abs(inverse_x))

# def Gaussian(sigma):
#     return np.exp(-((inverse_x * sigma * np.pi) ** 2) * 2)
