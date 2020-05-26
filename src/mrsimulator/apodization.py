# -*- coding: utf-8 -*-
import numpy as np
from numpy.fft import fft
from numpy.fft import fftshift
from numpy.fft import ifft
from numpy.fft import ifftshift


class Apodization_deprecated:
    """
    The Apodization class.

    Attributes:
        dimension: List of dimensions.
        fft: An FFT algorithm may be applied to the data
        x: array of frequency domain coordinates
        y: array of amplitudes
        inverse_x: array of time domain coordinates
    """

    def __init__(self, sim, dimension=0):
        self.dimension = dimension
        self.fft = True
        dim = sim.dimensions[dimension]
        self.x = dim.coordinates
        self.y = sim.dependent_variables[0].components[0]
        self.inverse_x = reciprocal_coordinates(dim).to("s").value
        # self.inverse_x = reciprocal_coordinates(dim.coordinates).to("s").value

    def Lorentzian(self, sigma):
        """Lorentzian apodization function.

        Args:
            Self: simulation object
            sigma: The full-width-half-max in Hz of the Lorentzian function

        Returns:
            A Numpy array
        """
        return np.exp(-sigma * np.pi * np.abs(self.inverse_x))

    def Gaussian(self, sigma):
        """Gaussian apodization function.

        Args:
            Self: simulation object
            sigma: standard deviation in Hz of the Gaussian function

        Returns:
            A Numpy array
        """
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
        """Returns the result of passing the selected apodization function .

        Args:
            Self: simulation object
            fn: The apodization function. See the Apodization member functions for a
                list of available functions

        Returns:
            A Numpy array
        """
        axis = -self.dimension - 1
        fapp = fn(self, **kwargs)
        TimeDomain = ifft(ifftshift(self.y, axes=axis), axis=axis)
        TimeDomain = np.roll(TimeDomain, int(self.x.size / 2), axis=axis)

        appodized = np.roll(TimeDomain * fapp, -int(self.x.size / 2), axis=axis)

        return fftshift(fft(appodized, axis=axis), axes=axis).real


def reciprocal_coordinates(self):
    """Return reciprocal coordinates assuming Nyquist-shannan theorem."""
    count = self.count
    increment = 1.0 / (count * self.increment)
    coordinates_offset = self.reciprocal.coordinates_offset
    coordinates = np.arange(count) * increment + coordinates_offset
    if self.complex_fft:
        return coordinates
    return coordinates - int(count / 2) * increment


# def reciprocal_coordinates(x):
#     """Returns the result of passing the selected apodization function .

#     Args:
#         x: array of frequency domain coordinates

#     Returns:
#         A Numpy array
#     """
#     delta_x = abs(x[1] - x[0]).to("Hz", 'nmr_frequency_ratio')
#     t_indices = np.arange(x.size) - int(x.size / 2)
#     t_increment = (1 / (x.size * delta_x)).to("s").value
#     return t_indices * t_increment


# def Lorentzian( sigma):
#     return np.exp(-sigma * np.pi * np.abs(inverse_x))

# def Gaussian(sigma):
#     return np.exp(-((inverse_x * sigma * np.pi) ** 2) * 2)
