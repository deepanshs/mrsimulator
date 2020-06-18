# -*- coding: utf-8 -*-
"""The Event class."""
from typing import List
from typing import Optional

import numpy as np
from mrsimulator.utils.parseable import Parseable
from numpy.fft import fft
from numpy.fft import fftshift
from numpy.fft import ifft
from numpy.fft import ifftshift


class Apodization(Parseable):
    """The Apodization class."""

    function: str
    dimension: int = 0
    args: List[float] = [0]
    fraction: Optional[float] = 1

    class Config:
        validate_assignment = True

    def apodize(self, csdm, **kwargs):
        """Returns the result of passing the selected apodization function .

        Args:
            csdm: simulation object

        Returns:
            A Numpy array
        """
        function_mapping = {"Gaussian": Gaussian, "Lorentzian": Lorentzian}

        y = csdm.dependent_variables[0].components[0]
        # x = csdm.dimensions[self.dimension].coordinates
        axis = -self.dimension - 1
        fapp = function_mapping[self.function](
            csdm, self.dimension, self.args, **kwargs
        )

        recip_dim = csdm.dimensions[0].reciprocal
        coord = csdm.dimensions[0].coordinates
        phase = np.exp(2j * np.pi * recip_dim.coordinates_offset * coord)

        TimeDomain = ifft(ifftshift(y * phase, axes=axis), axis=axis)

        appodized = TimeDomain * fapp
        fft_output = fftshift(fft(appodized, axis=axis), axes=axis).real

        # csdm.dependent_variables[index].components[0] = (
        #     self.fraction * phase.conj() * fft_output
        # )

        return self.fraction * phase.conj() * fft_output


def Lorentzian(csdm, axis, arg, **kwargs):
    """Lorentzian apodization function.

    Args:
        Self: simulation object
        arg: list with one entry. The full-width-half-max in Hz of the Lorentzian function

    Returns:
        A Numpy array
    """
    inv_x = csdm.dimensions[axis].reciprocal_coordinates().to("s").value
    return np.exp(-arg[0] * np.pi * np.abs(inv_x))


def Gaussian(csdm, axis, arg, **kwargs):
    """Gaussian apodization function.

    Args:
        Self: simulation object
        sigma: list with one entry. standard deviation in Hz of the Gaussian function

    Returns:
        A Numpy array
    """
    inv_x = csdm.dimensions[axis].reciprocal_coordinates().to("s").value
    return np.exp(-((inv_x * arg[0] * np.pi) ** 2) * 2)


class PostSimulator(Parseable):
    """scale: scaling factor
    """

    scale: Optional[float] = 1.0
    apodization: List[Apodization] = []
    truncation_artifact: Optional[bool] = False
    dead_time: Optional[float] = 0.0

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True
