from dataclasses import dataclass
from dataclasses import field

import numpy as np
from mrsimulator import Method
from mrsimulator import Simulator
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.simulator import ConfigSimulator
from mrsimulator.spin_system.isotope import Isotope

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


def get_principal_components(zeta, eta):
    """Return the principal components of a traceless second-rank symmetric
    Cartesian tensor.

    Args:
        zeta: The zeta parameter in PAS, according to the Haeberlen convention.
        eta: The eta parameter in PAS, according to the Haeberlen convention.
    """
    xx = -0.5 * zeta * (eta + 1.0)
    yy = 0.5 * zeta * (eta - 1.0)
    zz = zeta

    return [xx, yy, zz]


def get_Haeberlen_components(tensors):
    """Return zeta and eta parameters of the tensor using the Haeberlen convention.

    Args:
        ndarray tensors: A `N x 3 x 3` ndarray of `N` traceless symmetric second-rank
            Cartesian tensors.
    """
    n = tensors.shape[0]
    eig_val = np.linalg.eigvalsh(tensors)
    eig_val_sort_ = np.argsort(np.abs(eig_val), axis=1, kind="quicksort")
    eig_val_sort_ = (eig_val_sort_.T + 3 * np.arange(n)).T.ravel()
    eig_val_sorted = eig_val.ravel()[eig_val_sort_]
    eig_val_sorted.shape = (n, 3)

    zeta = eig_val_sorted[:, -1]
    eta = (eig_val_sorted[:, 0] - eig_val_sorted[:, 1]) / zeta

    return zeta, eta


def zeta_eta_to_x_y(zeta, eta):
    """Convert a set of (zeta, eta) coordinates from the Haeberlen convention to a set
    of (x, y) coordinates.
    """
    xa = np.empty(zeta.size)
    ya = np.empty(zeta.size)

    index = np.where(zeta >= 0)
    temp = np.tan(0.7853981634 * eta[index])
    ya[index] = np.sqrt(zeta[index] * zeta[index] / (temp * temp + 1.0))
    xa[index] = temp * ya[index]

    index = np.where(zeta < 0)
    temp = np.tan(0.7853981634 * eta[index])
    xa[index] = np.sqrt(zeta[index] * zeta[index] / (temp * temp + 1.0))
    ya[index] = temp * xa[index]

    zeta = eta = None
    del zeta, eta

    return xa, ya


def x_y_to_zeta_eta(x, y):
    """Convert a set of (x, y) coordinates defined by two numpy arrays into equivalent
    (zeta, eta) coordinates.
    """
    x = np.abs(x)
    y = np.abs(y)
    zeta = np.sqrt(x**2 + y**2)  # + offset
    eta = np.ones(zeta.shape)
    index = np.where(x > y)
    zeta[index] *= -1
    eta[index] = (4.0 / np.pi) * np.arctan(y[index] / x[index])

    index = np.where(x < y)
    eta[index] = (4.0 / np.pi) * np.arctan(x[index] / y[index])

    return zeta, eta


def _simulate_spectra_over_zeta_and_eta(ZZ, ee, method, tensor_type, config):
    """Helper function to generate the kernel"""
    isotope = Isotope.parse(
        method.channels[0]
    ).symbol  # Grab isotope from Method object

    spin_systems = [
        (
            SpinSystem(sites=[Site(isotope=isotope, quadrupolar=dict(Cq=Z, eta=e))])
            if tensor_type == "quadrupolar"
            else SpinSystem(
                sites=[Site(isotope=isotope, shielding_symmetric=dict(zeta=Z, eta=e))]
            )
        )
        for Z, e in zip(ZZ.ravel(), ee.ravel())
    ]
    sim = Simulator(spin_systems=spin_systems, methods=[method])
    sim.config = config
    sim.config.decompose_spectrum = "spin_system"
    sim.run(pack_as_csdm=False)

    amp = sim.methods[0].simulation.real
    return amp


@dataclass
class LineShapeKernel:
    """lineshape kernel object

    Arguments:
            (tuple) pos: A tuple of numpy arrays defining the 2D grid space.
            (bool) polar: If true, the grid is defined in polar coordinates.
            (mrsimulator.Method) method: The :py:class:`~mrsimulator.method.Method`
                used to simulate the spectra.
            (ConfigSimulator) config: Simulator config to be used in simulation.
    """

    pos: list
    method: Method
    kernel: np.ndarray = None
    polar: bool = False
    config: ConfigSimulator = field(default_factory=ConfigSimulator())

    def generate_lineshape(self, tensor_type: str = "shielding") -> np.ndarray:
        """Pre-compute a lineshape kernel to use for the least-squares fitting of an
        experimental spectrum. The isotope for the spin system is the isotope
        at index zero of the method's channel.

        Note: The lineshape kernel is currently limited to simulating the spectrum of
        an uncoupled spin system for a single Method over a range of nuclear shielding
        or quadrupolar tensors. Functionality may expand in the future.

        Arguments:

            (str) tensor_type: A string enumeration literal describing the type of
                tensor to use in the simulation. The allowed values are `shielding` and
                `quadrupolar`.

        Returns:
            (np.ndarray) kernel: The simulated lineshape kernel
        """
        if tensor_type not in {"shielding", "quadrupolar"}:
            raise ValueError(f"Unrecognized value of {tensor_type} for `tensor_type.")

        # If in polar coordinates, then (ZZ, ee) right now is really (xx, yy)
        ZZ, ee = np.meshgrid(self.pos[0], self.pos[1], indexing="xy")

        # Convert polar coords to cartesian coords
        if self.polar:
            ZZ, ee = x_y_to_zeta_eta(ZZ.ravel(), ee.ravel())

        self.kernel = _simulate_spectra_over_zeta_and_eta(
            ZZ, ee, self.method, tensor_type, self.config
        )
