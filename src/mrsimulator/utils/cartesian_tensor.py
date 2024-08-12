from typing import List
from typing import Tuple
from typing import Union

import numpy as np
from mrsimulator.spin_system.isotope import Isotope
from mrsimulator.spin_system.tensors import SymmetricTensor
from scipy.spatial.transform import Rotation as R

__author__ = "Philip Grandinetti"
__email__ = "grandinetti.1@osu.edu"


def to_mehring_params(tensor: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
    """Cartesian 3x3 tensor to Mehring parameters. Note, for non-symmetric
    tensors, the conversion is applied by first symmetrizing the tensor
    using (tensor + tensor.T) / 2.

    Args:
        tensor:
            A 3x3 np.ndarray Cartesian tensor.
    Returns:
        A tuple of (euler_angles, eigenvalues) where ``euler_angles`` is an
        ndarray  of three Euler angles [``alpha``, ``beta``, ``gamma``] using
        the  "zyz" convention, and ``eigenvalues`` are the corresponding
        ndarray of three Eigenvalues.
    """
    tensor = (tensor + tensor.T) / 2  # Make sure the tensor is symmetric
    # Calculate the eigenvalues and eigenvectors of the traceless tensor
    eigenvalues, L = np.linalg.eigh(tensor)
    idx = np.argsort(eigenvalues)[::-1]  # Sort the eigenvalues in descending order
    eigenvalues = eigenvalues[idx]  # Swap the eigenvalues
    L = L[:, idx]  # Swap the columns of the eigenvector matrix
    # Flip the sign of the eigenvectors if the determinant is negative
    if np.linalg.det(L) < 0:
        L = -L

    r = R.from_matrix(L)
    euler_angles = r.as_euler("zyz")

    return euler_angles, eigenvalues


def from_mehring_params(
    euler_angles: List[float], eigenvalues: List[float]
) -> np.ndarray:
    """Mehring parameters to a 3x3 symmetric Cartesian tensor.

    Args:
        euler_angles:
            An ndarray of three Euler angles [``alpha``, ``beta``, ``gamma``]
            using the "zyz" convention.
        eigenvalues:
            The corresponding ndarray of three Eigenvalues.

    Returns:
        tensor: A 3x3 np.ndarray of symmetric Cartesian tensor.
    """
    # Assuming alpha_s, beta_s, gamma_s are your Euler angles
    r = R.from_euler("zyz", euler_angles)
    L = r.as_matrix()

    diag_matrix = np.diag(eigenvalues)
    tensor = L @ diag_matrix @ L.T
    return tensor


def to_haeberlen_params(
    tensor: np.ndarray,
) -> Tuple[np.ndarray, float, float, float]:
    """Cartesian 3x3 tensor to Haeberlen parameters. Note, for non-symmetric
    tensors, the conversion is applied by first symmetrizing the tensor
    using (tensor + tensor.T) / 2.

    Args:
        tensor:
            A 3x3 np.ndarray Cartesian tensor.
    Returns:
        A tuple of (euler_angles, zeta_sigma, eta_sigma, isotropic_component)
        where ``euler_angles`` is an ndarray  of three the Euler angles
        [``alpha``, ``beta``, ``gamma``] using the  "zyz" convention,
        ``zeta_sigma``, ``eta_sigma``, and ``isotropic_component`` are the
        corresponding Haeberlen anisotropy, asymmetry, and isotropic
        parameters.
    """
    tensor = (tensor + tensor.T) / 2  # Make sure the tensor is symmetric
    isotropic_component = np.trace(tensor) / 3
    traceless_tensor = tensor - isotropic_component * np.eye(3)
    # Calculate the eigenvalues and eigenvectors of the traceless tensor
    eigenvalues, L = np.linalg.eigh(traceless_tensor)
    idx = np.argsort(
        np.abs(eigenvalues)
    )  # Sort the eigenvalues by their absolute values
    idx[0], idx[1] = idx[1], idx[0]  # Swap the first two indices
    eigenvalues = eigenvalues[idx]  # Swap the eigenvalues
    L = L[:, idx]  # Swap the columns of the eigenvector matrix
    # Flip the sign of the eigenvectors if the determinant is negative
    if np.linalg.det(L) < 0:
        L = -L

    r = R.from_matrix(L)
    euler_angles = r.as_euler("zyz")

    zeta_sigma = eigenvalues[2]
    eta_sigma = (eigenvalues[1] - eigenvalues[0]) / eigenvalues[2]

    return euler_angles, zeta_sigma, eta_sigma, isotropic_component


def from_haeberlen_params(
    euler_angles: List[float],
    zeta_sigma: float,
    eta_sigma: float,
    isotropic_component: float,
) -> np.ndarray:
    """Haeberlen parameters to a 3x3 symmetric Cartesian tensor.

    Args:
        euler_angles:
            An ndarray of three Euler angles [``alpha``, ``beta``, ``gamma``]
            using the "zyz" convention.
        zeta_sigma:
            Anisotropy parameter.
        eta_sigma:
            Asymmetry parameter
        isotropic_component:
            Isotropic parameter

    Returns:
        tensor: A 3x3 np.ndarray of symmetric Cartesian tensor.
    """
    eigenvalues = np.array(
        [
            -zeta_sigma * (1 + eta_sigma) / 2,
            -zeta_sigma * (1 - eta_sigma) / 2,
            zeta_sigma,
        ]
    )

    # Assuming alpha_s, beta_s, gamma_s are your Euler angles
    r = R.from_euler("zyz", euler_angles)
    L = r.as_matrix()
    diag_matrix = np.diag(eigenvalues)

    tensor = L @ diag_matrix @ L.T

    tensor += isotropic_component * np.eye(3)
    return tensor


def mehring_principal_components_to_maryland(
    lambdas: List[float],
) -> Tuple[float, float, float]:
    """Mehring principal components to  Maryland.

    Args:
        lambdas:
            Eigenvalues from Mehring convention.

    Returns:
        A tuple of (span, skew, isotropic) in Maryland convention where
        ``span`` is the width, ``skew`` is the asymmetry, and ``isotropic``
        is the isotropic components.
    """
    isotropic = np.mean(lambdas)
    span = lambdas[2] - lambdas[0]
    skew = 3 * (isotropic - lambdas[1]) / span
    return span, skew, isotropic


def maryland_to_mehring_principal_components(
    isotropic: float, span: float, skew: float
) -> np.ndarray:
    """Maryland components to Mehring principal components.

    Args:
        span:
            Maryland anisotropy component.
        skew:
            Maryland asymmetry component.
        isotropic:
            Maryland isotropic component.

    Returns:
        A ndarray of Mehring Eigenvalues.
    """
    lambda2 = isotropic - skew * span / 3
    lambda3 = (span + 2 * isotropic + skew * span / 3) / 2
    lambda1 = lambda3 - span
    return np.array([lambda1, lambda2, lambda3])


def maryland_to_haeberlen_params(
    isotropic: float, span: float, skew: float
) -> Tuple[float, float, float]:
    lambdas = maryland_to_mehring_principal_components(isotropic, span, skew)
    (
        euler_angles,
        zeta_sigma,
        eta_sigma,
        isotropic_component,
    ) = to_haeberlen_params(np.diag(lambdas))
    return zeta_sigma, eta_sigma, isotropic_component


def haeberlen_params_to_maryland(
    zeta_sigma: float, eta_sigma: float, isotropic_component: float
) -> Tuple[float, float, float]:
    tensor = from_haeberlen_params(
        np.array([0, 0, 0]), zeta_sigma, eta_sigma, isotropic_component
    )
    mehring_eigenvals = to_mehring_params(tensor)[1]
    return mehring_principal_components_to_maryland(mehring_eigenvals)


# def zeta_and_eta_to_xy(zeta: float, eta: float) -> Tuple[float, float]:
#     r = abs(zeta)
#     if zeta <= 0:
#         theta = np.pi / 4 * eta
#     else:
#         theta = np.pi / 2 * (1 - eta / 2)

#     x = r * np.cos(theta)
#     y = r * np.sin(theta)
#     return x, y


# Dipolar coupling constant is given by:
# D = -𝛾_1 𝛾_2 µ_0 ℏ/(8*π^2 r^3) in Hz
# where r is the distance in m.
# Here, we use
# D = -𝛾_1' * 𝛾_2' * µ_0*ℏ*1E12*1E30/(2*R^3) in Hz
# where 𝛾_1' and 𝛾_2' are the reduced gyromagnetic ratios
# of the two isotopes in MHz/T, and R is the distance in Å.
# -µ_0*ℏ*1E12*1E30/2 = -66.2607015 m^3•kg^2/(s^3•A^2)
def dipolar_coupling_constant(
    isotope_symbol_1: str, isotope_symbol_2: str, distance: float
):
    """Dipolar coupling constant between two isotopes a distance apart

    Args:
        isotope_symbol_1:
            A string of any of the mrsimulator allowed isotopes for isotope 1.
        isotope_symbol_2
            A string of any of the mrsimulator allowed isotopes for isotope 2.
        distance:
            Distance between the isotopes in units of Angstrom.

    Return:
        Dipolar coupling constant in units of Hz.
    """
    isotope_1 = Isotope(symbol=isotope_symbol_1).gyromagnetic_ratio
    isotope_2 = Isotope(symbol=isotope_symbol_2).gyromagnetic_ratio
    return -66.2607015 * isotope_1 * isotope_2 / (distance) ** 3


def calculate_D_tensor(r1: Union[list, np.ndarray], r2: Union[list, np.ndarray]):
    # Compute the vector r
    r1 = np.array(r1)
    r2 = np.array(r2)
    r = r1 - r2
    x, y, z = r
    r_magnitude = np.sqrt(x**2 + y**2 + z**2)
    r_squared = r_magnitude**2
    r_cubed = r_magnitude**3

    # Initialize the D tensor
    D = np.zeros((3, 3))
    r_components = np.array([x, y, z])

    # Compute the components of the D tensor
    for i in range(3):
        for k in range(3):
            D[i, k] = (1 / r_cubed) * (
                3 * r_components[i] * r_components[k] / r_squared - (1 if i == k else 0)
            )

    return D


def dipolar_tensor(site_1: list, site_2: list) -> np.ndarray:
    """Generates a 3x3 symmetric cartesian tensor from isotope site coordinates.

    Args:
        site_1:
            A list of (isotope_symbol, site_coordiantes) for site 1, where
            ``isotope_symbol`` is a string of any of the mrsimulator allowed
            isotopes and ``site_coordiantes`` is a list or ndarray of
            (x, y, z) site coordinates in units of Angstrom.
        site_2:
            A list of (isotope_symbol, site_coordiantes) for site 2, where
            ``isotope_symbol`` is a string of any of the mrsimulator allowed
            isotopes and ``site_coordiantes`` is a list or ndarray of
            (x, y, z) site coordinates in units of Angstrom.

    Returns:
        A 3x3 ndarray of symmetric dipolar Cartesian tensor in units of Hz.

    Example:
        >>> site1_coords = [2.1, 3.1, 1.3]  # coords for 1H in A
        >>> site2_coords = [0, 0, 1.2]  # coords for 13C in A
        >>> dipole_tensor = dipolar_tensor(
        ...     site_1=['1H', site1_coords],
        ...     site_2=['13C', site2_coords]
        ... )
    """
    isotope_1 = Isotope(symbol=site_1[0]).gyromagnetic_ratio
    isotope_2 = Isotope(symbol=site_2[0]).gyromagnetic_ratio
    D_tensor = calculate_D_tensor(site_1[1], site_2[1])
    D_tensor *= -66.2607015 * isotope_1 * isotope_2 / 2
    return D_tensor


def to_symmetric_tensor(tensor: np.ndarray, type: str = "shielding") -> SymmetricTensor:
    """Cartesian 3x3 tensor to mrsimulator SymmetricTensor object. Note that
    only the traceless symmetric part of the tensor is converted to the
    SymmetricTensor object.

    Args:
        tensor:
            A 3x3 np.ndarray Cartesian tensor
        type:
            String with any one of the allowed listings. [``"shielding"``,
            ``"j_coupling"``, ``"quadrupolar"``, ``"dipolar"``]

    Returns:
        A :ref:`sy_api` object.

    Example:
        >>> tensor = np.array([
        ...     [-8.05713333, -1.4523, 35.7252],
        ...     [-5.5725, 26.38916667, -5.2804],
        ...     [33.1405, -0.6241, -18.33203333],
        ... ])
        >>> symmetric_tensor = to_symmetric_tensor(tensor, type="shielding")
    """
    euler_angles, zeta, eta, _ = to_haeberlen_params(tensor)

    if type == "shielding" or type == "j_coupling":
        return SymmetricTensor(
            zeta=zeta,
            eta=eta,
            alpha=euler_angles[0],
            beta=euler_angles[1],
            gamma=euler_angles[2],
        )
    elif type == "quadrupolar":
        return SymmetricTensor(
            Cq=zeta,
            eta=eta,
            alpha=euler_angles[0],
            beta=euler_angles[1],
            gamma=euler_angles[2],
        )
    elif type == "dipolar":
        return SymmetricTensor(
            D=zeta,
            eta=eta,
            alpha=euler_angles[0],
            beta=euler_angles[1],
            gamma=euler_angles[2],
        )
    else:
        raise ValueError(f"Unknown tensor type: {type}")
