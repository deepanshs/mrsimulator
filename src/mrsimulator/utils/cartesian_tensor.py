from typing import List
from typing import Tuple

import numpy as np
from mrsimulator.spin_system.isotope import Isotope
from scipy.spatial.transform import Rotation as R

__author__ = "Philip Grandinetti"
__email__ = "grandinetti.1@osu.edu"


def to_mehring_parameters(tensor: np.ndarray) -> Tuple[np.ndarray, np.ndarray]:
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


def from_mehring_parameters(
    euler_angles: List[float], eigenvalues: List[float]
) -> np.ndarray:
    # Assuming alpha_s, beta_s, gamma_s are your Euler angles
    r = R.from_euler("zyz", euler_angles)
    L = r.as_matrix()

    diag_matrix = np.diag(eigenvalues)
    tensor = L @ diag_matrix @ L.T
    return tensor


def to_haeberlen_parameters(
    tensor: np.ndarray,
) -> Tuple[np.ndarray, float, float, float]:
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


def from_haeberlen_parameters(
    euler_angles: List[float],
    zeta_sigma: float,
    eta_sigma: float,
    isotropic_component: float,
) -> np.ndarray:
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
    isotropic = np.mean(lambdas)
    span = lambdas[2] - lambdas[0]
    skew = 3 * (isotropic - lambdas[1]) / span
    return span, skew, isotropic


def maryland_to_mehring_principal_components(
    isotropic: float, span: float, skew: float
) -> np.ndarray:
    lambda1 = isotropic - span
    lambda2 = isotropic - skew * span / 3
    lambda3 = (span + 2 * isotropic + skew * span / 3) / 2
    return np.array([lambda1, lambda2, lambda3])


def maryland_to_haeberlen_parameters(
    isotropic: float, span: float, skew: float
) -> Tuple[float, float, float]:
    lambdas = maryland_to_mehring_principal_components(isotropic, span, skew)
    (
        euler_angles,
        zeta_sigma,
        eta_sigma,
        isotropic_component,
    ) = to_haeberlen_parameters(lambdas)
    return zeta_sigma, eta_sigma, isotropic_component


def haeberlen_parameters_to_maryland(
    zeta_sigma: float, eta_sigma: float, isotropic_component: float
) -> Tuple[float, float, float]:
    tensor = from_haeberlen_parameters(
        np.array([0, 0, 0]), zeta_sigma, eta_sigma, isotropic_component
    )
    return mehring_principal_components_to_maryland(to_mehring_parameters(tensor))


def zeta_and_eta_to_xy(zeta: float, eta: float) -> Tuple[float, float]:
    r = abs(zeta)
    if zeta <= 0:
        theta = np.pi / 4 * eta
    else:
        theta = np.pi / 2 * (1 - eta / 2)

    x = r * np.cos(theta)
    y = r * np.sin(theta)
    return x, y


# Dipolar coupling constant is given by:
# D = -𝛾_1 𝛾_2 µ_0 ℏ/(8*π^2 r^3) in Hz
# where r is distance in m.
# Here, we use
# D = -𝛾_1' * 𝛾_2' * µ_0*ℏ*1E12*1E30/(2*R^3) in Hz
# where 𝛾_1' and 𝛾_2' are the reduced gyromagnetic ratios
# of the two isotopes in MHz/T, and R is the distance in Å.
# -µ_0*ℏ*1E12*1E30/2 = -66.2607015 m^3•kg^2/(s^3•A^2)
def dipolar_coupling_constant(isotope_symbol_1, isotope_symbol_2, distance: float):
    isotope_1 = Isotope(symbol=isotope_symbol_1).gyromagnetic_ratio
    isotope_2 = Isotope(symbol=isotope_symbol_2).gyromagnetic_ratio
    return -66.2607015 * isotope_1 * isotope_2 / (distance) ** 3
