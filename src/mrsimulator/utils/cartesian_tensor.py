from typing import List
from typing import Tuple

import numpy as np
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
