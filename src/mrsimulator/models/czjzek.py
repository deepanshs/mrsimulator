"""
Analytical czjzek distribution on polar and non-polar grid

__author__ = "Deepansh J. Srivastava"
__email__ = "dsrivastava@hyperfine.io"
"""
import csdmpy as cp
import mrsimulator.models.analytical_distributions as analytical_dist
import numpy as np
from mrsimulator.base_model import get_Haeberlen_components
from mrsimulator.clib import histogram2d
from mrsimulator.spin_system.tensors import SymmetricTensor

from .utils import get_principal_components
from .utils import zeta_eta_to_x_y

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"

ANALYTICAL_AVAILABLE = {"czjzek": analytical_dist.czjzek}


def _czjzek_random_distribution_tensors(sigma, n):
    r"""Czjzek random distribution model.

    Args:
        float sigma: The standard deviation of the five-dimensional multi-variate normal
            distribution.
        int n: Number of samples drawn from the Czjzek random distribution model.

    Description
    -----------

    U is an array of the coordinates randomly drawn from an uncorrelated
    five-dimensional multivariate normal distribution with standard deviation `sigma`
    and zero mean.

    The components of the traceless second-rank symmetric cartesian tensor, S_ij,
    follows

    Sxx = sqrt(3) * U5 - U1

    Sxy = Syx = sqrt(3) * U4

    Syy = -sqrt(3) * U5 - U1

    Sxz = Szx = sqrt(3) * U2

    Szz = 2 * U1

    Syz = Szy = sqrt(3) * U3
    """
    u1, u2, u3, u4, u5 = czjzek_random_components(sigma, n)

    # Create N random tensors
    tensors = np.empty((n, 3, 3))  # n x 3 x 3 tensors

    tensors[:, 0, 0] = u5 - u1  # xx
    # tensors[:, 0, 1] = U4  # xy
    # tensors[:, 0, 2] = U2  # xz

    tensors[:, 1, 0] = u4  # yx
    tensors[:, 1, 1] = -u5 - u1  # yy
    # tensors[:, 1, 2] = U3  # yz

    tensors[:, 2, 0] = u2  # zx
    tensors[:, 2, 1] = u3  # zy
    tensors[:, 2, 2] = 2 * u1  # zz

    return tensors


def czjzek_random_components(sigma: float, n: int):
    """Five dimensional random components of a 2nd rank symmetric tensor of size n.

    Args:
        float sigma: The standard deviation of the five-dimensional multi-variate normal
            distribution.
        int n: Number of samples drawn from the Czjzek random distribution model.
    """
    # The random sampling U1, U2, ... U5
    u1 = np.random.normal(0.0, sigma, n)

    sqrt_3_sigma = np.sqrt(3) * sigma
    u2 = np.random.normal(0.0, sqrt_3_sigma, n)
    u3 = np.random.normal(0.0, sqrt_3_sigma, n)
    u4 = np.random.normal(0.0, sqrt_3_sigma, n)
    u5 = np.random.normal(0.0, sqrt_3_sigma, n)
    return u1, u2, u3, u4, u5


def get_p_q_basis(u1, u2, u3, u4, u5):
    """Evalute the one expression base"""
    temp_q = u2**2 + u3**2 - 2 * u4**2 - 2 * u5**2
    temp_q2 = u2**2 - u3**2

    q_basis = np.array(
        [
            -u1 * (2 * u1**2 + temp_q) - u5 * temp_q2 - 2 * u2 * u3 * u4,  # rho^3
            -3 * u1**2 - 0.5 * temp_q,  # zeta rho^2
            -2.0 * u1 * u5 + 0.5 * temp_q2,  # zeta eta rho^2
            0.5 * u1,  # zeta^2 rho (-3 + eta^2)
            -0.25 * np.ones(u1.size),  # zeta^3 (1 - eta^2)
            -u5,  # zeta^2 eta rho
        ]
    )

    p_basis = np.array(
        [
            -3 * u1**2 - u2**2 - u3**2 - u4**2 - u5**2,  # rho^2
            -3 * u1,  # zeta rho
            u5,  # zeta eta rho
            -0.25 * np.ones(u1.size),  # zeta^2 (3.0 + eta^2)
        ]
    )
    return -p_basis / 3.0, -q_basis / 2.0


class AbstractDistribution:
    """Abstract distribution"""

    model_name = "base"

    def __init__(
        self,
        mean_isotropic_chemical_shift=0.0,
        abundance=100.0,
        polar=False,
        cache_tensors=False,
    ):
        """Basic class attributes for distributions"""
        self._cache_tensors = cache_tensors
        self._basis_p_q = None
        self.mean_isotropic_chemical_shift = mean_isotropic_chemical_shift
        self.abundance = abundance
        self.polar = polar

    def pdf(
        self,
        pos,
        size: int = 400000,
        analytical: bool = True,
        pack_as_csdm: bool = False,
    ):
        """Generates a probability distribution function by binning the random
        variates of length size onto the given grid system.

        Args:
            pos: A list of coordinates along the two dimensions given as NumPy arrays.
            size: The number of random variates drawn in generating the pdf. The default
                is 400000.
            pack_as_csdm: If true, returns as csdm object.

        Returns:
            A list of x and y coordinates and the corresponding amplitudes if not packed
            as csdm object, else a csdm object.

        Example:
            >>> import numpy as np
            >>> cq = np.arange(50) - 25
            >>> eta = np.arange(21)/20
            >>> amp = cz_model.pdf(pos=[cq, eta])  # returns amp as a CSDM object.
        """
        if analytical and self.model_name in ANALYTICAL_AVAILABLE:
            analytical_model = ANALYTICAL_AVAILABLE[self.model_name]
            pos_a, pos_b, data_ = analytical_model(self.sigma, pos, self.polar)
        else:
            pos_a, pos_b, data_ = self.pdf_numerical(pos, size)

        if pack_as_csdm:
            if len(pos_a.shape) == 1:
                return self.pack_csdm_object([pos_a, pos_b], data_)
            else:
                return self.pack_csdm_object([pos_a[0, :], pos_b[:, 0]], data_)
        return pos_a, pos_b, data_

    def pdf_numerical(self, pos, size: int = 400000):
        """Generate distribution numerically"""
        delta_z = (pos[0][1] - pos[0][0]) / 2
        delta_e = (pos[1][1] - pos[1][0]) / 2
        x = [pos[0][0] - delta_z, pos[0][-1] + delta_z]
        y = [pos[1][0] - delta_e, pos[1][-1] + delta_e]

        x_size = pos[0].size
        y_size = pos[1].size
        zeta, eta = self.rvs(size)

        _, _, hist = histogram2d(
            sample_x=zeta,
            sample_y=eta,
            weights=np.ones(zeta.size, dtype=float),
            x_count=x_size,
            y_count=y_size,
            x_min=x[0],
            x_max=x[1],
            y_min=y[0],
            y_max=y[1],
        )

        hist /= hist.sum()
        return pos[0], pos[1], hist.T

    def pack_csdm_object(self, pos, data):
        """Pack data and coordinates as csdm objects"""
        dims = [cp.as_dimension(item) for item in pos]
        dvs = cp.as_dependent_variable(data)
        return cp.CSDM(dimensions=dims, dependent_variables=[dvs])

    def add_lmfit_params(self, params, i):
        """Add lmfit params for base class"""
        prefix = self.model_name
        params.add(
            f"{prefix}_{i}_mean_isotropic_chemical_shift",
            value=self.mean_isotropic_chemical_shift,
        )
        params.add(f"{prefix}_{i}_abundance", value=self.abundance, min=0, max=100)

    def update_lmfit_params(self, params, i):
        """Update lmfit params for base class"""
        prefix = self.model_name
        self.mean_isotropic_chemical_shift = params[
            f"{prefix}_{i}_mean_isotropic_chemical_shift"
        ].value
        self.abundance = params[f"{prefix}_{i}_abundance"].value


class CzjzekDistribution(AbstractDistribution):
    r"""A Czjzek distribution model class.

    The Czjzek distribution model is a random sampling of second-rank traceless
    symmetric tensors whose explicit matrix form follows

    .. math::
        {\bf S} = \left[
        \begin{array}{l l l}
        \sqrt{3} U_5 - U_1   & \sqrt{3} U_4          & \sqrt{3} U_2 \\
        \sqrt{3} U_4         & -\sqrt{3} U_5 - U_1   & \sqrt{3} U_3 \\
        \sqrt{3} U_2         & \sqrt{3} U_3          & 2 U_1
        \end{array}
        \right],

    where the components, :math:`U_i`, are randomly drawn from a five-dimensional
    multivariate normal distribution. Each component, :math:`U_i`, is a dimension of
    the five-dimensional uncorrelated multivariate normal distribution with the mean
    of :math:`<U_i>=0` and the variance :math:`<U_iU_i>=\sigma^2`.

    .. math::
        S_T = S_C(\sigma),

    Args:
        float sigma: The Gaussian standard deviation.

    .. note:: In the original Czjzek paper, the parameter :math:`\sigma` is given as
        two times the standard deviation of the multi-variate normal distribution used
        here.

    Example:
        >>> from mrsimulator.models import CzjzekDistribution
        >>> cz_model = CzjzekDistribution(0.5)
    """

    model_name = "czjzek"

    def __init__(
        self,
        sigma: float,
        mean_isotropic_chemical_shift: float = 0.0,
        abundance: float = 100.0,
        polar=False,
        cache=True,
    ):
        super().__init__(
            cache_tensors=cache,
            polar=polar,
            mean_isotropic_chemical_shift=mean_isotropic_chemical_shift,
            abundance=abundance,
        )
        self.sigma = sigma

    def rvs(self, size: int):
        """Draw random variates of length `size` from the distribution.

        Args:
            size: The number of random points to draw.

        Returns:
            A list of two NumPy array, where the first and the second array are the
            anisotropic/quadrupolar coupling constant and asymmetry parameter,
            respectively.

        Example:
            >>> Cq_dist, eta_dist = cz_model.rvs(size=1000000)
        """
        if self._cache_tensors:
            if self._basis_p_q is None:
                self._basis_p_q = get_p_q_basis(
                    *czjzek_random_components(self.sigma, size)
                )
            basis_p_q = self._basis_p_q
        else:
            basis_p_q = get_p_q_basis(*czjzek_random_components(self.sigma, size))

        zeta, eta = get_Haeberlen_components(*basis_p_q, 0, 0, 1.0)
        if not self.polar:
            return zeta, eta
        return zeta_eta_to_x_y(zeta, eta)

    def add_lmfit_params(self, params, i):
        """Create lmfit params for index i"""
        prefix = self.model_name
        params.add(f"{prefix}_{i}_sigma", value=self.sigma, min=0)
        super().add_lmfit_params(params, i)
        return params

    def update_lmfit_params(self, params, i):
        """Update lmfit params for index i"""
        prefix = self.model_name
        self.sigma = params[f"{prefix}_{i}_sigma"].value
        super().update_lmfit_params(params, i)


class ExtCzjzekDistribution(AbstractDistribution):
    r"""An extended Czjzek distribution distribution model.

    The extended Czjzek random distribution [#f1]_ model is an extension of the Czjzek
    model, given as

    .. math::
        S_T = S(0) + \rho S_C(\sigma=1),

    where :math:`S_T` is the total tensor, :math:`S(0)` is the dominant tensor,
    :math:`S_C(\sigma=1)` is the Czjzek random model attributing to the random
    perturbation of the tensor about the dominant tensor, :math:`S(0)`, and
    :math:`\rho` is the size of the perturbation. Note, in the above equation, the
    :math:`\sigma` parameter from the Czjzek random model, :math:`S_C`, has no meaning
    and is set to one. The factor, :math:`\rho`, is defined as

    .. math::
        \rho = \frac{||S(0)|| \epsilon}{\sqrt{30}},

    where :math:`\|S(0)\|` is the 2-norm of the dominant tensor, and :math:`\epsilon`
    is a fraction.

    .. [#f1] Gérard Le Caër, Bruno Bureau, and Dominique Massiot,
        An extension of the Czjzek model for the distributions of electric field
        gradients in disordered solids and an application to NMR spectra of 71Ga in
        chalcogenide glasses. Journal of Physics: Condensed Matter, 2010, 22, 065402.
        DOI: 10.1088/0953-8984/22/6/065402

    Args:
        SymmetricTensor symmetric_tensor: A shielding or quadrupolar symmetric tensor
            or equivalent dict object.
        float eps: A fraction determining the extent of perturbation.

    Example
    -------

    >>> from mrsimulator.models import ExtCzjzekDistribution
    >>> S0 = {"Cq": 1e6, "eta": 0.3}
    >>> ext_cz_model = ExtCzjzekDistribution(S0, eps=0.35)
    """

    model_name = "ext_czjzek"

    def __init__(
        self,
        symmetric_tensor: SymmetricTensor,
        eps: float,
        mean_isotropic_chemical_shift: float = 0.0,
        abundance: float = 100.0,
        polar=False,
        cache=True,
    ):
        super().__init__(
            cache_tensors=cache,
            polar=polar,
            mean_isotropic_chemical_shift=mean_isotropic_chemical_shift,
            abundance=abundance,
        )
        if isinstance(symmetric_tensor, dict):
            self.symmetric_tensor = SymmetricTensor(**symmetric_tensor)
        else:
            self.symmetric_tensor = symmetric_tensor
        self.eps = eps

    def rvs(self, size: int):
        """Draw random variates of length `size` from the distribution.

        Args:
            size: The number of random points to draw.

        Returns:
            A list of two NumPy array, where the first and the second array are the
            anisotropic/quadrupolar coupling constant and asymmetry parameter,
            respectively.

        Example:
            >>> Cq_dist, eta_dist = ext_cz_model.rvs(size=1000000)
        """

        # czjzek_random_distribution model
        if self._cache_tensors:
            if self._basis_p_q is None:
                self._basis_p_q = get_p_q_basis(*czjzek_random_components(1, size))
            basis_p_q = self._basis_p_q
        else:
            basis_p_q = get_p_q_basis(*czjzek_random_components(1, size))

        symmetric_tensor = self.symmetric_tensor

        zeta = symmetric_tensor.zeta or symmetric_tensor.Cq
        eta = symmetric_tensor.eta

        # the traceless second-rank symmetric Cartesian tensor in PAS
        T0 = [0.0, 0.0, 0.0]
        if zeta != 0 and eta != 0:
            T0 = get_principal_components(zeta, eta)

        # 2-norm of the tensor
        norm_T0 = np.linalg.norm(T0)

        # the perturbation factor # np.sqrt(30) = 5.477225575
        rho = self.eps * norm_T0 / 5.4772255751

        zeta, eta = get_Haeberlen_components(*basis_p_q, zeta, eta, rho)
        if not self.polar:
            return zeta, eta
        return zeta_eta_to_x_y(zeta, eta)

    def add_lmfit_params(self, params, i):
        """Create lmfit params for index i"""
        prefix = self.model_name
        if self.symmetric_tensor.zeta is not None:
            zeta = self.symmetric_tensor.zeta
        else:
            zeta = self.symmetric_tensor.Cq
        params.add(f"{prefix}_{i}_symmetric_tensor_zeta", value=zeta)
        params.add(
            f"{prefix}_{i}_symmetric_tensor_eta",
            value=self.symmetric_tensor.eta,
            min=0,
            max=1,
        )
        params.add(f"{prefix}_{i}_eps", value=self.eps, min=1e-6)
        super().add_lmfit_params(params, i)
        return params

    def update_lmfit_params(self, params, i):
        """Create lmfit params for index i"""
        prefix = self.model_name

        zeta = params[f"{prefix}_{i}_symmetric_tensor_zeta"].value
        if self.symmetric_tensor.zeta is not None:
            self.symmetric_tensor.zeta = zeta
        else:
            self.symmetric_tensor.Cq = zeta

        self.symmetric_tensor.eta = params[f"{prefix}_{i}_symmetric_tensor_eta"].value
        self.eps = params[f"{prefix}_{i}_eps"].value
        super().update_lmfit_params(params, i)
