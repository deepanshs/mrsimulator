from typing import ClassVar
from typing import Dict

import numpy as np

from ._base import ModuleOperation

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


class Baseline(ModuleOperation):
    module_name: ClassVar[str] = __name__

    @property
    def function(self):
        return "baseline"


class ConstantOffset(Baseline):
    r"""Add an offset to the dependent variables (y) of the CSDM object.

    .. math::
        y += \text{offset}

    where :math:`y` is the CSDM dependent variable.

    Args:
        float offset: The offset factor. The default value is 0.

    Example
    -------

    >>> from mrsimulator import signal_processor as sp
    >>> operation1 = sp.baseline.ConstantOffset(offset=20)
    """

    offset: float = 0

    def operate(self, dataset):
        """Applies the operation for which the class is named for.

        Args:
            dataset: CSDM object
        """
        dataset += self.offset
        return dataset


class Polynomial(Baseline):
    r"""Add a baseline polynomial to all dependent variables (y) in the CSDM object.

    The baseline polynomial function is

    .. math::
            f(x) = \sum_{i=0}^n c_i \times x^i,

    where :math:`c_i` are the coefficients corresponding to :math:`x^i`.

    Args:
        Dict polynomial_dictionary: A dictionary of the form {'ci': coef}, where i
            represents the i-th order polynomial term and 'coef' is the leading
            coefficient for the i-th term. For example :math:`4x^2 + 5` would be
            supplied as {'c2': 4, 'c0': 5}

        int dim_index: The index of the CSDM dimension along which the operation is
            applied. The default is the dimension at index 0.

    Example
    -------

    >>> from mrsimulator import signal_processor as sp
    >>> operation1 = sp.baseline.Polynomial(polynomial_dictionary = {'c0':10, 'c2':2})
    """

    polynomial_dictionary: Dict = {}
    dim_index: int = 0
    # property_units: Dict = {"default": "Hz"}

    def operate(self, dataset):
        """Applies the operation.

        Args:
            dataset: CSDM object
        """
        x = dataset.dimensions[self.dim_index].coordinates
        d1 = x.value
        fn = np.zeros(len(x))
        for key, val in self.polynomial_dictionary.items():
            exponent = key.split("c")[-1]
            fn += float(val) * d1 ** int(exponent)
        for item in dataset.y:
            item.components += fn
        return dataset
