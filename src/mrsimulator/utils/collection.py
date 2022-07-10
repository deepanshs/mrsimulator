from typing import Dict
from typing import List
from typing import Union

import numpy as np
from mrsimulator import Site
from mrsimulator import SpinSystem


__author__ = ["Deepansh Srivastava", "Matthew D. Giammar"]
__email__ = ["srivastava.89@osu.edu", "giammar.7@buckeyemail.osu.edu"]

SHIELDING_SYM_PARAMS = ["zeta", "eta", "alpha", "beta", "gamma"]
QUADRUPOLAR_PARAMS = ["Cq", "eta", "alpha", "beta", "gamma"]
LIST_LEN_ERROR_MSG = (
    "All arguments must be the same size. If one attribute is a type list of length n, "
    "then all attributes with list types must also be of length n, and all remaining "
    "attributes must be scalar (singular float, int, or str)."
)


def single_site_system_generator(
    isotope: Union[str, List[str]],
    isotropic_chemical_shift: Union[float, List[float], np.ndarray] = 0,
    shielding_symmetric: Dict = None,
    shielding_antisymmetric: Dict = None,
    quadrupolar: Dict = None,
    abundance: Union[float, List[float], np.ndarray] = None,
    site_name: Union[str, List[str]] = None,
    site_label: Union[str, List[str]] = None,
    site_description: Union[str, List[str]] = None,
    rtol: float = 1e-3,
) -> List[SpinSystem]:
    r"""Generate and return a list of single-site spin systems from the input parameters

    Args:
        isotope:
            A required string or a list of site isotopes.
        isotropic_chemical_shift:
            A float or a list/ndarray of isotropic chemical shifts per site per spin
            system. The default is 0.
        shielding_symmetric:
            A shielding symmetric dict object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        shielding_antisymmetric:
            A shielding antisymmetric dict object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``alpha``, and ``beta``.
        quadrupolar:
            A quadrupolar dict object, where the keyword value can either be a float or
            a list/ndarray of floats. The default value is None. The allowed keywords
            are ``Cq``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        abundance:
            A float or a list/ndarray of floats describing the abundance of each spin
            system.
        site_name:
            A string or a list of strings with site names per site per spin system. The
            default is None.
        site_label:
            A string or a list of strings with site labels per site per spin system. The
            default is None.
        site_description:
            A string or a list of strings with site descriptions per site per spin
            system. The default is None.
        rtol:
            The relative tolerance used in determining the cutoff abundance, given as,
            :math:`\tt{abundance}_{\tt{cutoff}} = \tt{rtol} * \tt{max(abundance)}.`
            The spin systems with abundance below this threshold are ignored.

    Returns:
        List of :ref:`spin_sys_api` objects with a single :ref:`site_api`

    Example:
        **Single spin system:**

        >>> sys1 = single_site_system_generator(
        ...     isotope=["1H"],
        ...     isotropic_chemical_shift=10,
        ...     site_name="Single Proton",
        ... )
        >>> print(len(sys1))
        1

        **Multiple spin system:**

        >>> sys2 = single_site_system_generator(
        ...     isotope="1H",
        ...     isotropic_chemical_shift=[10] * 5,
        ...     site_name="5 Protons",
        ... )
        >>> print(len(sys2))
        5

        **Multiple spin system with dictionary arguments:**

        >>> Cq = [4.2e6] * 12
        >>> sys3 = single_site_system_generator(
        ...     isotope="17O",
        ...     isotropic_chemical_shift=60.0,  # in ppm,
        ...     quadrupolar={"Cq": Cq, "eta": 0.5},  # Cq in Hz
        ... )
        >>> print(len(sys3))
        12

    .. note::
        The parameter value can either be a float or a list/ndarray. If the parameter
        value is a float, the given value is assigned to the respective parameter in all
        the spin systems. If the parameter value is a list or ndarray, its `ith` value
        is assigned to the respective parameter of the `ith` spin system. When multiple
        parameter values are given as lists/ndarrays, the length of all the lists must
        be the same.
    """
    isotope = _flatten_item(isotope)
    isotropic_chemical_shift = _flatten_item(isotropic_chemical_shift)
    shielding_symmetric = _flatten_item(shielding_symmetric)
    shielding_antisymmetric = _flatten_item(shielding_antisymmetric)
    quadrupolar = _flatten_item(quadrupolar)
    site_name = _flatten_item(site_name)
    site_label = _flatten_item(site_label)
    site_description = _flatten_item(site_description)
    abundance = _flatten_item(abundance)
    args = [
        isotope,
        isotropic_chemical_shift,
        shielding_symmetric,
        shielding_antisymmetric,
        quadrupolar,
        site_name,
        site_label,
        site_description,
        abundance,
    ]

    n_sites = _check_lengths_of_args(*args)
    sites = _site_generator(n_sites, *args[:-1])  # Don't pass abundance

    if abundance is None:
        abundance = np.asarray([1 / n_sites] * n_sites)
    if isinstance(abundance, (int, float, np.floating)):
        abundance = np.asarray([abundance] * n_sites)

    keep_idxs = np.asarray(abundance > rtol * abundance.max()).nonzero()[0]

    return [
        SpinSystem(sites=[site], abundance=abd)
        for i, (site, abd) in enumerate(zip(sites, abundance))
        if i in keep_idxs
    ]


def site_generator(
    isotope: Union[str, List[str]],
    isotropic_chemical_shift: Union[float, List[float], np.ndarray] = 0,
    shielding_symmetric: Dict = None,
    shielding_antisymmetric: Dict = None,
    quadrupolar: Dict = None,
    name: Union[str, List[str]] = None,
    label: Union[str, List[str]] = None,
    description: Union[str, List[str]] = None,
) -> List[Site]:
    r"""Generate a list of Site objects from lists of site attributes.

    Args:
        isotope:
            A required string or a list of site isotopes.
        isotropic_chemical_shift:
            A float or a list/ndarray of isotropic chemical shifts per site. The default
            is 0.
        shielding_symmetric:
            A shielding symmetric dict object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        shielding_antisymmetric:
            A shielding antisymmetric dict object, where the keyword value can either
            be a float or a list/ndarray of floats. The default value is None. The
            allowed keywords are ``zeta``, ``alpha``, and ``beta``.
        quadrupolar:
            A quadrupolar dict object, where the keyword value can either be a float or
            a list/ndarray of floats. The default value is None. The allowed keywords
            are ``Cq``, ``eta``, ``alpha``, ``beta``, and ``gamma``.
        name:
            A string or a list of strings with site names per site. The default is None.
        label:
            A string or a list of strings with site labels per site. The default is
            None.
        description:
            A string or a list of strings with site descriptions per site. The default
            is None.

    Returns:
        sites: List of :ref:`site_api` objects

    Example:
        **Generating 10 hydrogen sites:**

        >>> sites1 = site_generator(
        ...     isotope=["1H"] * 10,
        ...     isotropic_chemical_shift=-15,
        ...     name="10 Protons",
        ... )
        >>> print(len(sites1))
        10

        **Generating 10 hydrogen sites with different shifts:**

        >>> shifts = np.arange(-10, 10, 2)
        >>> sites2 = site_generator(
        ...     isotope=["1H"] * 10,
        ...     isotropic_chemical_shift=shifts,
        ...     name="10 Proton",
        ... )
        >>> print(len(sites2))
        10

        **Generating multiple sites with dictionary arguments:**

        >>> Cq = [4.2e6] * 12
        >>> sys3 = site_generator(
        ...     isotope="17O",
        ...     isotropic_chemical_shift=60.0,  # in ppm,
        ...     quadrupolar={"Cq": Cq, "eta": 0.5},  # Cq in Hz
        ... )
        >>> print(len(sys3))
        12
    """
    isotope = _flatten_item(isotope)
    isotropic_chemical_shift = _flatten_item(isotropic_chemical_shift)
    shielding_symmetric = _flatten_item(shielding_symmetric)
    shielding_antisymmetric = _flatten_item(shielding_antisymmetric)
    quadrupolar = _flatten_item(quadrupolar)
    name = _flatten_item(name)
    label = _flatten_item(label)
    description = _flatten_item(description)
    args = [
        isotope,
        isotropic_chemical_shift,
        shielding_symmetric,
        shielding_antisymmetric,
        quadrupolar,
        name,
        label,
        description,
    ]

    n_sites = _check_lengths_of_args(*args)

    return list(_site_generator(n_sites, *args))


def _site_generator(  # noqa: C901
    n_sites: int,
    isotope: Union[str, np.ndarray],
    isotropic_chemical_shift: Union[float, np.ndarray] = 0,
    shielding_symmetric: Dict = None,
    shielding_antisymmetric: Dict = None,
    quadrupolar: Dict = None,
    name: Union[str, np.ndarray] = None,
    label: Union[str, np.ndarray] = None,
    description: Union[str, np.ndarray] = None,
):
    r"""A generator function which returns :ref:`site_api` objects in a memory efficient
    manner.

    When the next site is requested from the generator, each site argument is
    checked to see if the argument is an array or not a list (will be float or str). If
    the argument is an array, the item at the next index in the array is passed,
    otherwise the constant float/str is passed.

    Although there are many conditional statements in the method, this approach ensures
    multiple arrays of length n_sites aren't temporarily needed reducing peak memory
    usage when building large numbers of sites.
    """
    for index in range(n_sites):
        yield Site(
            isotope=isotope[index] if isinstance(isotope, np.ndarray) else isotope,
            isotropic_chemical_shift=isotropic_chemical_shift[index]
            if isinstance(isotropic_chemical_shift, np.ndarray)
            else isotropic_chemical_shift,
            name=name[index] if isinstance(name, np.ndarray) else name,
            label=label[index] if isinstance(label, np.ndarray) else label,
            description=description[index]
            if isinstance(description, np.ndarray)
            else description,
            shielding_symmetric={
                key: val[index] if isinstance(val, np.ndarray) else val
                for key, val in shielding_symmetric.items()
            }
            if shielding_symmetric is not None
            else None,
            shielding_antisymmetric={
                key: val[index] if isinstance(val, np.ndarray) else val
                for key, val in shielding_antisymmetric.items()
            }
            if shielding_antisymmetric is not None
            else None,
            quadrupolar={
                key: val[index] if isinstance(val, np.ndarray) else val
                for key, val in quadrupolar.items()
            }
            if quadrupolar is not None
            else None,
        )


def _flatten_item(item):
    """Flattens multi-dimensional arrays to long array"""
    if isinstance(item, dict):
        return {
            k: np.ravel(np.asarray(v)) if isinstance(v, (list, np.ndarray)) else v
            for k, v in item.items()
        }
    if isinstance(item, (list, np.ndarray)):
        return np.ravel(np.asarray(item))
    return item


def _check_lengths_of_args(*args):
    """Raises error when not all lists are same length. Returns length on success"""
    args = [arg for arg in args if arg is not None]
    dicts = [
        list(args.pop(i).values())
        for i, arg in reversed(list(enumerate(args)))
        if isinstance(arg, dict)
    ]
    dicts = dicts[0] if dicts != [] else []
    length = {len(lst) for lst in args + dicts if isinstance(lst, (list, np.ndarray))}
    if len(length) != 1:
        raise ValueError(
            f"Not all arrays/lists passed were of the same length. {LIST_LEN_ERROR_MSG}"
        )
    return length.pop()
