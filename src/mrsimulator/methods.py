# -*- coding: utf-8 -*-
from copy import deepcopy

from .method import Event
from .method import Method
from .method import SpectralDimension

# program to create class dynamically


def generate_method_from_template(template, __doc__):

    # constructor
    def constructor(self, dimensions, parse=False, **kwargs):
        prep = {
            "channels": kwargs["channels"],
            "name": template["name"],
            "description": template["description"],
        }

        kwargs.pop("channels")
        global_events = template["global_event_attributes"]
        ge = set(global_events)
        kw = set(kwargs)
        common = kw.intersection(ge)

        if common != set():
            raise AttributeError(
                f"The attribute(s) {common} is fixed for {prep['name']} method."
            )

        dim = []
        for i, s in enumerate(template["spectral_dimensions"]):
            events = []
            for e in s["events"]:
                ew = set(e)
                intersection = ew.intersection(kw)
                kw = kwargs.copy()
                for item in intersection:
                    kw.pop(item)
                if parse:
                    events.append({**e, **kw, **global_events})
                else:
                    events.append(Event(**e, **kw, **global_events))
            if parse:
                dim.append({**dimensions[i], "events": events})
            else:
                dim.append(SpectralDimension(**dimensions[i], events=events))

        if parse:
            method = {**prep, "spectral_dimensions": dim}
            return Method.parse_dict_with_units(method)
        else:
            return Method(**prep, spectral_dimensions=dim)

    @classmethod
    def parse_dict_with_units(cls, py_dict):
        """
        Parse the physical quantities of the method object from as a python dictionary.

        Args:
            py_dict: Dict object.
        """
        dict_copy = deepcopy(py_dict)
        dimensions = dict_copy["dimensions"]
        dict_copy.pop("dimensions")
        return cls.__new__(0, dimensions=dimensions, parse=True, **dict_copy)

    # creating class dynamically
    method = type(
        template["name"],
        (object,),
        {
            # constructor
            "__new__": constructor,
            "__str__": template["description"],
            # "__repr__": __repr__,
            "__doc__": __doc__,
            # data members
            # member functions
            "parse_dict_with_units": parse_dict_with_units,
            # "class_func": classMethod,
        },
    )
    return method


Bloch_decay = {
    "name": "BlochDecaySpectrum",
    "description": "A Bloch decay Spectrum.",
    "global_event_attributes": {},
    "spectral_dimensions": [
        {"events": [{"transition_query": {"P": {"channel-1": [[-1]]}}}]}
    ],
}

Bloch_decay_central_transition = {
    "name": "BlochDecaySpectrum",
    "description": "A central transition selective Bloch decay Spectrum.",
    "isotopomer_query": {
        "site_count": None,  # applies to any number of sites
        "site_isotope_is_quadrupole": None,  # applies to all isotopes
    },
    "global_event_attributes": {},
    "spectral_dimensions": [
        {
            "events": [
                {
                    "transition_query": {
                        "P": {"channel-1": [[-1]]},
                        "D": {"channel-1": [[0]]},
                    }
                }
            ]
        }
    ],
}

Bloch_decay_spectrum_doc = r"""
Method for simulating Bloch decay spectrum.

Args:
    dimensions: A list of python dict. Each dict is contains keywords that describe
        the coordinates along a spectral dimension. The keywords along with its
        definition are:

        count:
            An optional integer with the number of points, :math:`N`,
            along the dimension. The default value is ``1024``.
        spectral_width:
            A `required` float with the spectral width,
            :math:`\Delta x`, along the dimension in units of Hz.
        reference_offset:
            An `optional` float with the reference offset, :math:`x_0`
            along the dimension in units of Hz. The default value is ``0``.
        origin_offset:
            An `optional` float with the origin offset (Larmor frequency)
            along the dimension in units of Hz. The default value is ``0``.

    channel: Isotope on which the method will be applied.
    rotor_frequency: An `optional` float containing the sample spinning frequency
        :math:`\nu_r`, in units of Hz. The default value is ``0``.
    rotor_angle: An `optional` float containing the angle between the sample rotation
        axis and the applied external magnetic field, :math:`\theta`, in units of rad.
        The default value is ``0.9553166``, i.e. the magic angle.
    magetic_flux_density: An `optional` float containing the macroscopic magnetic
        flux density, :math:`H_0`, of the applied external magnetic field in units of
        T. The default value is ``9.4``.
"""
BlochDecaySpectrum = generate_method_from_template(
    Bloch_decay, Bloch_decay_spectrum_doc
)


Bloch_decay_ct_spectrum_doc = """
Method for simulating central transition selective Bloch decay spectrum.

Args:
    dimensions: A list of Dimension objects or equivalent python dict.
    channel: Isotope on which the method will be applied.
    rotor_frequency: An `optional` float containing the sample spinning frequency
        :math:`\nu_r`, in units of Hz. The default value is ``0``.
    rotor_angle: An `optional` float containing the angle between the sample rotation
        axis and the applied external magnetic field, :math:`\theta`, in units of rad.
        The default value is ``0.9553166``, i.e. the magic angle.
    magetic_flux_density: An `optional` float containing the macroscopic magnetic
        flux density, :math:`H_0`, of the applied external magnetic field in units of
        T. The default value is ``9.4``.
"""
BlochDecayCentralTransitionSpectrum = generate_method_from_template(
    Bloch_decay_central_transition, Bloch_decay_ct_spectrum_doc
)
