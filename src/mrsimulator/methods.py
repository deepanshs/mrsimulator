# -*- coding: utf-8 -*-
from .method import Event
from .method import Method
from .method import Sequence


# program to create class dynamically


def generate_method_from_template(template):

    # constructor
    def constructor(self, dimensions, **kwargs):
        global_events = template["global_event_attributes"]
        ge = set(global_events)
        kw = set(kwargs)
        common = kw.intersection(ge)
        name = template["name"]
        if common != set():
            raise AttributeError(
                f"The attribute(s) {common} is fixed for {name} method."
            )

        seq = []
        for i, s in enumerate(template["sequences"]):
            events = []
            for e in s["events"]:
                ew = set(e)
                intersection = ew.intersection(kw)
                kw = kwargs.copy()
                for item in intersection:
                    kw.pop(item)
                events.append(Event(**e, **kw, **global_events))
            seq.append(Sequence(**dimensions[i], events=events))

        m = Method(
            name=name,
            isotope=kwargs["isotope"],
            sequences=seq,
            description=template["description"],
        )
        return m

    # class method
    @classmethod
    def classMethod(cls):
        return template["description"]

    # creating class dynamically
    method = type(
        template["name"],
        (object,),
        {
            # constructor
            "__new__": constructor,
            "__str__": template["description"],
            #         "__repr__": __repr__,
            "__doc__": template["description"],
            # data members
            # member functions
            "class_func": classMethod,
        },
    )
    return method


Bloch_decay = {
    "name": "BlochDecayFT",
    "description": "Bloch decay Fourier Transform",
    "isotopomer_query": {
        "site_count": None,  # applies to any number of sites
        "site_isotope_is_quadrupole": None,  # applies to all isotopes
    },
    "global_event_attributes": {},
    "sequences": [{"events": [{"transition_query": {"P": [-1]}}]}],
}


def BlochDecayFT(dimensions, isotope, **kwargs):
    """A Bloch decay Fourier Transform method.

    Args:
        dimensions: A list of Dimension objects or equivalent python dict.
        isotope: The isotope on which the method will be applied.
        rotor_frequency: The sample rotation frequency in Hz.
        rotor_angle: The sample holder angle (inn radians) with respect to lab
                frame z-axis.
        magetic_flux_density: The magnetic flux density of macroscopic external
                magetic field.
    """
    return generate_method_from_template(Bloch_decay)(
        dimensions=dimensions, isotope=isotope, **kwargs
    )
