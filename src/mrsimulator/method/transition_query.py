# -*- coding: utf-8 -*-
from typing import Dict
from typing import List
from typing import Optional

from mrsimulator.transition import Transition
from pydantic import BaseModel
from pydantic import Field

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


class TransitionQuery(BaseModel):
    """Base TransitionQuery class.

    Attributes
    ----------

    P:
        A dict of channels, with each channel as a list of p symmetry functions per
        site. Here p = Δm is the difference between spin quantum numbers of the final
        and initial states.

        Example
        -------

        >>> method = Method2D()
        >>> method.spectral_dimensions[0].events[0].transition_query.P = {
        ...     'channel-1': [[-1]]
        ... }

    D:
        A dict of channels, with each channel as a list of d symmetry functions per
        site. Here p = Δm is the difference between spin quantum numbers of the final
        and initial states.

        Example
        -------

        >>> method.spectral_dimensions[0].events[0].transition_query.D = {
        ...     'channel-1': [[0]]
        ... }
    """

    P: Optional[Dict] = {"channel-1": [[-1.0]]}
    D: Optional[Dict] = Field(default=None)
    F: Optional[Dict] = Field(default=None)
    transitions: List[Transition] = None

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True

    def json(self) -> dict:
        """Parse the class object to a JSON compliant python dictionary object, where
        the attribute value with physical quantity is expressed as a string with a
        value and a unit."""
        return {k: v for k, v in self.dict().items() if v is not None}
