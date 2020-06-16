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

    Attributes:
        P: A list of p symmetry transition, where p = Δm and Δm is the difference
                between spin quantum numbers of the final and initial states.
    """

    P: Optional[Dict] = {"channel-1": [[-1.0]]}
    D: Optional[Dict] = Field(default=None)
    f: Optional[Dict] = Field(default=None)
    transitions: List[Transition] = None

    class Config:
        validate_assignment = True
        arbitrary_types_allowed = True

    def to_dict_with_units(self):
        return {k: v for k, v in self.dict().items() if v is not None}
