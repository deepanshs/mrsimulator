# -*- coding: utf-8 -*-
# import csdmpy as cp
from pydantic import BaseModel

__author__ = "Deepansh Srivastava"
__email__ = "srivastava.89@osu.edu"


class Base(BaseModel):
    def json(self):
        return Base.fullsimplify(super().dict())

    @staticmethod
    def simplify(val):
        """Remove value if it is None."""
        return {k: v for k, v in val.items() if v is not None}

    @staticmethod
    def fullsimplify(val):
        """Iteratively remove None values from a nested dict."""
        initial = {
            k: Base.simplify(Base.fullsimplify(v)) if isinstance(v, dict) else v
            for k, v in val.items()
        }
        return Base.simplify(initial)


class ChemicalShiftInfo(Base):
    isotropic: str
    zeta: str = None
    eta: float = None
    alpha: str = None
    beta: str = None
    gamma: str = None


class QuadrupolarInfo(Base):
    Cq: str
    eta: float
    alpha: str = None
    beta: str = None
    gamma: str = None


class SiteInfo(Base):
    isotope: str
    ChemicalShift: ChemicalShiftInfo = None
    Quadrupolar: QuadrupolarInfo = None


class MethodInfo(Base):
    larmorFrequency: str
    spinningFrequency: str
    spectralWidth: str
    rotorAngle: str


class SimulatorInfo(Base):
    experiment: str = None
    simulation: str = None
    site: SiteInfo = None
    method: MethodInfo = None

    class config:
        validate_assignment = True
        arbitrary_types_allowed = True


class ContribInfo(Base):
    data: SimulatorInfo
    project: str
    identifier: str
    composition: str = None
