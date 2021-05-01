# -*- coding: utf-8 -*-
# import csdmpy as cp
from typing import Any
from typing import List

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


class ChemicalShiftSchema(Base):
    isotropic: str
    zeta: str = None
    eta: float = None
    alpha: str = None
    beta: str = None
    gamma: str = None


class QuadrupolarSchema(Base):
    Cq: str
    eta: float
    alpha: str = None
    beta: str = None
    gamma: str = None


class SiteSchema(Base):
    isotope: str
    ChemicalShift: ChemicalShiftSchema = None
    Quadrupolar: QuadrupolarSchema = None


class MethodSchema(Base):
    larmorFrequency: str
    spinningFrequency: str
    spectralWidth: str
    rotorAngle: str


class SimulatorSchema(Base):
    site: SiteSchema = None
    method: MethodSchema = None

    class config:
        validate_assignment = True
        arbitrary_types_allowed = True


class ContribSchema(Base):
    project: str
    identifier: str
    formula: str = None
    is_public: bool = None
    data: SimulatorSchema
    structures: list = None
    tables: list = None
    notebook: dict = None
    attachments: List[Any] = None
