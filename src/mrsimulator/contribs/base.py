# -*- coding: utf-8 -*-
# import csdmpy as cp
from pathlib import Path
from typing import List
from typing import Union

from mpcontribs.client import Attachment
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
    """Schema for Chemical shift parameters.

    Args:
        str isotropic: Isotropic chemical shift.
        str zeta: Chemical shift anisotropy.
        float eta: Shift asymmetry.
        str alpha: Euler angle alpha.
        str beta: Euler angle beta.
        str gamma: Euler angle gamma.
    """

    isotropic: str
    zeta: str = None
    eta: float = None
    alpha: str = None
    beta: str = None
    gamma: str = None


class QuadrupolarSchema(Base):
    """Schema for Quadrupolar parameters.

    Args:
        str Cq: Quadrupolar coupling constant.
        float eta: Quadrupolar asymmetry.
        str alpha: Euler angle alpha.
        str beta: Euler angle beta.
        str gamma: Euler angle gamma.
    """

    Cq: str
    eta: float
    alpha: str = None
    beta: str = None
    gamma: str = None


class SiteSchema(Base):
    """Schema for Site parameters.

    Args:
        str isotope: The isotope given as atomic number followed by atomic symbol.
        ChemicalShiftSchema ChemicalShift: The chemical shift scheme.
        QuadrupolarSchema Quadrupolar: The quadrupolar scheme.
    """

    isotope: str
    ChemicalShift: ChemicalShiftSchema = None
    Quadrupolar: QuadrupolarSchema = None


class MethodSchema(Base):
    """Schema for Method parameters.

    Args:
        str larmorFrequency: The larmor frequency of observed isotope.
        str spinningFrequency: The rotor spinning frequency.
        str spectralWidth: The spectral width of the spectrum.
        str rotorAngle: Angle of sample rotation axis wrt the z-axis.
    """

    larmorFrequency: str
    spinningFrequency: str
    spectralWidth: str
    rotorAngle: str


class SimulatorSchema(Base):
    """Schema for Simulator parameters.

    Args:
        SiteSchema site: The site object.
        MethodSchema method: The method object.
    """

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
    attachments: List[Union[Attachment, Path]] = None
