"""Qchem Output File Parser."""
from __future__ import annotations

from operator import itemgetter
from re import match

import numpy as np
from mrsimulator import Coupling
from mrsimulator import Site
from mrsimulator import SpinSystem
from mrsimulator.spin_system.isotope import ISOTOPE_DATA
from mrsimulator.spin_system.tensors import SymmetricTensor
from pymatgen.analysis.nmr import ChemicalShielding
from pymatgen.io.qchem import check_for_structure_changes
from pymatgen.io.qchem.outputs import QCOutput
from pymatgen.io.qchem.utils import read_pattern
from pymatgen.io.qchem.utils import read_table_pattern


__author__ = "Maxwell C. Venetos"
__email__ = "mvenetos@berkeley.edu"


class QChem_NMR(QCOutput):
    """
    Class to parse QChem output files.
    """

    def __init__(self, filename: str):
        """
        Args:
            filename (str): Filename to parse
        """
        super().__init__(filename)

        self.data["nmr"] = read_pattern(
            self.text, {"key": r"(?i)\s*job(?:_)*type\s*(?:=)*\s*[(?i)nmr]"}
        ).get("key")
        if self.data.get("nmr", []):
            self._read_nmr_data()

        self.data["issc"] = read_pattern(
            self.text, {"key": r"(?i)\s*job(?:_)*type\s*(?:=)*\s*[(?i)issc]"}
        ).get("key")
        if self.data.get("issc", []):
            self._read_issc_data()

    def _read_nmr_data(self):
        header_pattern = r"\s*total shielding tensor\s*\n \
        \s*Trace\s*=\s*[0-9]*\.[0-9]+\n\s*Full Tensor:"
        table_pattern = (
            r"\s*(-?[0-9]*\.[0-9]*)\s*(-?[0-9]*\.[0-9]*)\s*(-?[0-9]*\.[0-9]*)\s*"
        )
        footer_pattern = r"\s*\-+"
        temp_tensors = read_table_pattern(
            self.text, header_pattern, table_pattern, footer_pattern, postprocess=float
        )
        temp_species = read_pattern(
            self.text, {"key": r"\s*\-*\s*ATOM\s*([a-zA-Z]+)\s*([0-9]+)\s*\-*\s*"}
        ).get("key")

        if temp_tensors is None or len(temp_tensors) == 0:
            self.data["nmr"] = None
        else:
            d = []
            for species, tensor in zip(temp_species, temp_tensors):
                CS_tensor = ChemicalShielding(tensor[0])
                Haeberlen_values = CS_tensor.haeberlen_values
                temp_dict = {
                    "species": species[0],
                    "isotropic_chemical_shift": Haeberlen_values[0],
                    "zeta": Haeberlen_values[2],
                    "eta": Haeberlen_values[3],
                }
                d.append(temp_dict)
            self.data["nmr"] = d

    def _read_issc_data(self):
        header_pattern = (
            r"\s*Total Spin-Spin Coupling Tensor J \(Hz\):\s*\n\s*1\s*2\s*3\s*"
        )
        table_pattern = r"\s*([1-3])\s*(-?[0-9]*\.[0-9]*) \
        \s*(-?[0-9]*\.[0-9]*)\s*(-?[0-9]*\.[0-9]*)\s*"
        footer_pattern = r"\s*Isotropic:\s*(-?[0-9]*\.[0-9]*)\s*"
        temp_tensors = read_table_pattern(
            self.text, header_pattern, table_pattern, footer_pattern, postprocess=float
        )
        temp_coupling = read_pattern(
            self.text,
            {
                "key": r"\s*Atoms\s*([a-zA-Z]+)\s*\(\#([0-9]+)\)\s*and \
                \s*([a-zA-Z]+)\s*\(\#([0-9]+)\)\s*with g-factors\s \
                *(-?[0-9]*\.[0-9]*)\s*and\s*(-?[0-9]*\.[0-9]*)\s*"
            },
        ).get("key")
        if temp_tensors is None or len(temp_tensors) == 0:
            self.data["issc"] = None
        else:
            d = []
            for coupling_pair, tensor in zip(temp_coupling, temp_tensors):
                tensor = np.asarray(tensor)
                temp_dict = {
                    "coupling": [coupling_pair[1], coupling_pair[3]],
                    "j_tensor": tensor[:, 1:],
                }
                d.append(temp_dict)
            self.data["issc"] = d

    def _get_sites(self, ignore_species=None, isotope_map={}):
        """
        ingore_species: list of species to ignore
        isotope map: dictionary with preferred isotope mappings ex: {'Na' : '29Na'}.
        Default behaviour is assign most abundant isotope
        """
        site_collection = []
        index_map = []
        for idx, base_atom in enumerate(self.data["species"]):

            if base_atom not in ignore_species:
                if base_atom in isotope_map.keys():
                    isotope_symbol = isotope_map[base_atom]
                else:
                    atom_types = []
                    for isotope in ISOTOPE_DATA.keys():
                        species = match(r"(\d+)\s*(\w+)", isotope)
                        if species.group(2).casefold() == base_atom.casefold():
                            atom_types.append(
                                [isotope, ISOTOPE_DATA[isotope]["natural_abundance"]]
                            )

                    isotope_symbol = sorted(
                        atom_types, key=itemgetter(1), reverse=True
                    )[0][0]

                temp_site = Site(
                    isotope=isotope_symbol,
                    isotropic_chemical_shift=self.data["nmr"][idx][
                        "isotropic_chemical_shift"
                    ],
                    shielding_symmetric=SymmetricTensor(
                        zeta=self.data["nmr"][idx]["zeta"],
                        eta=self.data["nmr"][idx]["eta"],
                    ),
                )
                index_map.append([idx, len(site_collection)])
                site_collection.append(temp_site)

        return site_collection, np.asarray(index_map)

    def _get_couplings(self, index_map):
        coupling_list = []

        index_list = index_map[:, 0]
        for coupled_pair in self.data["issc"]:
            if (
                coupled_pair["coupling"][0] in index_list
                and coupled_pair["coupling"][1] in index_list
            ):
                index_0 = index_list.where(coupled_pair["coupling"][0])
                index_1 = index_list.where(coupled_pair["coupling"][1])
                isotropic_j = (
                    coupled_pair["j_tensor"][0][0]
                    + coupled_pair["j_tensor"][1][1]
                    + coupled_pair["j_tensor"][2][2]
                ) / 3
                coupling_list.append(
                    Coupling(
                        site_index=[index_map[index_0, 1], index_map[index_1, 1]],
                        isotropic_j=isotropic_j,
                    )
                )
        return coupling_list

    def as_SpinSystem(self, ignore_species=None, isotope_map={}):
        """
        ingore_species: list of species to ignore
        isotope map: dictionary with preferred isotope mappings ex: {'Na' : '29Na'}.
         Default behaviour is assign most abundant isotope
        """
        site_collection, index_map = self._get_sites(ingore_species, isotope_map)

        if self.data["issc"] is not None:
            coupling = self._get_couplings(index_map)
        else:
            coupling = None

        return SpinSystem(sites=site_collection, couplings=coupling)

    def marge_data(self, QChem_to_SpinSystem_):
        check = check_for_structure_changes(
            self.data["initial_molecule"],
            QChem_to_SpinSystem_["data"]["initial_molecule"],
        )
        if check != "no change":
            raise ValueError("ERROR: Molecules are not consistent between files")

        for data_check in ["nmr", "issc"]:
            if self[data_check] is None:
                QChem_to_SpinSystem_["data"][data_check]
