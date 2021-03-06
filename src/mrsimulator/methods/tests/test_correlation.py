# -*- coding: utf-8 -*-
# import pytest
# from mrsimulator.method.transition_query import TransitionQuery
# from mrsimulator.methods import Cosy
# from mrsimulator.methods import Inadequate
#
# __author__ = "Deepansh J. Srivastava"
# __email__ = "srivastava.89@osu.edu"
#
# methods = [Cosy, Inadequate]
# names = ["Cosy", "Inadequate"]
#
#
# def test_coorelation_rotor_freq():
#     def error(name):
#         return f"`rotor_frequency` value cannot be modified for {name} method."
#
#     for name, method in zip(names, methods):
#         e = error(name)
#         with pytest.raises(ValueError, match=f".*{e}.*"):
#             method(rotor_frequency=10, spectral_dimensions=[{}, {}])
#
#
# def test_coorelation_spectral_dimension_count():
#     e = "Method requires exactly 2 spectral dimensions, given 1."
#     for _, method in zip(names, methods):
#         with pytest.raises(ValueError, match=f".*{e}.*"):
#             method(spectral_dimensions=[{}])
#
#
# def test_coorelation_setting_transition_query():
#     def error(name):
#         return f"`transition_query` value cannot be modified for {name} method."
#
#     for name, method in zip(names, methods):
#         e = error(name)
#         with pytest.raises(ValueError, match=f".*{e}.*"):
#             method(
#                 spectral_dimensions=[
#                     {"events": [{"transition_query": {"P": [-1]}}]},
#                     {},
#                 ],
#             )
#
#
# def test_Cosy_general():
#     """Inner satellite-transition variable-angle spinning method"""
#     mth = Cosy(
#         channels=["13C"],
#         magnetic_flux_density=9.4,  # in T
#         spectral_dimensions=[
#             {
#                 "count": 1024,
#                 "spectral_width": 5e4,  # in Hz
#                 "reference_offset": 0,  # in Hz
#             },
#             {
#                 "count": 1024,
#                 "spectral_width": 5e4,  # in Hz
#                 "reference_offset": 0,  # in Hz
#             },
#         ],
#     )
#     assert mth.name == "Cosy"
#
#     des = "Simulate an infinite spinning COrrelation SpectroscopY spectrum."
#     assert mth.description == des
#     assert mth.spectral_dimensions[0].events[0].transition_query == TransitionQuery(
#         P={"channel-1": [[-1]]}
#     )
#     assert mth.spectral_dimensions[1].events[0].transition_query == TransitionQuery(
#         P={"channel-1": [[-1]]}
#     )
#     assert Cosy.parse_dict_with_units(mth.json()) == mth
#
#     serialize = mth.json()
#     assert serialize == {
#         "channels": ["13C"],
#         "description": des,
#         "magnetic_flux_density": "9.4 T",
#         "name": "Cosy",
#         "rotor_angle": "0.955316618 rad",
#         "rotor_frequency": "1000000000000.0 Hz",
#         "spectral_dimensions": [
#             {
#                 "count": 1024,
#                 "reference_offset": "0.0 Hz",
#                 "spectral_width": "50000.0 Hz",
#             },
#             {
#                 "count": 1024,
#                 "reference_offset": "0.0 Hz",
#                 "spectral_width": "50000.0 Hz",
#             },
#         ],
#     }
#
#
# def test_Inadequate_general():
#     """Second to inner satellite-transition variable-angle spinning method"""
#     mth = Inadequate(
#         channels=["17O"],
#         magnetic_flux_density=9.4,  # in T
#         spectral_dimensions=[
#             {
#                 "count": 1024,
#                 "spectral_width": 5e4,  # in Hz
#                 "reference_offset": 0,  # in Hz
#             },
#             {
#                 "count": 1024,
#                 "spectral_width": 5e4,  # in Hz
#                 "reference_offset": 0,  # in Hz
#             },
#         ],
#     )
#     assert mth.name == "Inadequate"
#
#     des = (
#         "Simulate an infinite spinning Incredible Natural Abundance DoublE QUAntum "
#         "Transfer Experiment spectrum."
#     )
#     assert mth.description == des
#     assert mth.spectral_dimensions[0].events[0].transition_query == TransitionQuery(
#         P={"channel-1": [[-1, -1]]}
#     )
#     assert mth.spectral_dimensions[1].events[0].transition_query == TransitionQuery(
#         P={"channel-1": [[-1]]}
#     )
#     assert Inadequate.parse_dict_with_units(mth.json()) == mth
#
#     serialize = mth.json()
#     assert serialize == {
#         "channels": ["17O"],
#         "description": des,
#         "magnetic_flux_density": "9.4 T",
#         "name": "Inadequate",
#         "rotor_angle": "0.955316618 rad",
#         "rotor_frequency": "1000000000000.0 Hz",
#         "spectral_dimensions": [
#             {
#                 "count": 1024,
#                 "reference_offset": "0.0 Hz",
#                 "spectral_width": "50000.0 Hz",
#             },
#             {
#                 "count": 1024,
#                 "reference_offset": "0.0 Hz",
#                 "spectral_width": "50000.0 Hz",
#             },
#         ],
#     }
