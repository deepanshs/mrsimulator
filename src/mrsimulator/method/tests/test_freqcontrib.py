from mrsimulator.method.frequency_contrib import FrequencyEnum

__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


def test_01():
    a = FrequencyEnum("Quad2_4")
    assert a.value == "Quad2_4"
