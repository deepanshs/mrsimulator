# -*- coding: utf-8 -*-
from copy import deepcopy

import pytest
from mrsimulator.utils.utils import convert_transition_query


def template(tq):
    return {"spectral_dimensions": [{"events": [{"transition_query": tq}]}]}


def test_convert_transition_query():
    tq = {"P": -1, "D": -2}
    mth = template(tq)
    res = deepcopy(mth)
    convert_transition_query(res)
    assert res == template([{"ch1": {"P": [-1], "D": [-2]}}])

    tq = {"P": -1}
    mth = template(tq)
    res = deepcopy(mth)
    convert_transition_query(res)
    assert res == template([{"ch1": {"P": [-1]}}])

    tq = {"D": 2}
    mth = template(tq)
    res = deepcopy(mth)
    convert_transition_query(res)
    assert res == template([{"ch1": {"D": [2]}}])

    tq = {"P": -1, "D": [-2, -1]}
    mth = template(tq)
    res = deepcopy(mth)
    error = "Ambiguous definition for transition queries"
    with pytest.raises(Exception, match=f".*{error}.*"):
        convert_transition_query(res)

    # assert res == template([{"P": [-1], "D": [-2]}, {"P": [-1], "D": [-1]}])
    # tq = {"D": [-2, -1]}
    # mth = template(tq)
    # res = deepcopy(mth)
    # convert_transition_query(res)
    # assert res == template([{"D": [-2]}, {"D": [-1]}])
    # tq = {"P": [-2], "D": [-4, 4]}
    # mth = template(tq)
    # res = deepcopy(mth)
    # convert_transition_query(res)
    # assert res == template([{"P": [-2], "D": [-4]}, {"P": [-2], "D": [4]}])
    # tq = {"P": [-2, 2], "D": [-4, 4]}
    # mth = template(tq)
    # res = deepcopy(mth)
    # convert_transition_query(res)
    # assert res == template(
    #     [
    #         {"P": [-2], "D": [-4]},
    #         {"P": [-2], "D": [4]},
    #         {"P": [2], "D": [-4]},
    #         {"P": [2], "D": [4]},
    #     ]
    # )
    # tq = {"P": [[1, 1]], "D": [-4, 4]}
    # mth = template(tq)
    # res = deepcopy(mth)
    # convert_transition_query(res)
    # assert res == template([{"P": [1, 1], "D": [-4]}, {"P": [1, 1], "D": [4]}])

    tq = {"P": [[1, 1]]}
    mth = template(tq)
    res = deepcopy(mth)
    convert_transition_query(res)
    assert res == template([{"ch1": {"P": [1, 1]}}])

    mth = {"spectral_dimensions": [{}]}
    res = deepcopy(mth)
    convert_transition_query(res)
    assert res == mth

    mth = {"spectral_dimensions": [{"events": [{}]}]}
    res = deepcopy(mth)
    convert_transition_query(res)
    assert res == mth
