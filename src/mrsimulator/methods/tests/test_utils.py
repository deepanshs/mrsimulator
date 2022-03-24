# -*- coding: utf-8 -*-
from copy import deepcopy

import pytest
from mrsimulator.utils.utils import map_transition_query_object_to_v_7


def template(tq):
    return {"spectral_dimensions": [{"events": [{"transition_query": tq}]}]}


def test_map_transition_query_object_to_v_7():
    tq = {"P": -1, "D": -2}
    mth = template(tq)
    res = deepcopy(mth)
    map_transition_query_object_to_v_7(res)
    assert res == template([{"P": [-1], "D": [-2]}])

    tq = {"P": -1}
    mth = template(tq)
    res = deepcopy(mth)
    map_transition_query_object_to_v_7(res)
    assert res == template([{"P": [-1]}])

    tq = {"D": 2}
    mth = template(tq)
    res = deepcopy(mth)
    map_transition_query_object_to_v_7(res)
    assert res == template([{"D": [2]}])

    tq = {"P": -1, "D": [-2, -1]}
    mth = template(tq)
    res = deepcopy(mth)
    error = "Ambiguous definition for transition queries"
    with pytest.raises(Exception, match=f".*{error}.*"):
        map_transition_query_object_to_v_7(res)
    # assert res == template([{"P": [-1], "D": [-2]}, {"P": [-1], "D": [-1]}])

    # tq = {"D": [-2, -1]}
    # mth = template(tq)
    # res = deepcopy(mth)
    # map_transition_query_object_to_v_7(res)
    # assert res == template([{"D": [-2]}, {"D": [-1]}])

    # tq = {"P": [-2], "D": [-4, 4]}
    # mth = template(tq)
    # res = deepcopy(mth)
    # map_transition_query_object_to_v_7(res)
    # assert res == template([{"P": [-2], "D": [-4]}, {"P": [-2], "D": [4]}])

    # tq = {"P": [-2, 2], "D": [-4, 4]}
    # mth = template(tq)
    # res = deepcopy(mth)
    # map_transition_query_object_to_v_7(res)
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
    # map_transition_query_object_to_v_7(res)
    # assert res == template([{"P": [1, 1], "D": [-4]}, {"P": [1, 1], "D": [4]}])

    tq = {"P": [[1, 1]]}
    mth = template(tq)
    res = deepcopy(mth)
    map_transition_query_object_to_v_7(res)
    assert res == template([{"P": [1, 1]}])

    mth = {"spectral_dimensions": [{}]}
    res = deepcopy(mth)
    map_transition_query_object_to_v_7(res)
    assert res == mth

    mth = {"spectral_dimensions": [{"events": [{}]}]}
    res = deepcopy(mth)
    map_transition_query_object_to_v_7(res)
    assert res == mth

    tq = {"ch1": {"P": -1, "D": -2}}
    mth = template(tq)
    res = deepcopy(mth)
    map_transition_query_object_to_v_7(res)
    assert res == template([{"ch1": {"P": -1, "D": -2}}])
