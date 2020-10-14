# -*- coding: utf-8 -*-
"""Utility function module."""
import json
from urllib.parse import urlparse

from csdmpy.dependent_variables.download import download_file_from_url


__author__ = "Deepansh J. Srivastava"
__email__ = "srivastava.89@osu.edu"


def import_json(filename):
    res = urlparse(filename)
    if res[0] not in ["file", ""]:
        filename = download_file_from_url(filename)
    with open(filename, "rb") as f:
        content = f.read()
        return json.loads(str(content, encoding="UTF-8"))
