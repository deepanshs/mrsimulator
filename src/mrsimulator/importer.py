# -*- coding: utf-8 -*-
"""Utility function module."""
import json
import sys
from os import path
from urllib.parse import urlparse

import requests


__author__ = "Deepansh J. Srivastava"
__email__ = "deepansh2012@gmail.com"


def download_file_from_url(url):
    """Download the file from the given url.

    Args:
        str url: The url address.
    """
    res = urlparse(url)
    filename = path.split(res[2])[1]
    name, extension = path.splitext(filename)
    original_name = name
    i = 0
    while path.isfile(filename):
        filename = "{0}_{1}{2}".format(original_name, str(i), extension)
        i += 1

    with open(filename, "wb") as f:
        response = requests.get(url, stream=True)
        total = response.headers.get("content-length")

        if total is None:
            f.write(response.content)
        else:
            downloaded = 0
            total = int(total)
            sys.stdout.write(
                "Downloading '{0}'\nfrom '{1}' to file '{2}'.".format(
                    res[2], res[1], filename
                )
            )
            for data in response.iter_content(
                chunk_size=max(int(total / 1000), 1024 * 1024)
            ):
                downloaded += len(data)
                f.write(data)
                done = int(4 * downloaded / total)
                sys.stdout.write("\n[{}{}]".format("█" * done, "." * (4 - done)))
                sys.stdout.flush()

    sys.stdout.write("\n")

    return filename


def import_json(filename):
    res = urlparse(filename)
    if res[0] not in ["file", ""]:
        filename = download_file_from_url(filename)
    with open(filename, "rb") as f:
        content = f.read()
        return json.loads(str(content, encoding="UTF-8"))
