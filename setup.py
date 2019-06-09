# -*- coding: utf-8 -*-
import io
from os import listdir
from os import path
from os.path import abspath
from os.path import dirname
from os.path import join

from setuptools import Extension
from setuptools import find_packages
from setuptools import setup

# from mrsimulator.__version__ import __version__

# Package meta-data.
NAME = "mrsimulator"
DESCRIPTION = "A python toolbox for simulating NMR spectra"
URL = "https://github.com/DeepanshS/MRsimulator/"
EMAIL = "srivastava.89@osu.edu"
AUTHOR = "Deepansh J. Srivastava"
REQUIRES_PYTHON = ">=3.0"
VERSION = "0.1.1"


# What packages are required for this module to be executed?
REQUIRED = [
    "numpy>=1.13.3",
    "astropy>=3.0",
    "mkl",
    "mkl-include",
    "pytest-runner>=5.0",
    "pytest",
]

# What packages are optional?
EXTRAS = {
    "fancy feature": [
        "matplotlib>=3.0.2",
        "plotly>=3.6",
        "dash>=0.40",
        "dash_daq>=0.1",
    ]
}


here = abspath(dirname(__file__))
# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(join(here, "README.md"), encoding="utf-8") as f:
        long_description = "\n" + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(join(here, project_slug, "__version__.py")) as f:
        exec(f.read(), about)
else:
    about["__version__"] = VERSION


# Cython =========================================================== #
cmdclass = {}
ext_modules = []


# scr folder
_list = ["mrsimulator", "scr", "lib"]
nmr_lib_source_dir = path.join(*_list)
_source_files = []
for _file in listdir(nmr_lib_source_dir):
    if _file.endswith(".c"):
        _source_files.append(path.join(nmr_lib_source_dir, _file))


_list = ["mrsimulator", "scr", "include"]
include_nmr_lib_directories = [path.join(*_list)]


def get_include_and_lib_paths():
    """Get paths to include and lib folders on windows, linux and mac os."""
    import numpy

    # numpy include
    path_ = numpy.get_include()
    print("numpy include:", path_)
    include_nmr_lib_directories.append(path_)

    # conda lib
    for _ in range(5):
        path_ = path.split(path_)[0]
    path_lib = path_
    print("conda lib:")
    print("exist:", path.exists(path_lib))
    print("path:", path_lib)
    if path.exists(path_lib):
        include_nmr_lib_directories.append(path_lib)

    # conda include
    path_include_conda = path.join(path.split(path_)[0], "include")
    print("conda include:")
    print("exist:", path.exists(path_include_conda))
    print("path:", path_include_conda)
    if path.exists(path_include_conda):
        include_nmr_lib_directories.append(path_include_conda)

    # other include
    path_include_other = path.join(path_lib, "Library", "include")
    print("other include:")
    print("exist:", path.exists(path_include_other))
    print("path:", path_include_other)
    if path.exists(path_include_other):
        include_nmr_lib_directories.append(path_include_other)


get_include_and_lib_paths()

nmr_function_source_file = _source_files[:]

_list = ["mrsimulator", "scr", "mrmethods"]
nmr_function_source_dir = path.join(*_list)

for _file in listdir(nmr_function_source_dir):
    if _file.endswith(".c"):
        nmr_function_source_file.append(
            path.join(nmr_function_source_dir, _file)
        )


ext_modules = [
    Extension(
        name=NAME + ".methods",
        sources=nmr_function_source_file,
        include_dirs=include_nmr_lib_directories,
        language="c",
        extra_compile_args="-O1".split(),
        extra_link_args="-g -lfftw3 -lmkl_intel_lp64 -lmkl_intel_thread \
                        -lmkl_core -ldl -liomp5 -lm -W".split(),
    )
]

ext = ext_modules

setup(
    name=NAME,
    version=about["__version__"],
    description=DESCRIPTION,
    long_description=long_description,
    author=AUTHOR,
    author_email=EMAIL,
    python_requires=REQUIRES_PYTHON,
    url=URL,
    packages=find_packages(),
    include_package_data=True,
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    setup_requires=["pytest-runner", "numpy>=1.13.3", "mkl", "mkl-include"],
    tests_require=["pytest"],
    package_data={"mrsimulator": ["tests/*.*"]},
    cmdclass=cmdclass,
    ext_modules=ext,
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
    ],
)
