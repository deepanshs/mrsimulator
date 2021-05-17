# -*- coding: utf-8 -*-
import platform
import sys
from os.path import abspath
from os.path import dirname
from os.path import join
from os.path import split

from setuptools import Extension
from setuptools import find_packages
from setuptools import setup

import numpy as np

try:
    from Cython.Build import cythonize

    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

# get the version from file
python_version = sys.version_info
py_version = ".".join([str(i) for i in python_version[:3]])
print("Using python version", py_version)
if python_version.major != 3 and python_version.minor < 6:
    print(
        f"Python version 3.6 and higher is required for the setup. You are using "
        f"version {py_version}"
    )
    sys.exit(1)

with open("src/mrsimulator/__init__.py", "r") as f:
    for line in f.readlines():
        if "__version__" in line:
            before_keyword, keyword, after_keyword = line.partition("=")
            version = after_keyword.strip()[1:-1]
            print("mrsimulator version ", version)

module_dir = dirname(abspath(__file__))


libraries = []
include_dirs = []
library_dirs = []
extra_compile_args = []
extra_link_args = []
data_files = []

numpy_include = np.get_include()

conda_location = numpy_include
for _ in range(5):
    conda_location = split(conda_location)[0]

if platform.system() == "Windows":
    # windows system lib and include path
    include_dirs += [join(conda_location, "Library", "include", "fftw")]
    include_dirs += [join(conda_location, "Library", "include", "openblas")]
    include_dirs += [join(conda_location, "Library", "include")]
    include_dirs += [join(conda_location, "include")]
    library_dirs += [join(conda_location, "Library", "lib")]
    extra_compile_args += ["-DFFTW_DLL", "/DUSE_OPENBLAS"]

else:
    # unix system lib and include path
    conda_location = split(conda_location)[0]
    include_dirs += [join(conda_location, "include")]
    library_dirs += [join(conda_location, "lib")]
    extra_compile_args = ["-O3", "-ffast-math", "-DUSE_OPENBLAS"]

libraries += ["fftw3", "openblas"]
extra_link_args += ["-lm"]

include_dirs = list(set(include_dirs))
library_dirs = list(set(library_dirs))
libraries = list(set(libraries))

# other include paths
include_dirs += ["src/c_lib/include/", numpy_include]

# print info
print(include_dirs)
print(library_dirs)
print(libraries)
print(extra_compile_args)
print(extra_link_args)

source = [
    "src/c_lib/lib/angular_momentum.c",
    "src/c_lib/lib/interpolation.c",
    "src/c_lib/lib/mrsimulator.c",
    "src/c_lib/lib/octahedron.c",
    "src/c_lib/lib/simulation.c",
    "src/c_lib/lib/frequency_averaging.c",
    "src/c_lib/lib/schemes.c",
    "src/c_lib/lib/method.c",
]

ext = ".pyx" if USE_CYTHON else ".c"

# method
ext_modules = [
    Extension(
        name="mrsimulator.base_model",
        sources=[*source, "src/c_lib/base/base_model" + ext],
        include_dirs=include_dirs,
        language="c",
        libraries=libraries,
        library_dirs=library_dirs,
        # data_files=data_files,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )
]

# tests
ext_modules += [
    Extension(
        name="mrsimulator.tests.tests",
        sources=[*source, "src/c_lib/test/test" + ext],
        include_dirs=include_dirs,
        language="c",
        libraries=libraries,
        library_dirs=library_dirs,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )
]

# sandbox
ext_modules += [
    Extension(
        name="mrsimulator.sandbox",
        sources=[*source, "src/c_lib/sandbox/sandbox" + ext],
        include_dirs=include_dirs,
        language="c",
        libraries=libraries,
        library_dirs=library_dirs,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )
]

if USE_CYTHON:
    ext_modules = cythonize(ext_modules, language_level=3)

extras = {"lmfit": ["lmfit>=1.0.2"]}

description = "A python toolbox for simulating fast real-time solid-state NMR spectra."
setup(
    name="mrsimulator",
    version=version,
    description=description,
    long_description=open(join(module_dir, "README.md")).read(),
    author="Deepansh J. Srivastava",
    author_email="srivastava.89@osu.edu",
    python_requires=">=3.6",
    url="https://github.com/DeepanshS/MRsimulator/",
    packages=find_packages("src"),
    package_dir={"": "src"},
    extras_require=extras,
    ext_modules=ext_modules,
    include_package_data=True,
    zip_safe=False,
    license="BSD-3-Clause",
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "Intended Audience :: Science/Research",
        "Intended Audience :: Education",
        "Intended Audience :: Developers",
        "Operating System :: OS Independent",
        "Development Status :: 2 - Pre-Alpha",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: C",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Topic :: Education",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
