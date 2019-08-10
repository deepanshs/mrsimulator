# -*- coding: utf-8 -*-

import os
import numpy
import numpy.distutils.system_info as sysinfo
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

module_dir = os.path.dirname(os.path.abspath(__file__))

mkl_info = sysinfo.get_info("mkl")

include_dirs = ["src/c_lib/include", numpy.get_include()]

if "include_dirs" in mkl_info:
    include_dirs += mkl_info["include_dirs"]

extra_compile_args = "-O1 -std=c11"

ext_modules = [
    Extension(
        name="mrsimulator.methods",
        sources=[
            "src/c_lib/lib/c_array.c",
            "src/c_lib/lib/MRAngularMomentum.c",
            "src/c_lib/mrmethods/spinning_sidebands.c",
            "src/c_lib/mrmethods/powder_setup.c",
            "src/c_lib/mrmethods/nmr_methods.pyx",
        ],
        include_dirs=include_dirs,
        language="c",
        libraries=["fftw3","mkl_intel_lp64","mkl_core","mkl_intel_thread","pthread","iomp5","m","dl"],
        extra_compile_args=extra_compile_args.split(),
    )
]

setup(
    name="mrsimulator",
    version="0.1.0",
    description="A python toolbox for simulating NMR spectra",
    long_description=open(os.path.join(module_dir, "README.md")).read(),
    author="Deepansh J. Srivastava",
    author_email="srivastava.89@osu.edu",
    python_requires=">=3.0",
    url="https://github.com/DeepanshS/MRsimulator/",
    packages=find_packages("src"),
    package_dir={"": "src"},
    install_requires=[
        "numpy>=1.13.3",
        "astropy>=3.0",
        "pydantic==0.28",
        "requests>=2.21.0",
        "monty==2.0.4",
        "mkl==2019.0",
        "mkl-include==2019.0",
    ],
    tests_require=["pytest"],
    ext_modules=cythonize(ext_modules),
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
    ],
)
