# -*- coding: utf-8 -*-

import os
import numpy
from setuptools import setup, find_packages
from setuptools.extension import Extension
from Cython.Build import cythonize

module_dir = os.path.dirname(os.path.abspath(__file__))

ext_modules = [
    Extension(
        name="mrsimulator.methods",
        sources=[
            "mrsimulator/scr/lib/c_array.c",
            "mrsimulator/scr/lib/MRAngularMomentum.c",
            "mrsimulator/scr/mrmethods/spinning_sidebands.c",
            "mrsimulator/scr/mrmethods/powder_setup.c",
            "mrsimulator/scr/mrmethods/nmr_methods.pyx",
        ],
        include_dirs=["mrsimulator/scr/include", numpy.get_include()],
        libraries=["fftw3"],
        language="c",
        extra_compile_args="-O1".split(),
        extra_link_args="-g -lfftw3 -lmkl_intel_lp64 -lmkl_intel_thread \
                        -lmkl_core -ldl -liomp5 -lm -W".split(),
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
    packages=find_packages(),
    package_data={},
    install_requires=["numpy>=1.13.3", "astropy>=3.0", "pydantic", "requests>=2.21.0", "monty==2.0.4"],
    extras_require={"fancy feature": [
        "matplotlib>=3.0.2",
        "plotly>=3.6",
        "dash>=0.40",
        "dash_daq>=0.1",
    ]},
    tests_require=["nose"],
    entry_points={"console_scripts": ["nmr_app = mrsimulator.web_interface:main"]},
    ext_modules=cythonize(ext_modules, annotate=True, language_level=3, gdb_debug=True),
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
    ],
)
