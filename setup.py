# -*- coding: utf-8 -*-

from setuptools import Extension
from setuptools import find_packages
from setuptools import setup

from Cython.Build import cythonize

from os.path import join, split
from os.path import abspath
from os.path import dirname
import platform

import numpy as np
import numpy.distutils.system_info as sysinfo

# get the version from file

with open("src/mrsimulator/__init__.py", "r") as f:
    for line in f.readlines():
        if "__version__" in line:
            before_keyword, keyword, after_keyword = line.partition("=")
            version = after_keyword.strip()[1:-1]
            print(version)

module_dir = dirname(abspath(__file__))
print(module_dir)

include_dirs = [
    "/usr/include/",
    "/usr/include/openblas",
    "/usr/include/x86_64-linux-gnu/",
]
include_dirs += [f"/app/.apt{item}" for item in include_dirs]

library_dirs = ["/usr/lib64/", "/usr/lib/", "/usr/lib/x86_64-linux-gnu/"]
library_dirs += [f"/app/.apt{item}" for item in library_dirs]

libraries = ["openblas", "fftw3", "fftw3_threads", "pthread"]
data_files = []

numpy_include = np.get_include()

if platform.system() == "Windows":
    conda_location = numpy_include
    for _ in range(5):
        conda_location = split(conda_location)[0]
    include_dirs += [join(conda_location, "Library", "include", "fftw")]
    include_dirs += [join(conda_location, "Library", "include", "openblas")]
    include_dirs += [join(conda_location, "Library", "include")]
    include_dirs += [join(conda_location, "include")]
    library_dirs += [join(conda_location, "Library", "lib")]
    libraries += ["fftw3", "openblas"]
    name = "openblas"

    extra_link_args = ["-lm"]
    extra_compile_args = ["-DFFTW_DLL"]

else:
    openblas_info = sysinfo.get_info("openblas")
    fftw3_info = sysinfo.get_info("fftw3")

    if openblas_info != {}:
        name = "openblas"
        library_dirs += openblas_info["library_dirs"]
        libraries += openblas_info["libraries"]
    fftw_keys = fftw3_info.keys()
    if "include_dirs" in fftw_keys:
        include_dirs += fftw3_info["include_dirs"]
    if "library_dirs" in fftw_keys:
        library_dirs += fftw3_info["library_dirs"]
    if "libraries" in fftw_keys:
        libraries += fftw3_info["libraries"]

    extra_link_args = ["-lm"]
    extra_compile_args = ["-g", "-O3"]
include_dirs = list(set(include_dirs))
library_dirs = list(set(library_dirs))
libraries = list(set(libraries))

# other include paths
include_dirs += [f"{module_dir}/src/c_lib/include/", numpy_include]
print(include_dirs)
print(library_dirs)
print(libraries)

print(extra_compile_args)
print(extra_link_args)

# method
ext_modules = [
    Extension(
        name="mrsimulator.methods",
        sources=[
            "src/c_lib/lib/angular_momentum.c",
            "src/c_lib/lib/interpolation.c",
            "src/c_lib/lib/mrsimulator.c",
            "src/c_lib/lib/octahedron.c",
            "src/c_lib/lib/spinning_sidebands.c",
            "src/c_lib/lib/powder_setup.c",
            "src/c_lib/lib/averaging_scheme.c",
            "src/c_lib/mrmethods/nmr_methods.pyx",
        ],
        include_dirs=include_dirs,
        language="c",
        libraries=libraries,
        library_dirs=library_dirs,
        data_files=data_files,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )
]

# tests

ext_modules += [
    Extension(
        name="mrsimulator.tests.tests",
        sources=[
            "src/c_lib/lib/angular_momentum.c",
            "src/c_lib/lib/interpolation.c",
            "src/c_lib/lib/mrsimulator.c",
            "src/c_lib/lib/octahedron.c",
            "src/c_lib/lib/spinning_sidebands.c",
            "src/c_lib/lib/powder_setup.c",
            "src/c_lib/lib/averaging_scheme.c",
            "src/c_lib/test/test.pyx",
        ],
        include_dirs=include_dirs,
        language="c",
        libraries=libraries,
        library_dirs=library_dirs,
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )
]

# # sandbox

# ext_modules += [
#     Extension(
#         name="mrsimulator.sandbox",
#         sources=[
#             "src/c_lib/lib/angular_momentum.c",
#             "src/c_lib/lib/interpolation.c",
#             "src/c_lib/lib/mrsimulator.c",
#             "src/c_lib/lib/octahedron.c",
#             "src/c_lib/lib/spinning_sidebands.c",
#             "src/c_lib/lib/powder_setup.c",
#             "src/c_lib/lib/averaging_scheme.c",
#             "src/c_lib/sandbox/sandbox.pyx",
#         ],
#         include_dirs=include_dirs,
#         language="c",
#         libraries=libraries,
#         library_dirs=library_dirs,
#         extra_compile_args=extra_compile_args,
#         extra_link_args=extra_link_args,
#     )
# ]

setup(
    name="mrsimulator",
    version=version,
    description="A python toolbox for simulating NMR spectra",
    long_description=open(join(module_dir, "README.md")).read(),
    author="Deepansh J. Srivastava",
    author_email="deepansh2012@gmail.com",
    python_requires=">=3.6",
    url="https://github.com/DeepanshS/MRsimulator/",
    packages=find_packages("src"),
    package_dir={"": "src"},
    setup_requires=["numpy>=1.13.3", "setuptools>=27.3", "cython>=0.29.11"],
    install_requires=[
        "numpy>=1.13.3",
        "setuptools>=27.3",
        "cython>=0.29.11",
        "astropy>=3.0",
        "pydantic>=0.28",
        "requests>=2.21.0",
        "monty==2.0.4",
        "csdmpy>=0.1.2",
    ],
    ext_modules=cythonize(ext_modules, language_level=3),
    include_package_data=True,
    zip_safe=False,
    license="BSD-3-Clause",
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
    ],
)
