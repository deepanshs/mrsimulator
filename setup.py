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

import json

# get the version from file

with open("src/mrsimulator/__init__.py", "r") as f:
    for line in f.readlines():
        if "__version__" in line:
            before_keyword, keyword, after_keyword = line.partition("=")
            version = after_keyword.strip()[1:-1]

module_dir = dirname(abspath(__file__))

include_dirs = []
library_dirs = []
libraries = []
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
    # data_files += [join(conda_location, "Library", "bin", "mkl_rt.dll")]
    libraries += ["fftw3", "openblas"]
    name = "openblas"

    extra_link_args = ["-lm"]
    extra_compile_args = ["-DFFTW_DLL"]

else:
    openblas_info = sysinfo.get_info("openblas")
    fftw3_info = sysinfo.get_info("fftw3")
    # mkl_info = sysinfo.get_info("mkl")

    # if mkl_info != {}:
    #     name = "mkl"
    #     include_dirs += mkl_info["include_dirs"]
    #     library_dirs += mkl_info["library_dirs"]
    #     libraries += mkl_info["libraries"]
    if openblas_info != {}:
        name = "openblas"
        library_dirs += openblas_info["library_dirs"]
        libraries += openblas_info["libraries"]
    # else:
    #     raise Exception("mkl blas or openblas library not found.")

    include_dirs += fftw3_info["include_dirs"]
    library_dirs += fftw3_info["library_dirs"]
    libraries += fftw3_info["libraries"]

    extra_link_args = ["-lm"]
    extra_compile_args = ["-g", "-O3"]

    libraries += ["fftw3_threads", "pthread"]

include_dirs = list(set(include_dirs))
library_dirs = list(set(library_dirs))
libraries = list(set(libraries))


blas_info = {
    "name": name,
    "library_dirs": library_dirs,
    "include_dirs": include_dirs,
    "libraries": libraries,
}

print(blas_info)

with open("src/mrsimulator/__config__.json", "w", encoding="utf8") as outfile:
    json.dump(blas_info, outfile, ensure_ascii=True, indent=2)

# other include paths
include_dirs += ["src/c_lib/include", numpy_include]

# system = platform.system()
# arch = platform.architecture()[0]
# compiler = platform.python_compiler()
# if system == 'Linux':
#     extra_link_args += ["-lm", "-ldl"]
#     if arch == '64bit':
#         extra_compile_args += ["-m64", "-DMKL_ILP64"]
#     if arch == '32bit':
#         extra_compile_args += ["-m32"]
# if system == 'Darwin':
#     extra_link_args += ["-Wl", "-lm", "-ldl"]
#     extra_compile_args += ["-m64"]

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
        "matplotlib>=3.0.2",
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
