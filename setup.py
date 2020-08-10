# -*- coding: utf-8 -*-
import platform
import sys
from os.path import abspath
from os.path import dirname
from os.path import exists
from os.path import join

import numpy as np
import numpy.distutils.system_info as sysinfo
from setuptools import Extension
from setuptools import find_packages
from setuptools import setup

from settings import use_accelerate
from settings import use_openblas

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
    print(f"Python>=3.6 is required for the setup. You are using version {py_version}")
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
extra_compile_args = [
    "-O3",
    "-ffast-math",
    # "-msse4.2",
    # "-ftree-vectorize",
    # "-fopt-info-vec-optimized",
    # "-mavx",
]
extra_link_args = []
data_files = []

numpy_include = np.get_include()

if sys.platform.startswith("win"):

    extra_link_args += ["-Wl,--allow-multiple-definition"]
    extra_compile_args = ["-DFFTW_DLL"]

    # FFTW3 info
    fftw3_info = sysinfo.get_info("fftw3")
    fftw_keys = fftw3_info.keys()
    if "include_dirs" in fftw_keys:
        include_dirs += fftw3_info["include_dirs"]
    if "library_dirs" in fftw_keys:
        library_dirs += fftw3_info["library_dirs"]
    if "libraries" in fftw_keys:
        libraries += fftw3_info["libraries"]

    # OpenBLAS info
    openblas_info = sysinfo.get_info("openblas")
    if openblas_info != {}:
        include_dirs += [join(fftw3_info["include_dirs"][0], "openblas")]
        library_dirs += openblas_info["library_dirs"]
        libraries += openblas_info["libraries"]


def message(lib):
    return f"Please install {lib} from homebrew with:\n\t$ brew install {lib}"


if platform.system() == "Darwin":  # OSX-specific tweaks:
    # BLAS framework

    # Apple's Accelerate framework for BLAS:
    if use_accelerate:
        acc_info = sysinfo.get_info("accelerate")
        if "extra_compile_args" in acc_info:
            extra_compile_args += acc_info["extra_compile_args"]
        if "extra_link_args" in acc_info:
            extra_link_args += acc_info["extra_link_args"]

    # OpenBLAS framework
    if use_openblas:
        BLAS_INCLUDE = "/usr/local/opt/openblas/include"
        BLAS_LIB = "/usr/local/opt/openblas/lib"
        libraries += ["openblas"]

        if not exists(BLAS_INCLUDE):
            print(message("openblas"))
            sys.exit(1)

        include_dirs += [BLAS_INCLUDE]
        library_dirs += [BLAS_LIB]

    # # MKL framework
    # if use_mkl:
    #     mkl_info = np.__config__.blas_mkl_info
    #     if mkl_info == {}:
    #         print("Please enable mkl for numpy before proceeding.")
    #         sys.exit(1)

    #     BLAS_INCLUDE = mkl_info["include_dirs"]
    #     BLAS_LIB = mkl_info["library_dirs"]
    #     libraries += mkl_info["libraries"]

    #     include_dirs += BLAS_INCLUDE
    #     library_dirs += BLAS_LIB

    # FFTW framework
    FFTW_INCLUDE = "/usr/local/opt/fftw/include"
    FFTW_LIB = "/usr/local/opt/fftw/lib"
    libraries += ["fftw3"]

    if not exists(FFTW_INCLUDE):
        print(message("fftw"))
        sys.exit(1)

    include_dirs += [FFTW_INCLUDE]
    library_dirs += [FFTW_LIB]

    # if USE_SSE_AVX:
    #     extra_compile_args += ["-Wa,-q"]

if platform.system() == "Linux":
    include_dirs += [
        "/usr/include/",
        "/usr/include/openblas",
        "/usr/include/x86_64-linux-gnu/",
    ]

    library_dirs += ["/usr/lib64/", "/usr/lib/", "/usr/lib/x86_64-linux-gnu/"]
    libraries += ["openblas", "fftw3"]
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

    extra_link_args += ["-lm"]
    extra_compile_args += ["-g"]

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
    "src/c_lib/lib/powder_setup.c",
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

extras = {"lmfit": ["lmfit>=1.0.0"], "matplotlib": ["matplotlib>=3.0"]}

description = "A python toolbox for simulating fast real-time solid-state NMR spectra."
setup(
    name="mrsimulator",
    version=version,
    description=description,
    long_description=open(join(module_dir, "README.md")).read(),
    long_description_content_type="text/markdown",
    author="Deepansh J. Srivastava",
    author_email="deepansh2012@gmail.com",
    python_requires=">=3.6",
    url="https://github.com/DeepanshS/MRsimulator/",
    packages=find_packages("src"),
    package_dir={"": "src"},
    setup_requires=["numpy>=1.17"],
    install_requires=[
        "numpy>=1.17",
        "csdmpy>=0.3.2",
        "pydantic>=1.0",
        "monty>=2.0.4",
        "typing-extensions>=3.7",
    ],
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
        "Development Status :: 4 - Beta",
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
