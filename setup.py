import platform
import sys
from os.path import abspath
from os.path import dirname
from os.path import exists
from os.path import join
from os import environ
import warnings
from setuptools import Extension
from setuptools import find_packages
from setuptools import setup

import numpy as np

from settings import use_accelerate
from settings import use_openblas

try:
    from Cython.Build import cythonize

    USE_CYTHON = True
except ImportError:
    USE_CYTHON = False

PY_MINOR_VERSION = 8


def message(lib, env, command, key):
    arg = f"{key} {lib}" if key != "" else f"{lib}"
    warning = (
        f"\nLibraries not found - {lib}.\n",
        "Use environ variable to add the path to the include and lib folders. ",
        "For example,\n",
        '\texport LDFLAGS="-L/usr/local/opt/openblas/lib"\n',
        '\texport CPPFLAGS="-I/usr/local/opt/openblas/include"\n',
        f"\nYou can also try installing '{lib}' from {env} with: ",
        f"\n\t{command} install {arg}\n",
    )
    warnings.warn("".join(warning))


class Setup:
    __slots__ = [
        "include_dirs",
        "library_dirs",
        "libraries",
        "extra_compile_args",
        "extra_link_args",
    ]

    def __init__(self):
        self.libraries = []
        self.include_dirs = []
        self.library_dirs = []
        self.extra_compile_args = []
        self.extra_link_args = []

    def check_valid_path(self, pathlist):
        return [pth for pth in pathlist if exists(pth)]

    def check_if_lib_exists(self, lib):
        return np.any([exists(join(pth, lib)) for pth in self.library_dirs])

    def check_if_header_exists(self, header):
        return np.any([exists(join(pth, header)) for pth in self.include_dirs])

    def conda_setup_for_windows(self):
        self.libraries += ["fftw3", "openblas"]
        self.extra_compile_args = ["/DUSE_OPENBLAS"]

        print(sys.version)
        loc = dirname(sys.executable)
        print("executable location", loc)
        if "conda" not in loc:
            return

        print("Found Python installation:", loc)
        self.include_dirs += self.check_valid_path(
            [
                join(loc, "Library", "include", "fftw"),
                join(loc, "Library", "include", "openblas"),
                join(loc, "Library", "include"),
                join(loc, "include"),
            ]
        )
        self.library_dirs += self.check_valid_path([join(loc, "Library", "lib")])
        environ["MRSIM_LIB"] = str(join(loc, "Library", "lib"))
        self.on_exit_message("openblas.lib", "fftw3.lib")

    def conda_setup_for_unix(self):
        loc = dirname(sys.executable)
        print("Found Python installation:", loc)

        self.include_dirs += self.check_valid_path([join(loc, "include")])
        self.library_dirs += self.check_valid_path([join(loc, "lib")])
        self.extra_compile_args = ["-O3", "-ffast-math", "-DUSE_OPENBLAS"]
        self.libraries += ["fftw3", "openblas"]

    def on_exit_message(self, blas_lib, fftw_lib):
        found_blas = self.check_if_lib_exists(blas_lib)
        found_fftw = self.check_if_lib_exists(fftw_lib)

        cmd_list = ["conda", "conda", "-c conda-forge"]
        if not found_blas and not found_fftw:
            message("openblas fftw", *cmd_list)

        if not found_blas:
            message("openblas", *cmd_list)

        if not found_fftw:
            message("fftw", *cmd_list)

    def mkl_blas_info(self):
        mkl_info = np.__config__.blas_mkl_info
        if mkl_info == {}:
            print("Please enable mkl for numpy before proceeding.")
            message("mkl mkl-include", "pip", "pip", "")

        self.include_dirs += mkl_info["include_dirs"]
        self.library_dirs += mkl_info["library_dirs"]
        self.libraries += mkl_info["libraries"]

        if not self.check_if_header_exists("mkl.h"):
            print("mkl header file not found.")
            message("mkl-include", "pip", "pip", "")

        print("Attempting to link mrsimulator with the mkl blas.")
        self.extra_compile_args += ["-DUSE_MKL", "/DUSE_MKL"]


class WindowsSetup(Setup):
    def __init__(self):
        super().__init__()

        self.extra_link_args += ["-Wl"]
        self.extra_compile_args = ["-DFFTW_DLL"]

        # if use_mkl:
        #     self.mkl_blas_info()

        self.conda_setup_for_windows()


class LinuxSetup(Setup):
    def __init__(self):
        super().__init__()
        self.extra_compile_args = [
            "-O3",
            "-ffast-math",
            "-fcommon",
            # "-msse4.2",
            # "-ftree-vectorize",
            # "-fopt-info-vec-all",
            # "-fopt-info-vec-optimized",
            # "-mavx",
            "-g",
            "-DUSE_OPENBLAS",
        ]
        self.extra_link_args += ["-lm"]
        self.include_dirs += [
            "/usr/include/",
            "/usr/include/openblas",
            "/usr/include/x86_64-linux-gnu/",
        ]

        self.library_dirs += ["/usr/lib64/", "/usr/lib/", "/usr/lib/x86_64-linux-gnu/"]
        self.libraries += ["openblas", "fftw3"]

    def message(self, lib_centos, lib_ubuntu):
        print(f"Warning: {lib_ubuntu} might not be installed. See below.\n")
        stat = f"yum install {lib_centos}"
        print(f"For CentOS users\n\t{stat}")
        stat = f"sudo apt-get update\n\tsudo apt-get install {lib_ubuntu}"
        print(f"For Ubuntu users\n\t{stat}")


class MacOSSetup(Setup):
    def __init__(self):
        super().__init__()
        self.extra_compile_args = [
            "-O3",
            "-ffast-math",
            # "-Rpass=loop-vectorize",
            # "-Rpass-missed=loop-vectorize",
            # "-Rpass-analysis=loop-vectorize",
            "-fvectorize",
            "-fcommon",
        ]
        self.extra_link_args += ["-lm"]

        # Blas
        if use_accelerate:
            self.accelerate_info()

        if use_openblas:
            self.openblas_info()

        # FFTW
        self.fftw_info()

    def accelerate_info(self):
        """Apple's Accelerate framework for BLAS"""
        self.extra_compile_args += [
            "-I/System/Library/Frameworks/vecLib.framework/Headers",
            "-Wl,-framework",
            "-Wl,Accelerate",
        ]
        print("Attempting to link mrsimulator with the Apple accelerate library.")
        self.extra_compile_args += ["-DUSE_ACCELERATE"]

    def openblas_info(self):
        """openblas includes and lib are for brew installation"""
        blas_include_dir = [
            "/usr/local/opt/openblas/include",
            "/opt/homebrew/opt/openblas/include",
        ]
        blas_library_dir = [
            "/usr/local/opt/openblas/lib",
            "/opt/homebrew/opt/openblas/lib",
        ]
        blas_library = "openblas"

        exists_all = [exists(item) for item in blas_include_dir]
        if not any(exists_all):
            message("openblas", "homebrew", "brew", "")

        print("Attempting to link mrsimulator with the openblas library.")
        self.include_dirs += blas_include_dir
        self.library_dirs += blas_library_dir
        self.libraries += [blas_library]
        self.extra_compile_args += ["-DUSE_OPENBLAS"]

    def fftw_info(self):
        """fftw includes and lib are for brew installation"""
        fftw_include_dir = [
            "/usr/local/opt/fftw/include",
            "/opt/homebrew/opt/fftw/include",
        ]
        fftw_library_dir = ["/usr/local/opt/fftw/lib", "/opt/homebrew/opt/fftw/lib"]
        fftw_library = "fftw3"

        exists_all = [exists(item) for item in fftw_include_dir]
        if not any(exists_all):
            message("fftw", "homebrew", "brew", "")

        print("Attempting to link mrsimulator with the fftw library.")
        self.include_dirs += fftw_include_dir
        self.library_dirs += fftw_library_dir
        self.libraries += [fftw_library]


# get the version from file
python_version = sys.version_info
py_version = ".".join([str(i) for i in python_version[:3]])
print("Using python version", py_version)
if python_version.major != 3 and python_version.minor < PY_MINOR_VERSION:
    print(f"Python>=3.{PY_MINOR_VERSION} is required, found {py_version}")
    sys.exit(1)

with open("src/mrsimulator/__init__.py") as f:
    for line in f.readlines():
        if "__version__" in line:
            before_keyword, keyword, after_keyword = line.partition("=")
            version = after_keyword.strip()[1:-1]
            print("mrsimulator version ", version)
            break

module_dir = dirname(abspath(__file__))
data_files = []
numpy_include = np.get_include()

if sys.platform.startswith("win"):
    win = WindowsSetup()

if platform.system() == "Darwin":
    win = MacOSSetup()

if platform.system() == "Linux":
    win = LinuxSetup()

extra_link_args = list(set(win.extra_link_args))
extra_compile_args = list(set(win.extra_compile_args))
library_dirs = list(set(win.library_dirs))
include_dirs = list(set(win.include_dirs))
libraries = list(set(win.libraries))

# other include paths
include_dirs += ["src/c_lib/include/", numpy_include]

# print info
print("include_dirs", include_dirs)
print("library_dirs", library_dirs)
print("libraries", libraries)
print("extra_compile_args", extra_compile_args)
print("extra_link_args", extra_link_args)

source = [
    "src/c_lib/lib/angular_momentum/wigner_element.c",
    "src/c_lib/lib/angular_momentum/wigner_matrix.c",
    "src/c_lib/lib/interpolation.c",
    "src/c_lib/lib/method.c",
    "src/c_lib/lib/mrsimulator.c",
    "src/c_lib/lib/octahedron.c",
    "src/c_lib/lib/frequency_averaging.c",
    "src/c_lib/lib/schemes.c",
    "src/c_lib/lib/simulation.c",
    "src/c_lib/lib/vm_linalg.c",
]

print("USE_CYTHON", USE_CYTHON)
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

# c-lib
ext_modules += [
    Extension(
        name="mrsimulator.clib",
        sources=["src/c_lib/lib/histogram.c", "src/c_lib/clib/clib" + ext],
        include_dirs=include_dirs,
        language="c",
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
    )
]

if USE_CYTHON:
    ext_modules = cythonize(ext_modules, language_level=3, gdb_debug=False)

extras = {}  # {"all": ["matplotlib>=3.3.4"]}

description = "A python toolbox for simulating fast real-time solid-state NMR spectra."
setup(
    name="mrsimulator",
    version=version,
    description=description,
    long_description=open(join(module_dir, "README.md")).read(),
    long_description_content_type="text/markdown",
    author="Deepansh J. Srivastava",
    author_email="srivastava.89@osu.edu",
    python_requires=">=3.10",
    url="https://github.com/deepanshs/mrsimulator/",
    packages=find_packages("src"),
    package_dir={"": "src"},
    setup_requires=["numpy>=1.20"],
    install_requires=[
        "numpy>=1.20",
        "csdmpy>=0.7",
        "pydantic>=2.8",
        "typing-extensions>=3.7",
        "psutil>=5.4.8",
        "joblib>=1.0.0",
        "pandas>=1.1.3",
        "lmfit>=1.0.2",
        "matplotlib>=3.3.4",
        "astropy>=6",
    ],
    entry_points={"console_scripts": ["mrsimulator=mrsimulator.__main__:run"]},
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
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "Programming Language :: Python :: 3.13",
        "Topic :: Education",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Topic :: Scientific/Engineering :: Physics",
    ],
)
