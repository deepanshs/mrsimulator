# -*- coding: utf-8 -*-
import io

from setuptools import Extension
from setuptools import find_packages
from setuptools import setup

try:
    from Cython.Distutils import build_ext
    from Cython.Build import cythonize

except ImportError:
    use_cython = False
else:
    use_cython = True

import platform
from os import path, listdir
from os.path import join, abspath, dirname


# Package meta-data.
NAME = "mrsimulator"
DESCRIPTION = "A python toolbox for simulating NMR spectra"
URL = "https://github.com/DeepanshS/mrsimulator/"
EMAIL = "srivastava.89@osu.edu"
AUTHOR = "Deepansh J. Srivastava"
REQUIRES_PYTHON = ">=3.0"
VERSION = "0.1.1"


# What packages are required for this module to be executed?
REQUIRED = ["numpy>=1.13.3", "astropy>=3.0", "mkl", "mkl-include"]

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


if use_cython:
    cmdclass.update({"build_ext": build_ext})


# scr_c folder
_list = ["mrsimulator", "scr", "lib"]
nmr_lib_source_dir = path.join(*_list)
_source_files = []
for _file in listdir(nmr_lib_source_dir):
    if _file.endswith(".c"):
        _source_files.append(path.join(nmr_lib_source_dir, _file))

_list = ["mrsimulator", "scr", "include"]
include_nmr_lib_directories = [path.join(*_list)]


def get_numpy_includes():
    import numpy

    path_ = numpy.get_include()
    if platform.system() == "Windows":
        for i in range(5):
            path_ = path.split(path_)[0]
        path_ = path.join(path_, "Library", "include")
        include_nmr_lib_directories.append(path_)
    return path_


include_nmr_lib_directories.append(get_numpy_includes())


nmr_function_source_file = _source_files[:]

_list = ["mrsimulator", "scr", "mrmethods"]
nmr_function_source_dir = path.join(*_list)

for _file in listdir(nmr_function_source_dir):
    if _file.endswith(".c") and _file != "nmr_methods.c":
        nmr_function_source_file.append(
            path.join(nmr_function_source_dir, _file)
        )


source2 = nmr_function_source_file[:]

if use_cython:
    _list = ["mrsimulator", "scr", "mrmethods", "nmr_methods.pyx"]
    nmr_function_source_file.append(path.join(*_list))

else:
    _list = ["mrsimulator", "scr", "mrmethods", "nmr_methods.c"]
    nmr_function_source_file.append(path.join(*_list))


print(include_nmr_lib_directories)
print("NMR method Source files----------------------------------")
for item in nmr_function_source_file:
    print(item)


ext_modules = [
    Extension(
        name=NAME + ".methods",
        sources=nmr_function_source_file,
        # include_dirs=[numpy.get_include()],
        # # if you compile fc.cpp separately
        # extra_objects= ['./mrlib'], # ["fc.o"],
        #  # .../site-packages/numpy/core/include
        include_dirs=include_nmr_lib_directories,
        language="c",
        # libraries=["./mrlib"],
        # -ffast-math -flax-vector-conversions -g -Ofast
        extra_compile_args="-O1".split(),  #
        extra_link_args="-g -lfftw3 -lmkl_intel_lp64 -lmkl_intel_thread \
                        -lmkl_core -ldl -liomp5 -lm -Wl".split(),  #
    )
]


# Sandbox

# sandbox_files = source2[:]

# _list = [
#     'mrsimulator', 'scr', 'sandbox', 'interpolation', 'interpolation.pyx'
# ]
# sandbox_files.append(path.join(*_list))

# ext_modules += [
#     Extension(
#         name=NAME+'.sandbox.interpolation',
#         sources=sandbox_files,
#         # # include_dirs=[numpy.get_include()],
#         # if you compile fc.cpp separately
#         # extra_objects= ['./mrlib'], # ["fc.o"],
#         # # .../site-packages/numpy/core/include
#         include_dirs=include_nmr_lib_directories,
#         language="c",
#         # libraries=["./mrlib"],
#         # -ffast-math -flax-vector-conversions -g -Ofast
#         extra_compile_args="-O1".split(),
#         extra_link_args="-g -lfftw3 -lmkl_intel_lp64 -lmkl_intel_thread \
#                         -lmkl_core -ldl -liomp5 -lm -Wl".split()
#     )
# ]

if use_cython:
    ext = cythonize(
        ext_modules, annotate=True, language_level=3, gdb_debug=True
    )
else:
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
    install_requires=REQUIRED,
    extras_require=EXTRAS,
    # data_files = ['mrsimulator/test/isotopomers.json'],
    cmdclass=cmdclass,
    # cythonize(ext_modules, annotate=True, language_level=3, gdb_debug=True),
    ext_modules=ext,
    # ext_modules = cythonize(ext_modules)  ? not in 0.14.1
    classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3",
    ],
)
