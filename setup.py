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

from os import listdir
from os.path import join
from os.path import abspath
from os.path import dirname
from os.path import split
from os.path import exists
import platform
import numpy as np

accelerate = False
mkl = True

if mkl:
    include_path = np.__config__.blas_mkl_info["include_dirs"]
    include_path = [abspath(_) for _ in include_path]
if accelerate:
    include_path = ["/System/Library/Frameworks/Accelerate.framework"]

# Package meta-data.
NAME = "mrsimulator"
DESCRIPTION = "A python toolbox for simulating NMR spectra"
URL = "https://github.com/DeepanshS/mrsimulator/"
EMAIL = "srivastava.89@osu.edu"
AUTHOR = "Deepansh J. Srivastava"
REQUIRES_PYTHON = ">=3.0"
VERSION = "0.1.1"


# What packages are required for this module to be executed?
REQUIRED = [
    "numpy>=1.13.3",
    "astropy>=3.0",
    "mkl>=2019",
    "mkl-include>=2019",
    "requests>=2.21.0",
    "matplotlib>=3.0.2",
    "numba",
]

# What packages are optional?
EXTRAS = {
    "fancy feature": ["matplotlib>=3.0.2", "plotly>=3.6", "dash>=0.40", "dash_daq>=0.1"]
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
nmr_lib_source_dir = join(*_list)
_source_files = []
for _file in listdir(nmr_lib_source_dir):
    if _file.endswith(".c"):
        _source_files.append(join(nmr_lib_source_dir, _file))

_list = ["mrsimulator", "scr", "include"]
include_path.append(join(*_list))
include_path.append(np.get_include())

# def get_include_and_lib_paths(include_path, lib_path):
#     """Get paths to include and lib folders on windows, linux and mac os."""
#     import numpy

#     # numpy include
#     path_ = numpy.get_include()
#     print("numpy include:", path_)
#     include_path.append(path_)

#     # conda lib
#     for _ in range(5):
#         path_ = split(path_)[0]
#     path_lib = path_
#     print("conda lib:")
#     print("exist:", exists(path_lib))
#     print("path:", path_lib)
#     if exists(path_lib):
#         lib_path.append(path_lib)

#     # conda include
#     path_include_conda = join(split(path_)[0], "include")
#     print("conda include:")
#     print("exist:", exists(path_include_conda))
#     print("path:", path_include_conda)
#     if exists(path_include_conda):
#         include_path.append(path_include_conda)

#     # other include
#     path_include_other = join(path_lib, "Library", "include")
#     print("other include:")
#     print("exist:", exists(path_include_other))
#     print("path:", path_include_other)
#     if exists(path_include_other):
#         include_path.append(path_include_other)

#     path_lib = join(path_lib, "Library", "bin")
#     print("other bin:")
#     print("exist:", exists(path_lib))
#     print("path:", path_lib)
#     if exists(path_lib):
#         lib_path.append(path_lib)

lib_path = []
# get_include_and_lib_paths(include_path, lib_path)

nmr_function_source_file = _source_files[:]

_list = ["mrsimulator", "scr", "mrmethods"]
nmr_function_source_dir = join(*_list)

for _file in listdir(nmr_function_source_dir):
    if _file.endswith(".c") and _file != "nmr_methods.c":
        nmr_function_source_file.append(join(nmr_function_source_dir, _file))


source2 = nmr_function_source_file[:]

if use_cython:
    _list = ["mrsimulator", "scr", "mrmethods", "nmr_methods.pyx"]
    nmr_function_source_file.append(join(*_list))

else:
    _list = ["mrsimulator", "scr", "mrmethods", "nmr_methods.c"]
    nmr_function_source_file.append(join(*_list))


print(include_path)
print("NMR method Source files----------------------------------")
for item in nmr_function_source_file:
    print(item)


def get_config(platform_):
    if platform_ == "Windows":
        return Windows()
    if platform_ == "Darwin":
        return Darwin()


def Darwin():
    config = {}
    if accelerate:
        # include_path.append('/System/Library/Frameworks/Accelerate.framework')
        config["libraries"] = []
        config["library_dirs"] = []
        config["extra_link_args"] = ["-Wl,-framework", "-Wl,Accelerate"]
        config["extra_compile_args"] = ["-msse3", "-O1"]
    if mkl:
        config = {}
        # include_path.append(*np.__config__.blas_mkl_info['include_dirs'])
        config["libraries"] = ["mkl_intel_lp64", "mkl_core", "pthread"]
        config["library_dirs"] = []
        config["extra_link_args"] = ["-lm", "-ldl", "-W"]
        config["extra_compile_args"] = ["-m64", "-g", "-msse3", "-O1"]
    return config


def Windows():
    lib = ["mkl_intel_lp64_dll", "mkl_core_dll", "mkl_intel_thread_dll"]
    lib_path = [abspath(_) for _ in np.__config__.blas_mkl_info["library_dirs"]]
    if mkl:
        # inc_path = np.__config__.blas_mkl_info['include_dirs']
        # [include_path.append(_) for _ in inc_path]
        config = {}
        config["libraries"] = lib
        config["library_dirs"] = lib_path
        config["extra_link_args"] = []
        config["extra_compile_args"] = ["-O1"]
    return config


# print('system', platform.system())
# if platform.system() == 'Darwin':
#     framework = '/System/Library/Frameworks/Accelerate.framework'
#     include_path.append(framework)
#     extra_link_args = ['-Wl,-framework', '-Wl,Accelerate']
#     extra_compile_args = ['-msse3', '-O3']

# if platform.system() == 'Windows':
#     extra_link_args = ['-lmkl_intel_lp64', '-lmkl_sequential', '-lmkl_core',
#                        '-lpthread', '-lm', '-ldl', '-W']
#     extra_compile_args = ['-m64', '-g', '-msse3', 'O3']

config = get_config(platform.system())
print(config)

ext_modules = [
    Extension(
        name=NAME + ".methods",
        sources=nmr_function_source_file,
        include_dirs=include_path,
        language="c",
        libraries=config["libraries"],
        library_dirs=config["library_dirs"],
        extra_compile_args=config["extra_compile_args"],
        extra_link_args=config["extra_link_args"]
        # extra_link_args="-g -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
        #                  -liomp5 -lpthread -lm -ldl -W".split(),
    )
]


# Sandbox

sandbox_files = source2[:]

_list = ["mrsimulator", "scr", "sandbox", "sandbox.pyx"]
sandbox_files.append(join(*_list))

ext_modules += [
    Extension(
        name=NAME + ".sandbox",
        sources=sandbox_files,
        include_dirs=include_path,
        language="c",
        extra_compile_args=config["extra_compile_args"],
        extra_link_args=config["extra_link_args"]
        # extra_link_args="-g -lmkl_intel_lp64 -lmkl_intel_thread \
        #                 -lmkl_core -ldl -liomp5 -lpthread -lm -W".split()
        # extra_link_args="-g -lmkl_intel_lp64 -lmkl_intel_thread \
        #                 -lmkl_core -ldl -liomp5 -lpthread -lm -Wl".split(),
    )
]

if use_cython:
    ext = cythonize(ext_modules, annotate=True, language_level=3, gdb_debug=True)
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
    license="BSD-3-Clause",
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
