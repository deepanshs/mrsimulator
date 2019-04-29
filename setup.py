import numpy
import io
from setuptools import setup, Extension, find_packages
from MRsimulator.__version__ import __version__ 
try:
      from Cython.Distutils import build_ext
      from Cython.Build import cythonize

      # from Cython.Compiler import Options
      # Options.docstrings = True

except ImportError:
    use_cython = False
else:
    use_cython = True

from os import path, listdir
from os.path import isfile, join, abspath, dirname


# Package meta-data.
NAME = 'MRsimulator'
DESCRIPTION = 'A python toolbox for simulating NMR spectra'
URL = 'https://github.com/DeepanshS/MRsimulator/'
EMAIL = 'srivastava.89@osu.edu'
AUTHOR = 'Deepansh J. Srivastava'
REQUIRES_PYTHON = '>=3.5'
VERSION = __version__


# What packages are required for this module to be executed?
REQUIRED = [
    'scipy>=0.16.0',
    'numpy>=1.10.1'
]

# What packages are optional?
EXTRAS = {
     'fancy feature': ['matplotlib>=3.0.2'],
}


here = abspath(dirname(__file__))
# Import the README and use it as the long-description.
# Note: this will only work if 'README.md' is present in your MANIFEST.in file!
try:
    with io.open(join(here, 'README.md'), encoding='utf-8') as f:
        long_description = '\n' + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

# Load the package's __version__.py module as a dictionary.
about = {}
if not VERSION:
    project_slug = NAME.lower().replace("-", "_").replace(" ", "_")
    with open(join(here, project_slug, '__version__.py')) as f:
        exec(f.read(), about)
else:
    about['__version__'] = VERSION





## Cython =========================================================== ##
cmdclass = {}
ext_modules = [ ]


if use_cython:
      cmdclass.update({ 'build_ext': build_ext })


# scr_c folder
nmr_lib_source_dir = 'MRsimulator/scr/mrlib/scr_c/'
_source_files = []
for _file in listdir(nmr_lib_source_dir):
      if _file.endswith(".c"):
            _source_files.append(nmr_lib_source_dir+_file)


# powder_averaging_scheme
nmr_lib_source_dir = 'MRsimulator/scr/mrlib/scr_c/powder_averaging_scheme/'
for _file in listdir(nmr_lib_source_dir):
      if _file.endswith(".c"):
            _source_files.append(nmr_lib_source_dir+_file)


include_nmr_lib_directories = ["MRsimulator/scr/mrlib/include/"]
include_nmr_lib_directories.append(numpy.get_include())

nmr_lib_source_file = _source_files[:]
if use_cython:
      # pyx
      nmr_lib_source_file.append("MRsimulator/scr/mrlib/mrlib.pyx")
else:
      nmr_lib_source_file.append("MRsimulator/scr/mrlib/mrlib.c")


print("NMR lib Source files----------------------------------")
for item in nmr_lib_source_file:
      print(item)


ext_modules += [
      Extension(
            name = NAME+'.lib',
            sources = nmr_lib_source_file,
            include_dirs = include_nmr_lib_directories,
            language="c",
            extra_compile_args = "-flax-vector-conversions -g -Ofast".split(),
            extra_link_args = "-g".split()
      )
]


nmr_function_source_file = _source_files[:]
nmr_function_source_dir = 'MRsimulator/scr/mrmethods/'

for _file in listdir(nmr_function_source_dir):
      if _file.endswith(".c") and _file != 'nmr_methods.c':
            nmr_function_source_file.append(nmr_function_source_dir+_file)

if use_cython:
      nmr_function_source_file.append("MRsimulator/scr/mrmethods/nmr_methods.pyx")
      nmr_function_source_file.append("MRsimulator/scr/mrlib/mrlib.pxd")
else:
      nmr_function_source_file.append("MRsimulator/scr/mrmethods/nmr_methods.c")
      nmr_function_source_file.append("MRsimulator/scr/mrlib/mrlib.c")

include_nmr_lib_directories.append("MRsimulator/scr/mrmethods/include")

print (include_nmr_lib_directories)
print("NMR method Source files----------------------------------")
for item in nmr_function_source_file:
      print(item)


ext_modules += [
      Extension(
            name=NAME+'.methods',
            sources=nmr_function_source_file,
            # include_dirs=[numpy.get_include()],
            # extra_objects= ['./mrlib'], # ["fc.o"],  # if you compile fc.cpp separately
            include_dirs = include_nmr_lib_directories,  # .../site-packages/numpy/core/include
            language="c",
            # libraries=["./mrlib"],
            # -ffast-math -flax-vector-conversions -g -Ofast
            extra_compile_args = "-O1".split(), # 
            extra_link_args = "-g -lfftw3 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
                              -ldl -liomp5 -lm -Wl".split() #  
      )
]


if use_cython:
      ext = cythonize(ext_modules, annotate=True, language_level=3, gdb_debug=True)
else:
      ext = ext_modules

setup(
      name = NAME,
      version=about['__version__'],
      description=DESCRIPTION,
      long_description=long_description,

      author=AUTHOR,
      author_email=EMAIL,
      python_requires=REQUIRES_PYTHON,
      url=URL,
      packages=find_packages(),

      install_requires=REQUIRED,
      extras_require=EXTRAS,

      cmdclass = cmdclass,
      ext_modules = ext, # cythonize(ext_modules, annotate=True, language_level=3, gdb_debug=True),
      # ext_modules = cythonize(ext_modules)  ? not in 0.14.1
      classifiers=[
        # Trove classifiers
        # Full list: https://pypi.python.org/pypi?%3Aaction=list_classifiers
        'License :: OSI Approved :: BSD License',
        'Programming Language :: Python :: 3',
    ],
)