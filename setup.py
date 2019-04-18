import numpy
from setuptools import setup, Extension
# from distutils.core import setup
# from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

from os import path
import os
from Cython.Compiler import Options

Options.docstrings = True


from os import listdir
from os.path import isfile, join


library_name = "MRsimulator"
# os.environ['CC'] = 'icc'

# scr_c folder
nmr_lib_source_dir = 'MRsimulator/mrlib/scr_c/'
_source_files = []
for _file in listdir(nmr_lib_source_dir):
      if _file.endswith(".c"):
            _source_files.append(nmr_lib_source_dir+_file)

# powder_averaging_scheme
nmr_lib_source_dir = 'MRsimulator/mrlib/scr_c/powder_averaging_scheme/'
for _file in listdir(nmr_lib_source_dir):
      if _file.endswith(".c"):
            _source_files.append(nmr_lib_source_dir+_file)

# pyx
# nmr_lib_source_dir = 'scr/mrlib/'
# for _file in listdir(nmr_lib_source_dir):
#       if _file.endswith(".pyx"):
#             _source_files.append(nmr_lib_source_dir+_file)
include_nmr_lib_directories = ["MRsimulator/mrlib/include/"]
include_nmr_lib_directories.append(numpy.get_include())

nmr_lib_source_file = _source_files[:]
nmr_lib_source_file.append("MRsimulator/mrlib/mrlib.pyx")
print("NMR lib Source files----------------------------------")
for item in nmr_lib_source_file:
      print(item)


ext_modules = [
      Extension(
            name = 'nmr.lib',
            sources = nmr_lib_source_file,
            include_dirs = include_nmr_lib_directories,
            language="c",
            extra_compile_args = "-flax-vector-conversions -g -Ofast".split(),
            extra_link_args = "-g".split()
      )
]

# setup(
#       # name = 'mrlib',
#       cmdclass = {'build_ext': build_ext},
#       ext_modules = cythonize(ext_modules, annotate=True, language_level="3"),
#       # ext_modules = cythonize(ext_modules)  ? not in 0.14.1
#       version = '0.0.1dev1',
#       description = 'NMR library',
#       author = 'Deepansh J. Srivastava',
#       author_email= 'srivastava.89@osu.edu'
#       )


nmr_function_source_file = _source_files[:]
nmr_function_source_dir = 'MRsimulator/mrmethods/'
for _file in listdir(nmr_function_source_dir):
      if _file.endswith(".c") and _file != 'nmr_methods.c':
            nmr_function_source_file.append(nmr_function_source_dir+_file)
nmr_function_source_file.append("MRsimulator/mrmethods/nmr_methods.pyx")
nmr_function_source_file.append("MRsimulator/mrlib/mrlib.pxd")
include_nmr_lib_directories.append("MRsimulator/mrmethods/include")

print (include_nmr_lib_directories)

print("NMR method Source files----------------------------------")
for item in nmr_function_source_file:
      print(item)



ext_modules.append(
      Extension(
            name='nmr.methods',
            sources=nmr_function_source_file,
            # include_dirs=[numpy.get_include()],
            # extra_objects= ['./mrlib'], # ["fc.o"],  # if you compile fc.cpp separately
            include_dirs = include_nmr_lib_directories,  # .../site-packages/numpy/core/include
            language="c",
            # libraries=["./mrlib"],
            extra_compile_args = "-ffast-math -flax-vector-conversions -g -Ofast".split(), # 
            extra_link_args = "-g -lfftw3 -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core \
                              -ldl -liomp5 -lm -Wl".split() #  
      )
)
 
setup(
      # name = 'MRmethods',
      cmdclass = {'build_ext': build_ext},
      ext_modules = cythonize(ext_modules, annotate=True, language_level=3, gdb_debug=True),
      # ext_modules = cythonize(ext_modules)  ? not in 0.14.1
      version = '0.0.1dev1',
      description = 'A toolbox for simulating NMR spectra',
      author = 'Deepansh J. Srivastava',
      author_email= 'srivastava.89@osu.edu'
      )