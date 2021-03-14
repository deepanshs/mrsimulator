#!/bin/bash

source /Users/deepansh/opt/anaconda3/bin/activate
export DYLD_FALLBACK_LIBRARY_PATH="/usr/local/opt/fftw/lib:/usr/local/opt/openblas/lib"
export PATH="/Users/deepansh/opt/anaconda3/bin:$PATH"
export MACOSX_DEPLOYMENT_TARGET=10.9

# python 3.7 - 3.9
for PYBIN in py37 py38 py39; do
    conda activate "${PYBIN}"
    python --version
    < requirements-dev.txt sed -e '/^\s*#.*$/d' -e '/^\s*$/d' | xargs -n 1 python -m pip install
    python setup.py develop bdist_wheel -d macwheels
    pytest
done

# Bundle external shared libraries into the wheels
echo $DYLD_FALLBACK_LIBRARY_PATH
for whl in macwheels/*.whl; do
    delocate-wheel "$whl" -w dist/ -v
done

# generate source dist
python setup.py develop sdist

# clean up
rm -r build
rm macwheels/*.whl
rm src/mrsimulator/*-darwin.so
rm src/mrsimulator/tests/*-darwin.so
