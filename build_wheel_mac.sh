#!/bin/bash

source /Users/deepansh/anaconda3/bin/activate

# python 3.6 - 3.9
for PYBIN in py36 py37 py38 py39; do
    conda activate "${PYBIN}"
    python --version
    cat requirements-dev.txt | sed -e '/^\s*#.*$/d' -e '/^\s*$/d' | xargs -n 1 pip install
    python setup.py develop bdist_wheel -d macwheels
    pytest
done


# Bundle external shared libraries into the wheels
for whl in macwheels/*.whl; do
    delocate-wheel "$whl" -w dist/
done

# generate source dist
python setup.py develop sdist

# clean up
rm -r build
rm macwheels/*.whl
rm src/mrsimulator/*-darwin.so
rm src/mrsimulator/tests/*-darwin.so
