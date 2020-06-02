#!/bin/bash

source /Users/deepansh/opt/anaconda3/bin/activate

# python 3.6, 3.7
for PYBIN in py36 py37 py38; do
    conda activate "${PYBIN}"
    python --version
    pip install -r requirements.txt
    python setup.py develop bdist_wheel -d macwheels
    pip install -r requirements-dev.txt
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
