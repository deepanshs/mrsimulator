#!/bin/bash

source /Users/deepansh/opt/anaconda3/bin/activate

# python 3.6, 3.7
for PYBIN in py36 py37; do
    conda activate "${PYBIN}"
    python --version
    pip install -r requirements.txt
    pip install -r requirements-dev.txt
    python setup.py develop
    pytest
    python setup.py bdist_wheel -d macwheels
    delocate-wheel -w dist/ macwheels/*.whl
    rm macwheels/*.whl
done

# clean up
rm -r build
rm src/mrsimulator/*-darwin.so
rm src/mrsimulator/tests/*-darwin.so
