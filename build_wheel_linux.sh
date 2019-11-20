#!/bin/bash
set -e -x
cd io

# used with docker image - quay.io/pypa/manylinux2010_x86_64
# Install a system package required by our library
yum install -y openblas-devel git fftw-devel

# Compile wheels
for PYBIN in /opt/python/cp3[6-7]*/bin; do
    "${PYBIN}/pip" install --upgrade pip
    "${PYBIN}/pip" install -r requirements.txt
    "${PYBIN}/pip" install -r requirements-dev.txt
    "${PYBIN}/python" setup.py develop
    "${PYBIN}/pytest"
    "${PYBIN}/python" setup.py bdist_wheel -d linuxwheels
    rm src/c_lib/mrmethods/nmr_methods.c
    rm src/c_lib/sandbox/sandbox.c
    rm src/c_lib/test/test.c
done

# Bundle external shared libraries into the wheels
for whl in linuxwheels/*.whl; do
    auditwheel repair "$whl" -w dist/
done

# clean up
rm src/mrsimulator/*-linux-gnu.so
rm src/mrsimulator/tests/*-linux-gnu.so
rm linuxwheels/*.whl
rm -r build
