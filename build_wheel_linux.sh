#!/bin/bash
set -e -x

# used with docker image - quay.io/pypa/manylinux2010_x86_64
# Install a system package required by our library
yum install -y openblas-devel git fftw-devel
mkdir mrsim
git clone --single-branch --branch master https://github.com/DeepanshS/mrsimulator.git mrsim
cd mrsim

# Compile wheels
for PYBIN in /opt/python/cp3[6-7]*/bin; do
    "${PYBIN}/pip" install -r requirements.txt
    "${PYBIN}/pip" install -r requirements-dev.txt
    "${PYBIN}/python" setup.py develop
    "${PYBIN}/pytest"
    "${PYBIN}/python" setup.py bdist_wheel
done

# Bundle external shared libraries into the wheels
for whl in dist/*.whl; do
    auditwheel repair "$whl" -w /io/dist/
done
