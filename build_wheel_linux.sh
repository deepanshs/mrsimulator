#!/bin/bash
#
# docker run -it -v $(pwd):/io quay.io/pypa/manylinux2010_x86_64 /io/build_wheel_linux.sh
#

set -e -x
cd io

# used with docker image - quay.io/pypa/manylinux2014_x86_64
# Install a system package required by our library
yum install -y openblas-devel git fftw-devel

# Compile wheels
for PYBIN in /opt/python/cp3[6-9]*/bin; do
    "${PYBIN}/python" --version
    "${PYBIN}/pip" install --upgrade pip
    "${PYBIN}/pip" install -r requirements.txt
    "${PYBIN}/python" setup.py develop bdist_wheel -d linuxwheels
    "${PYBIN}/pip" install -r requirements-dev.txt
    "${PYBIN}/pytest"
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
