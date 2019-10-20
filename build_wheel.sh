#!/bin/bash
set -e -x

# Install a system package required by our library
yum install -y openblas-devel git fftw-devel
mkdir mrsim
git clone --single-branch --branch master https://github.com/DeepanshS/mrsimulator.git mrsim
cd mrsim

# Compile wheels
for PYBIN in /opt/python/cp3[6-8]*/bin; do
    "${PYBIN}/pip" install -r requirements.txt
    "${PYBIN}/python" setup_unix.py bdist_wheel /io/ -w wheelhouse/
done

# Bundle external shared libraries into the wheels
for whl in wheelhouse/*.whl; do
    auditwheel repair "$whl" --plat $PLAT -w /io/wheelhouse/
done

# Install packages and test
for PYBIN in /opt/python/*/bin/; do
    "${PYBIN}/pip" install python-manylinux-demo --no-index -f /io/wheelhouse
    (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
done
