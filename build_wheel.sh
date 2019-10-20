#!/bin/bash
set -e -x

# Install a system package required by our library
apt-get install libopenblas-dev libfftw3-dev
# yum install -y openblas-devel git fftw-devel
mkdir mrsim
git clone --single-branch --branch master https://github.com/DeepanshS/mrsimulator.git mrsim
cd mrsim

# Compile wheels
for PYBIN in /opt/python/cp36*/bin; do
    "${PYBIN}/pip" install -r requirements.txt
    "${PYBIN}/python" setup_unix.py bdist_wheel
done

# Bundle external shared libraries into the wheels
for whl in dist/*.whl; do
    auditwheel repair "$whl" -w /io/wheelhouse/
done


# Install packages and test
for PYBIN in /opt/python/cp36*/bin; do
    "${PYBIN}/pip" install -r requirements-dev.txt
    "${PYBIN}/pip" install mrsimulator --no-index -f /io/wheelhouse
    (cd "$HOME"; "${PYBIN}/pytest mrsim/")
done
