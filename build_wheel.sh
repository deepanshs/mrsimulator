#!/bin/bash
set -e -x

# Install a system package required by our library
yum install -y openblas-devel git fftw-devel
mkdir io
git clone --single-branch --branch master https://github.com/DeepanshS/mrsimulator.git io
cd io

# Compile wheels
for PYBIN in /opt/python/cp3[6-8]*/bin; do
    "${PYBIN}/pip" install -r requirements.txt
    "${PYBIN}/python" setup_unix.py bdist_wheel --dist-dir=/io/dist/
done

# Bundle external shared libraries into the wheels
for whl in dist/*.whl; do
    auditwheel repair "$whl" -w /io/dist/
done


# Install packages and test
# for PYBIN in /opt/python/*/bin/; do
#     "${PYBIN}/pip" install python-manylinux-demo --no-index -f /io/dist
#     (cd "$HOME"; "${PYBIN}/nosetests" pymanylinuxdemo)
# done
