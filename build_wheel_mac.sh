#!/bin/bash

source /Users/deepansh/opt/anaconda3/bin/activate

conda activate py36
pip install -r requirements.txt
pip install -r requirements-dev.txt
python setup.py develop
pytest
python setup.py bdist_wheel -d macwheels

conda activate py37
pip install -r requirements.txt
pip install -r requirements-dev.txt
python setup.py develop
pytest
python setup.py bdist_wheel -d macwheels

conda activate py38
pip install -r requirements.txt
pip install -r requirements-dev.txt
python setup.py develop
pytest
python setup.py bdist_wheel -d macwheels

delocate-wheel -w dist macwheels/*.whl
