#!/bin/bash

set -ex

$PYTHON setup.py train
$PYTHON setup.py install
