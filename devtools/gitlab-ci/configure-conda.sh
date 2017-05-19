#!/bin/bash

set -ev

if [[ -z ${PACKAGE_VERSION} || \
      -z ${ANACONDA_TOKEN} ]] ; then
    echo 'Required environment variable not set!'
    exit -1
fi

conda config --add channels ostrokach-forge
conda config --append channels bioconda
conda config --append channels salilab
case "${PACKAGE_VERSION}" in
  *dev*)
    conda config --append channels kimlab/label/dev;
    conda config --append channels kimlab;
    conda config --append channels ostrokach/label/dev;
    conda config --append channels ostrokach;
  ;;
  *)
    conda config --append channels kimlab;
    conda config --append channels ostrokach;
  ;;
  esac
conda config --append channels https://conda.anaconda.org/t/$ANACONDA_TOKEN/ostrokach
conda install -y -q requests==2.11
