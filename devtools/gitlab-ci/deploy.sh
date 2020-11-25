#!/bin/bash

case "${PACKAGE_VERSION}" in
  *dev*)
    anaconda \
      -t $ANACONDA_TOKEN upload $CI_PROJECT_DIR/conda-bld/linux-64/*.tar.bz2 \
      -u kimlab \
      --label dev \
      --force
  ;;
  *)
    anaconda
      -t $ANACONDA_TOKEN upload $CI_PROJECT_DIR/conda-bld/linux-64/*.tar.bz2 \
      -u kimlab
  ;;
esac
