#!/bin/bash

if [[ -z $ANACONDA_TOKEN || \
      -z $KEY_MODELLER || \
      -z $1 ]] ; then
  echo "Required environment variables have not been set!"
  exit -1
fi

mkdir -p conda-bld/linux-64 conda-bld/noarch
chmod ugo+rwX -R conda-bld

gitlab-runner exec docker \
    --env ANACONDA_TOKEN="${ANACONDA_TOKEN}" \
    --env KEY_MODELLER="${KEY_MODELLER}" \
    --env CI_PROJECT_NAME=elaspic \
    --docker-volumes `pwd`/conda-bld/linux-64:/opt/conda/conda-bld/linux-64 \
    $1
