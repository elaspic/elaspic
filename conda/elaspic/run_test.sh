#!/bin/bash

if [[ "$CI" ]] ; then 
  echo "Going to run tests through continuous integration..."
  exit
fi

conda install -q --yes -n _test fake_provean fake_foldx

export TEST_SUITE='database_2'

. $SRC_DIR/scripts/configure_test.sh
. $SRC_DIR/scripts/run_test.sh
. $SRC_DIR/scripts/cleanup_test.sh
