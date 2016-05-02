#!/bin/bash

set -ev

if [[ "$CI" ]] ; then
  echo "Going to run tests through continuous integration..."
  exit 0
fi

conda env list
conda install -n root fake_provean fake_foldx
head $(which provean) -n 2
head $(which foldx) -n 2

export TEST_SUITE='test_database_pipeline_2'

# $SRC_DIR set by conda
export SCRIPTS_DIR="$SRC_DIR/devtools/travis-ci"
if [[ -z $TMP_DIR ]] ; then
    export TEST_DIR="/tmp/elaspic_test"
else
    export TEST_DIR="$TMP_DIR/elaspic_test"
fi

source $SCRIPTS_DIR/configure_tests.sh
cd "${TEST_DIR}" && ls

echo ${SRC_DIR}
echo ${SCRIPTS_DIR}
echo ${TEST_DIR}
echo ${PDB_DIR}
echo ${BLAST_DB_DIR}

source $SCRIPTS_DIR/${TEST_SUITE}.sh
source $SCRIPTS_DIR/cleanup_test.sh
