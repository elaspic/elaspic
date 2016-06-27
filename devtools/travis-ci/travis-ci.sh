#!/bin/bash

set -ev

# Sanity checks
if [[ -z ${TEST_SUITE} || \
      -z ${SCRIPTS_DIR} || \
      -z ${TEST_DIR} ]] ; then
    echo 'Required environment variable not set!'
    exit
fi

if [[ ${TEST_SUITE} = 'test_unittests' ]] ; then
    py.test
    "${SCRIPTS_DIR}/build-docs.sh"
else
    source "${SCRIPTS_DIR}/configure_tests.sh"
    cd "${TEST_DIR}"
    ls
    "${SCRIPTS_DIR}/${TEST_SUITE}.sh"
fi
