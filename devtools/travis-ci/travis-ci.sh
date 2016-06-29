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
elif [[ ${TEST_SUITE} = 'test_standalone*' ]] ; then
    source "${SCRIPTS_DIR}/configure_tests.sh"
    "${SCRIPTS_DIR}/${TEST_SUITE}.sh"
elif [[ ${TEST_SUITE} = 'test_database*' ]] ; then
    source "${SCRIPTS_DIR}/configure_tests.sh"
    source "${SCRIPTS_DIR}/configure_database.sh"
    "${SCRIPTS_DIR}/${TEST_SUITE}.sh"
else
    echo "Wrong TEST_SUITE: '$TEST_SUITE'"
    exit -1
fi
