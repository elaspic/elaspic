#!/bin/bash


cd "${TEST_DIR}"

# ====== Database ======
if [[ -z ${TEST_SUITE} || ${TEST_SUITE} == database* ]] ; then

# Clear the database
elaspic_database -c "$TEST_DIR/tests/travis_config_file.ini" delete

fi


# Unset environmental variables
unset TEST_DIR
