#!/bin/bash

set -ev


# ====== Database ======
if [[ -z ${TEST_SUITE} || ${TEST_SUITE} == database* ]] ; then

# Clear the database
# ..or keep it so you run tests if something goes wrong
# elaspic database -c "${TEST_DIR}/tests/travis_config_file.ini" delete

fi


# Unset environmental variables
unset TEST_DIR
