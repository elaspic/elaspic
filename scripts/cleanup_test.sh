#!/bin/bash

# Clear the database
elaspic_database -c "$TEST_DIR/tests/travis_config_file.ini" delete

# Unset environmental variables
unset TEST_DIR
