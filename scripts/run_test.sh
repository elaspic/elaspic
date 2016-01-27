#!/bin/bash

set -ev


# Sanity checks
if [[ -z ${TEST_DIR} ]] ; then
    echo 'Error! The ${TEST_DIR} environment variable must be set!'
    exit 1
fi


# ====== Local Test 1 ======
if [[ -z ${TEST_SUITE} || ${TEST_SUITE} == 'local_1' ]] ; then

CONFIG_FILE="${TEST_DIR}/local_1.ini"
ARCHIVE_DIR="${TEST_DIR}/archive"
mkdir -p "${ARCHIVE_DIR}"

# Update the configuration file
cp -f "${TEST_DIR}/tests/travis_config_file.ini" "${CONFIG_FILE}"
sed -i "s|^archive_type = .*|archive_type = directory|" "${CONFIG_FILE}"
sed -i "s|^archive_dir = .*|archive_dir = $ARCHIVE_DIR|" "${CONFIG_FILE}"

# Run tests
py.test "${TEST_DIR}" -vsx --cache-clear --config-file="${CONFIG_FILE}" 

fi


# ====== Database Test 1 ======
if [[ -z ${TEST_SUITE} || ${TEST_SUITE} == 'database_1' ]] ; then

CONFIG_FILE="$TEST_DIR/database_1.ini"
ARCHIVE_DIR="$TEST_DIR/archive"
mkdir -p "${ARCHIVE_DIR}"

7z x "${TEST_DIR}/elaspic.kimlab.org/provean/provean.7z" -o"${ARCHIVE_DIR}"
7z x "${TEST_DIR}/elaspic.kimlab.org/uniprot_domain/uniprot_domain.7z" -o"${ARCHIVE_DIR}"
7z x "${TEST_DIR}/elaspic.kimlab.org/uniprot_domain_pair/uniprot_domain_pair.7z" -o"${ARCHIVE_DIR}"

# Update the configuration file
cp -f "${TEST_DIR}/tests/travis_config_file.ini" "${CONFIG_FILE}"
sed -i "s|^archive_type = .*|archive_type = directory|" "${CONFIG_FILE}"
sed -i "s|^archive_dir = .*|archive_dir = $ARCHIVE_DIR|" "${CONFIG_FILE}"

# Run tests
py.test "${TEST_DIR}/tests/test_database_pipeline.py" -vsx --cache-clear \
    --config-file="${CONFIG_FILE}"

fi


# ====== Database Test 2 ======
if [[ -z ${TEST_SUITE} || ${TEST_SUITE} == 'database_2' ]] ; then

CONFIG_FILE="$TEST_DIR/database_2.ini"
ARCHIVE_DIR="$TEST_DIR/elaspic.kimlab.org"
mkdir -p "${ARCHIVE_DIR}"

# Update the configuration file
cp -f "${TEST_DIR}/tests/travis_config_file.ini" "${CONFIG_FILE}"
sed -i "s|^archive_type = .*|archive_type = 7zip|" "${CONFIG_FILE}"
sed -i "s|^archive_dir = .*|archive_dir = $ARCHIVE_DIR|" "${CONFIG_FILE}"

# Run tests
py.test "${TEST_DIR}/tests/test_database_pipeline.py" -vsx --cache-clear \
    --config-file="${CONFIG_FILE}"

fi
