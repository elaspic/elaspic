#!/bin/bash

set -ev

# Sanity checks
if [[ -z ${TEST_DIR} || \
      -z ${PDB_DIR} || \
      -z ${BLAST_DB_DIR} ]] ; then
    echo 'Required environment variable not set!'
    exit -1
fi

# Run tests
CONFIG_FILE="${TEST_DIR}/$(basename ${0%%.sh}).ini"
ARCHIVE_DIR="${TEST_DIR}/elaspic.kimlab.org"
mkdir -p "${ARCHIVE_DIR}"

# Update the configuration file
cp -f "${TEST_DIR}/database_pipeline.ini" "${CONFIG_FILE}"
sed -i "s|^pdb_dir = .*|pdb_dir = ${PDB_DIR}|" "${CONFIG_FILE}"
sed -i "s|^blast_db_dir = .*|blast_db_dir = ${BLAST_DB_DIR}|" "${CONFIG_FILE}"
sed -i "s|^archive_dir = .*|archive_dir = ${ARCHIVE_DIR}|" "${CONFIG_FILE}"
sed -i "s|^archive_type = .*|archive_type = 7zip|" "${CONFIG_FILE}"

# Run tests
py.test "${TEST_DIR}/test_database_pipeline.py" -vsx --cache-clear --cov=elaspic \
    --config-file="${CONFIG_FILE}"
