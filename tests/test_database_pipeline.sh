#!/bin/bash

set -ev

# Sanity checks
if [[ -z ${TEST_DIR} || \
      -z ${PDB_DIR} || \
      -z ${BLAST_DB_DIR} || \
      -z ${TEST_NUMBER} ]] ; then
    echo 'Required environment variable not set!'
    exit -1
fi

if [[ ${TEST_NUMBER} != "1" && ${TEST_NUMBER} != "2" ]] ; then
  echo "Incorrect TEST_NUMBER: '$TEST_NUMBER'"
  exit -1
fi

# Run tests
CONFIG_FILE="${TEST_DIR}/$(basename ${0%%.sh}).ini"

# Update the configuration file
cp -f "${TEST_DIR}/database_pipeline.ini" "${CONFIG_FILE}"
sed -i "s|^pdb_dir = .*|pdb_dir = ${PDB_DIR}|" "${CONFIG_FILE}"
sed -i "s|^blast_db_dir = .*|blast_db_dir = ${BLAST_DB_DIR}|" "${CONFIG_FILE}"

if [[ ${TEST_NUMBER} == "1" ]] ; then
  # Archive directory
  ARCHIVE_DIR="${TEST_DIR}/archive"
  sed -i "s|^archive_dir = .*|archive_dir = ${ARCHIVE_DIR}|" "${CONFIG_FILE}"
  sed -i "s|^archive_type = .*|archive_type = directory|" "${CONFIG_FILE}"
  # Extract archive contents
  7z -bd x "${TEST_DIR}/elaspic.kimlab.org/provean/provean.7z" -o"${ARCHIVE_DIR}"
  7z -bd x "${TEST_DIR}/elaspic.kimlab.org/uniprot_domain/uniprot_domain.7z" -o"${ARCHIVE_DIR}"
  7z -bd x "${TEST_DIR}/elaspic.kimlab.org/uniprot_domain_pair/uniprot_domain_pair.7z" -o"${ARCHIVE_DIR}"
  # Set archive type to directory
else
  ARCHIVE_DIR="${TEST_DIR}/elaspic.kimlab.org"
  sed -i "s|^archive_dir = .*|archive_dir = ${ARCHIVE_DIR}|" "${CONFIG_FILE}"
  sed -i "s|^archive_type = .*|archive_type = 7zip|" "${CONFIG_FILE}"
else

# Run tests
py.test "${TEST_DIR}/test_database_pipeline.py" -vsx --cache-clear --cov=elaspic \
    --config-file="${CONFIG_FILE}"
