#!/bin/bash

set -ev

# Sanity checks
if [[ -z ${SRC_DIR} || \
      -z ${SCRIPTS_DIR} || \
      -z ${TEST_DIR} ]] ; then
    echo 'Required environment variable not set!'
    exit -1
fi

if [[ ${SRC_DIR} = ${TEST_DIR} ]] ; then
    echo "SRC_DIR and TEST_DIR must be different!"
    exit -1
fi


# Copy package data
mkdir -p "${TEST_DIR}"
rsync -av "${SRC_DIR}" "${TEST_DIR}" --exclude='[._]*'


# Common directories
export PDB_DIR="${TEST_DIR}/pdb"
mkdir -p "${PDB_DIR}"

export BLAST_DB_DIR="${TEST_DIR}/blast/db"
mkdir -p "${BLAST_DB_DIR}"
touch "${BLAST_DB_DIR}/nr.pal"
touch "${BLAST_DB_DIR}/pdbaa.pal"


# ====== Database ======
if [[ -z ${TEST_SUITE} || ${TEST_SUITE} == database* ]] ; then
    ${SCRIPTS_DIR}/configure_database.sh
fi
