#!/bin/bash

set -ev

# Environment variables
if [[ $CONDA_BUILD ]] ; then
    echo 'CONDA'
    export TEST_DIR=`pwd`
elif [[ $TRAVIS ]] ; then
    echo 'TRAVIS'
    if [[ -z ${TEST_DIR} ]] ; then
        echo 'Error! The ${TEST_DIR} environment variable must be set when using travis-ci!'
        exit 1
    fi
    mkdir -p "${TEST_DIR}"
    SRC_DIR="$TRAVIS_BUILD_DIR"
else 
    echo 'Unknown environment!'
    exit
fi


# Copy test data
rsync -av "${SRC_DIR}/scripts" "${TEST_DIR}" --exclude='[._]*'
rsync -av "${SRC_DIR}/tests" "${TEST_DIR}" --exclude='[._]*'
rsync -av "${SRC_DIR}/setup.cfg" "${TEST_DIR}"


# Common directories
cd "${TEST_DIR}"
PDB_DIR="${TEST_DIR}/pdb"
BLAST_DB_DIR="${TEST_DIR}/blast/db"
ARCHIVE_DIR="${TEST_DIR}/archive"

mkdir -p "${PDB_DIR}"
mkdir -p "${BLAST_DB_DIR}"
mkdir -p "${ARCHIVE_DIR}"

touch "${BLAST_DB_DIR}/nr.pal"
touch "${BLAST_DB_DIR}/pdbaa.pal"

sed -i "s|^pdb_dir = .*|pdb_dir = $PDB_DIR|" "${TEST_DIR}/tests/travis_config_file.ini"
sed -i "s|^blast_db_dir = .*|blast_db_dir = $BLAST_DB_DIR|" "${TEST_DIR}/tests/travis_config_file.ini"
sed -i "s|^archive_dir = .*|archive_dir = $ARCHIVE_DIR|" "${TEST_DIR}/tests/travis_config_file.ini"


# ====== Database ======
if [[ -z ${TEST_SUITE} || ${TEST_SUITE} == database* ]] ; then

# Download external files
wget -P "${TEST_DIR}" \
    -r --no-parent --reject "index.html*" --cut-dirs=4  \
    http://elaspic.kimlab.org/static/download/current_release/Homo_sapiens_test/
wget -P "${TEST_DIR}/elaspic.kimlab.org" \
    http://elaspic.kimlab.org/static/download/current_release/domain.tsv.gz
wget -P "${TEST_DIR}/elaspic.kimlab.org" \
    http://elaspic.kimlab.org/static/download/current_release/domain_contact.tsv.gz

ls "${TEST_DIR}"


# Configure the database
mysql -u root -e 'drop database if exists travis_test';
mysql -u root -e 'create database travis_test';
mysql -u root -e "grant all on travis_test.* to 'travis'@'127.0.0.1'";
mysql -u root travis_test <<'EOF'
DROP FUNCTION IF EXISTS MUTATION_IN_DOMAIN;
DELIMITER $$
CREATE FUNCTION `MUTATION_IN_DOMAIN`(
  mutation VARCHAR(255),
  domain_def VARCHAR(255)
) RETURNS tinyint(1)
    READS SQL DATA
    DETERMINISTIC
BEGIN
DECLARE domain_start int;
DECLARE domain_end int;
DECLARE mutation_pos int;
DECLARE mutation_in_domain bool;
SET domain_start = CONVERT(SUBSTRING_INDEX(domain_def, ':', 1), unsigned integer);
SET domain_end = CONVERT(SUBSTRING_INDEX(domain_def, ':', -1), unsigned integer);
SET mutation_pos = CONVERT(SUBSTRING(mutation, 2, LENGTH(mutation) - 2), unsigned integer);
SET mutation_in_domain = 
    CASE 
        WHEN mutation_pos = 0 THEN null
        WHEN mutation_pos >= domain_start and mutation_pos <= domain_end THEN true 
        ELSE false
    END;
RETURN mutation_in_domain;
END$$
DELIMITER ;
EOF

mysql -u root travis_test <<'EOF'
DROP FUNCTION IF EXISTS  MUTATION_IN_INTERFACE;
DELIMITER $$
CREATE FUNCTION `MUTATION_IN_INTERFACE`(
  mutation VARCHAR(255),
  interacting_aa TEXT
) RETURNS tinyint(1)
    READS SQL DATA
    DETERMINISTIC
BEGIN
DECLARE mutation_pos int;
DECLARE mutation_in_interface int;
SET mutation_pos = SUBSTRING(mutation, 2, LENGTH(mutation) - 2);
SET mutation_in_interface = FIND_IN_SET(mutation_pos, interacting_aa);
RETURN mutation_in_interface != 0;
END$$
DELIMITER ;
EOF

# Load precalculated data to the database
elaspic_database -c "${TEST_DIR}/tests/travis_config_file.ini" create
elaspic_database -c "${TEST_DIR}/tests/travis_config_file.ini" load_data \
    --data_folder "${TEST_DIR}/elaspic.kimlab.org"

# Remove some rows from the database, so that we have something to calculate
mysql -u root travis_test -e "DELETE FROM provean LIMIT 100";
mysql -u root travis_test -e "DELETE FROM uniprot_domain_model LIMIT 100";
mysql -u root travis_test -e "DELETE FROM uniprot_domain_pair_model LIMIT 100";

fi
