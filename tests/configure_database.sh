#!/bin/bash

set -ev

# Sanity checks
if [[ -z ${SCRIPTS_DIR} ||
      -z ${TEST_DIR} ]] ; then
    echo 'Required environment variables have not been set!'
    exit 1
fi

# Download external files
wget --no-verbose -P "${TEST_DIR}" \
    -r --no-parent --reject "index.html*" --cut-dirs=4  \
    http://elaspic.kimlab.org/static/download/2015-11-01/Homo_sapiens_test/
wget --no-verbose -P "${TEST_DIR}/elaspic.kimlab.org" \
    http://elaspic.kimlab.org/static/download/2015-11-01/domain.tsv.gz
wget --no-verbose -P "${TEST_DIR}/elaspic.kimlab.org" \
    http://elaspic.kimlab.org/static/download/2015-11-01/domain_contact.tsv.gz

ls "${TEST_DIR}"


# Configure the database
mysql -u root -e 'drop database if exists travis_test';
mysql -u root -e 'create database travis_test';
mysql -u root -e "grant all on travis_test.* to 'travis'@'127.0.0.1'";
mysql -u root -e "set global max_allowed_packet=67108864;"
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
elaspic database -c "${SCRIPTS_DIR}/database_pipeline.ini" create
elaspic database -c "${SCRIPTS_DIR}/database_pipeline.ini" load_data \
    --data_folder "${TEST_DIR}/elaspic.kimlab.org"


# Remove some rows from the database, so that we have something to calculate in our tests
mysql -u root travis_test -e "DELETE FROM provean LIMIT 20";
mysql -u root travis_test -e "DELETE FROM uniprot_domain_model LIMIT 20";
mysql -u root travis_test -e "DELETE FROM uniprot_domain_pair_model LIMIT 20";
