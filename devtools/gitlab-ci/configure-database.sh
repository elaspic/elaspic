#!/bin/bash

set -ev

# Sanity checks
if [[ -z ${CI_PROJECT_DIR} ]] ; then
    echo 'Required environment variables have not been set!'
    exit 1
fi

# Download external files
wget --no-verbose -P "${CI_PROJECT_DIR}/tmp" \
    -r --no-parent --reject "index.html*" --cut-dirs=4  \
    http://elaspic.kimlab.org/static/download/2015-11-01/Homo_sapiens_test/
wget --no-verbose -P "${CI_PROJECT_DIR}/tmp/elaspic.kimlab.org" \
    http://elaspic.kimlab.org/static/download/2015-11-01/domain.tsv.gz
wget --no-verbose -P "${CI_PROJECT_DIR}/tmp/elaspic.kimlab.org" \
    http://elaspic.kimlab.org/static/download/2015-11-01/domain_contact.tsv.gz

ls "${CI_PROJECT_DIR}"


# Configure the database
mysql -u root -Prootpass -e 'drop database if exists elaspic_test';
mysql -u root -Prootpass -e 'create database elaspic_test';
mysql -u root -Prootpass -e "grant all on elaspic_test.* to 'root'@'127.0.0.1'";
mysql -u root -Prootpass -e "set global max_allowed_packet=67108864;"
mysql -u root -Prootpass elaspic_test <<'EOF'
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

mysql -u root -Prootpass elaspic_test <<'EOF'
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
elaspic database -c "${CI_PROJECT_DIR}/test_database_pipeline.ini" create
elaspic database -c "${CI_PROJECT_DIR}/test_database_pipeline.ini" load_data \
    --data_folder "${CI_PROJECT_DIR}/tmp/elaspic.kimlab.org"


# Remove some rows from the database, so that we have something to calculate in our tests
mysql -u root -Prootpass elaspic_test -e "DELETE FROM provean LIMIT 20";
mysql -u root -Prootpass elaspic_test -e "DELETE FROM uniprot_domain_model LIMIT 20";
mysql -u root -Prootpass elaspic_test -e "DELETE FROM uniprot_domain_pair_model LIMIT 20";
