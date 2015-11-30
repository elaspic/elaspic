#!/bin/bash

set -e

# Copy test data
rsync -av $SRC_DIR/tests ./ --exclude='[._]*'
rsync -av $SRC_DIR/setup.cfg ./

# Directories
PWD=`pwd`
PDB_DIR="$PWD/pdb"
BLAST_DB_DIR="$PWD/blast/db"
ARCHIVE_DIR="$PWD/elaspic.kimlab.org"

mkdir -p "$PDB_DIR"
mkdir -p "$BLAST_DB_DIR"
mkdir -p "$ARCHIVE_DIR"

touch "$BLAST_DB_DIR/nr.pal"
touch "$BLAST_DB_DIR/pdbaa.pal"

# Update the configuration file
sed -i "s|^pdb_dir = .*|pdb_dir = $PDB_DIR|" ./tests/travis_config_file.ini
sed -i "s|^blast_db_dir = .*|blast_db_dir = $BLAST_DB_DIR|" ./tests/travis_config_file.ini
sed -i "s|^archive_dir = .*|archive_dir = $ARCHIVE_DIR|" ./tests/travis_config_file.ini

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

# Download external files
wget -r --no-parent --reject "index.html*" --cut-dirs=4 \
    http://elaspic.kimlab.org/static/download/current_release/Homo_sapiens_test/
wget -P elaspic.kimlab.org \
    http://elaspic.kimlab.org/static/download/current_release/domain.tsv.gz
wget -P elaspic.kimlab.org \
    http://elaspic.kimlab.org/static/download/current_release/domain_contact.tsv.gz

# Extract precalculated data into archive folder
#tar -xzvf elaspic.kimlab.org/provean/provean.tar.gz -C archive/
#tar -xzvf elaspic.kimlab.org/uniprot_domain/uniprot_domain.tar.gz -C archive/
#tar -xzvf elaspic.kimlab.org/uniprot_domain_pair/uniprot_domain_pair.tar.gz -C archive/

# Load precalculated data to the database
elaspic_database -c tests/travis_config_file.ini create
elaspic_database -c tests/travis_config_file.ini load_data --data_folder elaspic.kimlab.org

# Remove some rows from the database
mysql -u root travis_test -e "DELETE FROM provean LIMIT 100";
mysql -u root travis_test -e "DELETE FROM uniprot_domain_model LIMIT 100";
mysql -u root travis_test -e "DELETE FROM uniprot_domain_pair_model LIMIT 100";

# Run tests
py.test -vsx --cache-clear --quick
py.test tests/test_database_pipeline.py -vsx --cache-clear --quick \
    --config-file="./tests/travis_config_file.ini"

# Clear the database
elaspic_database -c tests/travis_config_file.ini delete

