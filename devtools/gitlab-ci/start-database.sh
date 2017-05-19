#!/bin/bash

set -ev

DB_BASEDIR="$(conda env list | grep '*' | awk '{print $3}')"
DB_DATADIR="${DB_BASEDIR}/db"
mkdir -p "${DB_DATADIR}"

# Sanity checks
if [[ -z ${MYSQL_INIT_SCRIPT} ]] ; then
    echo 'Required environment variable not set!'
    exit -1
fi

rm -rf ${DB_DATADIR}
mkdir -p ${DB_BASEDIR} ${DB_DATADIR}

mysqld --initialize --basedir=${DB_BASEDIR} --datadir=${DB_DATADIR} || \
  echo -e "\033[31mLooks like MySQL database '${DB_DATADIR}' has already been initialized!\e[0m"

mysqld --secure-file-priv=/ --user=root --basedir=${DB_BASEDIR} --datadir=${DB_DATADIR} \
  --socket=${DB_BASEDIR}/mysql.sock --port=3306 --default_storage_engine=InnoDB \
  --init-file=${MYSQL_INIT_SCRIPT} &
