# -*- coding: utf-8 -*-
import os
import six
import random
import re
import sqlalchemy as sa

from elaspic import conf
from elaspic import helper_functions as hf
from elaspic import test_elaspic

if six.PY3:
    from imp import reload


#%%
DEBUG = False

config_filenames = {
    'mysql': 'config_file_mysql_test.ini',
    'postgresql': 'config_file_postgresql_test.ini',
}

try:
    base_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), '../')
except:
    base_path = os.path.join(os.getcwd(), '../')
code_path = os.path.join(base_path, 'elaspic/')

print('base_path: {}'.format(base_path))



#%%

###
mysql_load_table_template = (
    r"""mysql --local-infile --host={db_url} --user={db_username} --password={db_password} """
    r"""{db_schema} -e "{sql_command}" """
)

psql_load_table_template = (
    r"""PGPASSWORD={db_password} psql -h {db_url} -p {db_port} -U {db_username} """
    r"""-d {db_database} -c "{sql_command}" """
)

# Need to double up on '\\'
mysql_command_template = (
    r"""load data local infile '{organism_folder}/{table_name}.tsv' """
    r"""into table {db_schema}.{table_name} """
    r"""fields terminated by '\t' escaped by '\\\\' lines terminated by '\n'; """
)

psql_command_template = (
    r"""\\copy {db_schema}.{table_name} """
    r"""from '{organism_folder}/{table_name}.tsv' """
    r"""with csv delimiter E'\t' null '\N' escape '\\'; """
)

def _format_configs(configs, table_name):
    configs['table_name'] = table_name
    if table_name == 'uniprot_sequence':
        configs['db_schema'] = configs['db_schema_uniprot']
    return configs
    
    
def load_table_to_db(db_type, configs, table_name):
    configs = _format_configs(configs.copy(), table_name)
    if db_type == 'mysql':
        configs['sql_command'] = mysql_command_template.format(**configs)
        system_command = mysql_load_table_template.format(**configs)
    elif db_type == 'postgresql':
        configs['sql_command'] = psql_command_template.format(**configs)
        system_command = psql_load_table_template.format(**configs)
    else:
        raise Exception('Unsupported database type: {}'.format(db_type))
    if DEBUG:
        print(system_command)
    result, error_message, return_code = hf.popen(system_command)
    if return_code != 0:
        print(result)
        raise Exception(error_message)



###
table_names = [
    'domain', 'domain_contact', 
    'uniprot_sequence', 'provean', 
    'uniprot_domain', 'uniprot_domain_template', 'uniprot_domain_model', 'uniprot_domain_mutation',
    'uniprot_domain_pair', 'uniprot_domain_pair_template', 
    'uniprot_domain_pair_model', 'uniprot_domain_pair_mutation',
]

def load_all_tables_to_db(db_type, configs):
    for table_name in table_names:
        print(table_name)
        load_table_to_db(db_type, configs, table_name)
    


###
def delete_schema(configs, logger, clear_schema, keep_uniprot_sequence):
    from elaspic import sql_db
    
    if not clear_schema:
        return
    
    configs = configs.copy()
    my_db = sql_db.MyDatabase(configs, logger)
    
    for table_name in reversed(table_names):
        if table_name == 'uniprot_sequence':
            continue
        configs['table_name'] = table_name
        my_db.engine.execute('drop table {db_schema}.{table_name};'.format(**configs))

    def drop_uniprot_table_and_schema():
        my_db.engine.execute('drop table {db_schema_uniprot}.uniprot_sequence;'.format(**configs))
        my_db.engine.execute('drop schema {db_schema_uniprot};'.format(**configs))
        
    if configs['db_schema'] != configs['db_schema_uniprot']:
        my_db.engine.execute('drop schema {db_schema};'.format(**configs))
        if not keep_uniprot_sequence:
            drop_uniprot_table_and_schema()
    else:
        if not keep_uniprot_sequence:
            drop_uniprot_table_and_schema()



#%%
def test(DB_TYPE='mysql'):   
    clear_schema = True
    keep_uniprot_sequence = False
    config_filename = os.path.join(base_path, 'config/', config_filenames[DB_TYPE])
    logger = hf.get_logger()
    result_index = random.randint(0, 100) 
    organism_folder = os.path.join(base_path, 'database/Homo_sapiens_test')
    
    logger.info('*' * 80)
    logger.info('Reading the configuration file...')
    conf.read_configuration_file(config_filename)
    configs = conf.configs.copy()
    configs['organism_folder'] = organism_folder
    configs['result_index'] = result_index
    
    from elaspic import sql_db
    logger.info('*' * 80)
    logger.info('Creating an empty database schema...')
    reload(sql_db)
    sql_db.create_database_tables(
        configs=configs, clear_schema=clear_schema, keep_uniprot_sequence=keep_uniprot_sequence)

    logger.info('*' * 80)
    logger.info('Loading data from text files to the empty database...')
    load_all_tables_to_db(DB_TYPE, configs)
    
    logger.info('*' * 80)
    logger.info('Running a sample domain mutation...')
    test_uniprot_domain = test_elaspic.TestUniprotDomain()
    test_uniprot_domain.setup_class(configs)
    test_uniprot_domain.test(make_provean=False)
    
    logger.info('*' * 80)
    logger.info('Running a sample interface mutation...')
    test_uniprot_domain_pair = test_elaspic.TestUniprotDomainPair()
    test_uniprot_domain_pair.setup_class(configs)
    test_uniprot_domain_pair.test()
    
    logger.info('*' * 80)
    logger.info('Removing the schema that we created for testing')
    delete_schema(configs, logger, clear_schema, keep_uniprot_sequence)



#%%
if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--db_type', type=str, default='mysql')
    args = parser.parse_args()
    test(args.db_type)
    
    
