# -*- coding: utf-8 -*-
"""
Created on Sat Sep 26 13:39:23 2015

@author: strokach

TODO: Modify the export database scripts to only export proteins with <= 2 domains
      and <= 3 interactions.

TODO: Add cases for precalculated mutations (can do many of these :).

"""
import os
import os.path as op
import random
import logging
import pytest
import pandas as pd
from conftest import TESTS_BASE_DIR
from elaspic import conf, helper

# %% Constants
if hasattr(pytest, "config"):
    QUICK = pytest.config.getoption('--quick')
    CONFIG_FILE = (
        pytest.config.getoption('--config-file') or
        op.join(TESTS_BASE_DIR, 'test_database_pipeline.ini')
    )
else:
    QUICK = False
    CONFIG_FILE = op.join(TESTS_BASE_DIR, 'test_database_pipeline.ini')


print('Running quick: {}'.format(QUICK))
print('Config file: {}'.format(CONFIG_FILE))


# %%
conf.read_configuration_file(CONFIG_FILE)
from elaspic import database, database_pipeline


# %%
logger = logging.getLogger(__name__)
configs = conf.Configs()

try:
    base_dir = op.dirname(__file__)
except NameError:
    base_dir = os.getcwd()

db = database.MyDatabase()
engine = db.get_engine()
engine.execute("SET sql_mode = ''")


# %%
test_cases = []


def append_test_cases(df, num=3, num_mutations=3):
    """
    Parameters
    ----------
    df : DataFrame
        Contains the following columns:
          - `uniprot_id`
          - `uniprot_sequence`
          - `interacting_aa` OR `model_domain_def` OR `domain_def`

    """
    if QUICK:
        num = 1
    for i in range(num):
        row_idx = random.randint(0, len(df))
        if df.empty:
            raise Exception('empty dataframe supplied: {}'.format(df))
        row = df.iloc[row_idx]
        uniprot_id = row['uniprot_id']
        uniprot_sequence = row['uniprot_sequence']
        logger.debug('Protein ID: {}'.format(uniprot_id))
        for i in range(num_mutations):
            if 'interacting_aa_1' in row and row['interacting_aa_1'] != None:
                mutation_pos = random.choice([int(x) for x in row['interacting_aa_1'].split(',')])
                logger.debug('Selected interface AA: {}'.format(mutation_pos))
            elif 'model_domain_def' in row and row['model_domain_def'] != None:
                domain_start, domain_end = (int(x) for x in row['model_domain_def'].split(':'))
                mutation_pos = random.randint(domain_start, domain_end)
                logger.debug(
                    'Selected AA: {} falling inside model domain: {}'
                    .format(mutation_pos, row['model_domain_def'])
                )
            elif 'domain_def' in row and ['domain_def'] != None:
                domain_start, domain_end = (int(x) for x in row['domain_def'].split(':'))
                mutation_pos = random.randint(domain_start, domain_end)
                logger.debug(
                    'Selected AA: {} falling inside domain: {}'
                    .format(mutation_pos, row['domain_def'])
                )
            else:
                mutation_pos = random.randint(1, len(uniprot_sequence))
            mutation_from = uniprot_sequence[mutation_pos-1]
            mutation_to = random.choice('GVALICMFWPDESTYQNKRH')
            mutation = '{}{}{}'.format(mutation_from, mutation_pos, mutation_to)
            test_cases.append((uniprot_id, mutation,))


# %% Everything is missing
sql_query = """
select ud.uniprot_id, us.uniprot_sequence, udt.domain_def
from {db_schema}.uniprot_domain ud
join {db_schema}.uniprot_domain_template udt using (uniprot_domain_id)
join {db_schema}.uniprot_domain_pair udp on (udp.uniprot_domain_id_1 = ud.uniprot_domain_id)
join {db_schema}.uniprot_domain_pair_template udpt using (uniprot_domain_pair_id)
join {db_schema_uniprot}.uniprot_sequence us using (uniprot_id)
where uniprot_id not in
    (select uniprot_id from {db_schema}.provean)
and uniprot_domain_id not in
    (select uniprot_domain_id from {db_schema}.uniprot_domain_model)
and uniprot_domain_pair_id not in
    (select uniprot_domain_pair_id from {db_schema}.uniprot_domain_pair_model)
and CHAR_LENGTH(us.uniprot_sequence) < 1000
and db = 'sp'
limit 1000;
""".format(db_schema=configs['db_schema'], db_schema_uniprot=configs['db_schema_uniprot'])
df = pd.read_sql_query(sql_query, engine)

logger.debug("Everything is missing: {}".format(len(df)))
if df.empty:
    logger.error("Skipping...")
else:
    append_test_cases(df)


# %% Have provean and domain model but not interface model
sql_query = """
select ud.uniprot_id, us.uniprot_sequence, udm.model_domain_def
from {db_schema}.uniprot_domain ud
join {db_schema}.provean using (uniprot_id)
join {db_schema}.uniprot_domain_template using (uniprot_domain_id)
join {db_schema}.uniprot_domain_model udm using (uniprot_domain_id)
join {db_schema}.uniprot_domain_pair udp on (udp.uniprot_domain_id_1 = ud.uniprot_domain_id)
join {db_schema}.uniprot_domain_pair_template using (uniprot_domain_pair_id)
join {db_schema_uniprot}.uniprot_sequence us using (uniprot_id)
where uniprot_domain_pair_id not in
    (select uniprot_domain_pair_id from {db_schema}.uniprot_domain_pair_model)
and CHAR_LENGTH(us.uniprot_sequence) < 1000
and db = 'sp'
limit 1000;
""".format(db_schema=configs['db_schema'], db_schema_uniprot=configs['db_schema_uniprot'])
df = pd.read_sql_query(sql_query, engine)

logger.debug("Have provean and domain model but not interface model: {}".format(len(df)))
if df.empty:
    logger.error("Skipping...")
elif not QUICK:
    append_test_cases(df)


# %% Have provean and all models
sql_query = """
select ud.uniprot_id, us.uniprot_sequence, udpm.interacting_aa_1
from {db_schema}.uniprot_domain ud
join {db_schema}.provean using (uniprot_id)
join {db_schema}.uniprot_domain_template using (uniprot_domain_id)
join {db_schema}.uniprot_domain_model using (uniprot_domain_id)
join {db_schema}.uniprot_domain_pair udp on (udp.uniprot_domain_id_1 = ud.uniprot_domain_id)
join {db_schema}.uniprot_domain_pair_template using (uniprot_domain_pair_id)
join {db_schema}.uniprot_domain_pair_model udpm using (uniprot_domain_pair_id)
join {db_schema_uniprot}.uniprot_sequence us using (uniprot_id)
where CHAR_LENGTH(us.uniprot_sequence) < 1000
and db = 'sp'
and udpm.model_filename is not null
and udpm.model_errors is null
limit 1000;
""".format(db_schema=configs['db_schema'], db_schema_uniprot=configs['db_schema_uniprot'])
df = pd.read_sql_query(sql_query, engine)

logger.debug("Have provean and all models: {}".format(len(df)))
if df.empty:
    logger.error("Skipping...")
else:
    append_test_cases(df, 10)


# %% Fixtures
@pytest.fixture(scope='session', params=test_cases)
def uniprot_id_mutation(request):
    return request.param


# %%
@database.retry_database
def validate_mutation_1(uniprot_id, mutation):
    """Select Provean; assert length > 0
    """
    logger.debug(helper.underline('Validating that we have provean...'))
    sql_query = """\
select *
from {db_schema}.provean
where uniprot_id = '{uniprot_id}' and
provean_supset_filename is not null;
""".format(uniprot_id=uniprot_id, db_schema=configs['db_schema'])
    logger.debug(sql_query)
    df = pd.read_sql_query(sql_query, engine)
    logger.debug(df)
    assert len(df) >= 1


@database.retry_database
def validate_mutation_2(uniprot_id, mutation):
    """Select domains without models; assert length 0
    """
    logger.debug(helper.underline('Validating that we have domain models...'))
    sql_query = """\
select *
from {db_schema}.uniprot_domain
join {db_schema}.uniprot_domain_template using (uniprot_domain_id)
left join {db_schema}.uniprot_domain_model using (uniprot_domain_id)
where uniprot_id = '{uniprot_id}' and
model_filename is null and
model_errors is null;
""".format(uniprot_id=uniprot_id, db_schema=configs['db_schema'])
    logger.debug(sql_query)
    df = pd.read_sql_query(sql_query, engine)
    assert len(df) == 0


@database.retry_database
def validate_mutation_3(uniprot_id, mutation):
    """Select interfaces without models; assert length 0
    """
    logger.debug(helper.underline('Validating that we have domain pair models...'))
    sql_query = """\
select *
from {db_schema}.uniprot_domain_pair udp
join {db_schema}.uniprot_domain ud1 on (ud1.uniprot_domain_id = udp.uniprot_domain_id_1)
join {db_schema}.uniprot_domain ud2 on (ud2.uniprot_domain_id = udp.uniprot_domain_id_2)
join {db_schema}.uniprot_domain_pair_template udpt using (uniprot_domain_pair_id)
left join {db_schema}.uniprot_domain_pair_model udpm using (uniprot_domain_pair_id)
where (ud1.uniprot_id = '{uniprot_id}' or ud2.uniprot_id = '{uniprot_id}') and
model_filename is null and model_errors is null;
""".format(uniprot_id=uniprot_id, db_schema=configs['db_schema'])
    logger.debug(sql_query)
    df = pd.read_sql_query(sql_query, engine)
    assert len(df) == 0


@database.retry_database
def validate_mutation_4(uniprot_id, mutation):
    """Select domains where we don't have mutatons even though we should; assert length 0
    """
    logger.debug(helper.underline('Validating that we have domain mutations...'))
    sql_query = """\
select *
from {db_schema}.uniprot_domain ud
join {db_schema}.uniprot_domain_template using (uniprot_domain_id)
join {db_schema}.uniprot_domain_model udm using (uniprot_domain_id)
left join {db_schema}.uniprot_domain_mutation udmut on
    (udmut.uniprot_domain_id = ud.uniprot_domain_id and mutation = '{mutation}')
where ud.uniprot_id = '{uniprot_id}' and
    {db_schema}.mutation_in_domain('{mutation}', model_domain_def)
and model_filename_wt is null;
""".format(uniprot_id=uniprot_id, mutation=mutation, db_schema=configs['db_schema'])
    logger.debug(sql_query)
    df = pd.read_sql_query(sql_query, engine)
    assert len(df) == 0


@database.retry_database
def validate_mutation_5(uniprot_id, mutation):
    """Select domain pairs where we don't have mutatons even though we should; assert length 0
    """
    logger.debug(helper.underline('Validating that we have domain pair mutations...'))
    sql_query = """\
select *
from {db_schema}.uniprot_domain_pair udp
join {db_schema}.uniprot_domain ud1 on (ud1.uniprot_domain_id = udp.uniprot_domain_id_1)
join {db_schema}.uniprot_domain ud2 on (ud2.uniprot_domain_id = udp.uniprot_domain_id_2)
join {db_schema}.uniprot_domain_pair_template udpt using (uniprot_domain_pair_id)
join {db_schema}.uniprot_domain_pair_model udpm using (uniprot_domain_pair_id)
left join {db_schema}.uniprot_domain_pair_mutation udpmut on
    (udpmut.uniprot_domain_pair_id = udp.uniprot_domain_pair_id and
     udpmut.uniprot_id = '{uniprot_id}' and udpmut.mutation = '{mutation}')
where (ud1.uniprot_id = '{uniprot_id}' and
    {db_schema}.mutation_in_interface('{mutation}', udpm.interacting_aa_1))
or (ud2.uniprot_id = '{uniprot_id}' and
    {db_schema}.mutation_in_interface('{mutation}', udpm.interacting_aa_2))
and model_filename_wt is null;
""".format(uniprot_id=uniprot_id, mutation=mutation, db_schema=configs['db_schema'])
    logger.debug(sql_query)
    df = pd.read_sql_query(sql_query, engine)
    assert len(df) == 0


# %%
def test_database_pipeline(uniprot_id_mutation):
    #
    if len(uniprot_id_mutation) == 2:
        uniprot_id, mutation = uniprot_id_mutation
        uniprot_domain_pair_ids = []
    else:
        uniprot_id, mutation, uniprot_domain_pair_ids = uniprot_id_mutation
    #
    logger.debug('uniprot_id: {}'.format(uniprot_id))
    logger.debug('mutation: {}'.format(mutation))
    logger.debug('uniprot_domain_pair_ids: {}'.format(uniprot_domain_pair_ids))
    #
    sp = database_pipeline.DatabasePipeline(
        uniprot_id, mutation,
        run_type='5', uniprot_domain_pair_ids=uniprot_domain_pair_ids)
    sp.run()
    logger.debug('Finished running database pipeline. Now checking results...')
    validate_mutation_1(uniprot_id, mutation)
    validate_mutation_2(uniprot_id, mutation)
    validate_mutation_3(uniprot_id, mutation)
    validate_mutation_4(uniprot_id, mutation)
    validate_mutation_5(uniprot_id, mutation)


# %%
problematic_inputs = [
    ('P15153', 'P34S'),
    ('P15891', 'I305H'),  # mutation outside domain
    ('A1VXZ7', 'I92W'),  # mutation outside interface
    ('P24941', 'V197N'),  # chains not interacting
    ('O00522', 'K343Y'),  # model is already calculated; should not calculate again
] + [
    ('P0A910', m)
    for m in 'F144A,F191A,F72A,W164A,W28A,W36A,W78A,Y150A,Y162A,Y189A,Y64A,Y76A'.split(',')
]


# %%
if __name__ == '__main__':
    # import pytest
    # pytest.main([__file__, '-svx', '--quick'])

    # %%
    test_database_pipeline(problematic_inputs[0])

    # %%
    # uniprot_id_mutation = ('A1VXZ7', 'I92W')
    # test_database_pipeline(uniprot_id_mutation)


    # %%
#    uniprot_ids_mutations = [
##        ('P0A910', m)
##        for m in 'F144A,F191A,F72A,W164A,W28A,W36A,W78A,Y150A,Y162A,Y189A,Y64A,Y76A'.split(',')
##        ('P15891', 'I305H',),  # mutation outside domain
##        ('A1VXZ7', 'I92W',),  # mutation outside interface
##        ('P24941', 'V197N'),  # chains not interacting
##        ('O00522', 'K343Y'),  # model is already calculated; should not calculate again
#    ]
#    for uim in uniprot_ids_mutations:
#        test_database_pipeline(uim)
