"""
ELASPIC configs should be initialized with a configuration file
before sourcing this file, and should have an `engine` attribute.

The `logger` should also be preconfigured if you want to see something.
"""
import os.path as op
import logging
import shutil
import pandas as pd

from elaspic import (
    conf, helper, elaspic_sequence, structure_tools, standalone_pipeline
)

logger = logging.getLogger(__name__)


# %% Functions
PDB_URL = 'http://www.rcsb.org/pdb/files/{}.pdb'


def get_structure(pdb_id, input_folder, output_folder, use_remote=True):
    """Move PDB structure to the local working directory."""
    input_file = op.join(input_folder, pdb_id + '.pdb')
    output_file = op.join(output_folder, pdb_id + '.pdb')

    # If the PDB already exists, do nothing...
    if op.isfile(output_file):
        logger.debug('Structure file {} already exists!'.format(output_file))
        return output_file

    # Look for PDB file in the same folder
    if not op.isfile(input_file):
        if use_remote:
            input_file = structure_tools.download_pdb_file(pdb_id, input_folder)
        else:
            raise Exception('No PDB input file found!')

    logger.info('Copying {} to {}...'.format(input_file, output_file))
    shutil.copy(input_file, output_file)
    return output_file


def get_sequence(uniprot_id, input_dir, output_dir, use_remote=True):
    """Move PDB structure to the local working directory."""
    input_file = op.join(input_dir, uniprot_id + '.fasta')
    output_file = op.join(output_dir, uniprot_id + '.fasta')

    # If the PDB already exists, do nothing...
    if op.isfile(output_file):
        logger.debug('Sequence file {} already exists!'.format(output_file))
        return output_file

    # Look for PDB file in the same folder
    if not op.isfile(input_file):
        if use_remote:
            input_file = elaspic_sequence.download_uniport_sequence(uniprot_id, input_dir)
        else:
            raise Exception('No PDB input file found!')

    logger.info('Copying {} to {}...'.format(input_file, output_file))
    shutil.copy(input_file, output_file)
    return output_file


# %% Local tests
def run_pdb_mutation_pipeline(pdb_id, pdb_mutatations):
    pdb_file = structure_tools.download_pdb_file(pdb_id, conf.CONFIGS['unique_temp_dir'])
    for chain_id in pdb_mutatations[pdb_id]:
        for mutation in pdb_mutatations[pdb_id][chain_id]:
            mutation_pdb = '{}_{}'.format(chain_id, mutation)
            lp = standalone_pipeline.StandalonePipeline(
                pdb_file, mutations=mutation_pdb)
            lp.run_all_sequences()
            lp.run_all_models()
            lp.run_all_mutations()
    logger.info('Pass')


def run_sequence_mutation_pipeline(pdb_id_sequence, sequence_mutations):
    """
    Parameters
    ----------
    pdb_id_sequence : tuple
        (pdb_id, uniprot_id)
    """
    pdb_id, sequence_id = pdb_id_sequence
    pdb_file = structure_tools.download_pdb_file(pdb_id, conf.CONFIGS['unique_temp_dir'])
    sequence_file = (
        elaspic_sequence.download_uniport_sequence(sequence_id, conf.CONFIGS['unique_temp_dir'])
    )
    for chain_pos in sequence_mutations[pdb_id_sequence]:
        for mutation in sequence_mutations[pdb_id_sequence][chain_pos]:
            mutation_sequence = '{}_{}'.format(chain_pos, mutation)
            lp = standalone_pipeline.StandalonePipeline(
                pdb_file, sequence_file, mutations=mutation_sequence)
            lp.run_all_sequences()
            lp.run_all_models()
            lp.run_all_mutations()
    logger.info('Pass')


# %% Database tests
@helper.retry_database
def validate_mutation_1(uniprot_id, mutation):
    """Select Provean; assert length > 0
    """
    logger.debug(helper.underline('Validating that we have provean...'))
    sql_query = """\
select 1
from {db_schema}.provean
where uniprot_id = '{uniprot_id}' and
provean_supset_filename is not null;
""".format(uniprot_id=uniprot_id, db_schema=conf.CONFIGS['db_schema'])
    logger.debug(sql_query)
    df1 = pd.read_sql_query(sql_query, conf.CONFIGS['engine'])
    logger.debug(df1.head(2))
    #
    logger.debug(helper.underline('And that we have at least one domain with a template...'))
    sql_query = """\
select 1
from {db_schema}.uniprot_domain
join {db_schema}.uniprot_domain_template using (uniprot_domain_id)
where uniprot_id = '{uniprot_id}';
""".format(uniprot_id=uniprot_id, db_schema=conf.CONFIGS['db_schema'])
    logger.debug(sql_query)
    df2 = pd.read_sql_query(sql_query, conf.CONFIGS['engine'])
    logger.debug(df2.head(2))
    assert len(df1) >= 1 or len(df2) == 0


@helper.retry_database
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
""".format(uniprot_id=uniprot_id, db_schema=conf.CONFIGS['db_schema'])
    logger.debug(sql_query)
    df = pd.read_sql_query(sql_query, conf.CONFIGS['engine'])
    assert len(df) == 0


@helper.retry_database
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
""".format(uniprot_id=uniprot_id, db_schema=conf.CONFIGS['db_schema'])
    logger.debug(sql_query)
    df = pd.read_sql_query(sql_query, conf.CONFIGS['engine'])
    assert len(df) == 0


@helper.retry_database
def validate_mutation_4(uniprot_id, mutation):
    """Select domains where we don't have mutatons even though we should; assert length 0
    """
    logger.debug(helper.underline('Validating that we have domain mutations...'))
    sql_query = """\
SELECT *
FROM {db_schema}.uniprot_domain ud
JOIN {db_schema}.uniprot_domain_template using (uniprot_domain_id)
JOIN {db_schema}.uniprot_domain_model udm using (uniprot_domain_id)
LEFT JOIN {db_schema}.uniprot_domain_mutation udmut ON
    (udmut.uniprot_domain_id = ud.uniprot_domain_id AND mutation = '{mutation}')
WHERE
ud.uniprot_id = '{uniprot_id}' AND
{db_schema}.mutation_in_domain('{mutation}', model_domain_def) AND
udm.model_errors is null AND  -- we allow for some model errors
model_filename_wt is null;
""".format(uniprot_id=uniprot_id, mutation=mutation, db_schema=conf.CONFIGS['db_schema'])
    logger.debug(sql_query)
    df = pd.read_sql_query(sql_query, conf.CONFIGS['engine'])
    assert len(df) == 0


@helper.retry_database
def validate_mutation_5(uniprot_id, mutation):
    """Select domain pairs where we don't have mutatons even though we should; assert length 0
    """
    logger.debug(helper.underline('Validating that we have domain pair mutations...'))
    sql_query = """\
SELECT *
FROM {db_schema}.uniprot_domain_pair udp
JOIN {db_schema}.uniprot_domain ud1 on (ud1.uniprot_domain_id = udp.uniprot_domain_id_1)
JOIN {db_schema}.uniprot_domain ud2 on (ud2.uniprot_domain_id = udp.uniprot_domain_id_2)
JOIN {db_schema}.uniprot_domain_pair_template udpt using (uniprot_domain_pair_id)
JOIN {db_schema}.uniprot_domain_pair_model udpm using (uniprot_domain_pair_id)
LEFT JOIN {db_schema}.uniprot_domain_pair_mutation udpmut ON
    (udpmut.uniprot_domain_pair_id = udp.uniprot_domain_pair_id AND
     udpmut.uniprot_id = '{uniprot_id}' and udpmut.mutation = '{mutation}')
WHERE
((ud1.uniprot_id = '{uniprot_id}' AND
    {db_schema}.mutation_in_interface('{mutation}', udpm.interacting_aa_1))
OR (ud2.uniprot_id = '{uniprot_id}' AND
    {db_schema}.mutation_in_interface('{mutation}', udpm.interacting_aa_2))) AND
udpm.model_errors IS NULL AND  -- we allow for some model errors
model_filename_wt IS NULL;
""".format(uniprot_id=uniprot_id, mutation=mutation, db_schema=conf.CONFIGS['db_schema'])
    logger.debug(sql_query)
    df = pd.read_sql_query(sql_query, conf.CONFIGS['engine'])
    assert len(df) == 0


def run_database_pipeline(uniprot_id_mutation):
    from elaspic import database_pipeline
    if 'engine' not in conf.CONFIGS:
        raise Exception('You must add an `engine` to the configs to use this function!')
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
    logger.info('Pass')
