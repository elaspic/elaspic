# -*- coding: utf-8 -*-
"""
Created on Fri Feb 27 16:57:30 2015

@author: Alexey Strokach
"""
#%%
from __future__ import print_function
import random
import pandas as pd
import sqlalchemy as sa
from init_tests import Pipeline, default_config_file

# Randomly chose the result so that we don't always check on the same domains
result_index = random.randint(0, 100)
print('Result index: {}\n'.format(result_index))

engine = sa.create_engine('mysql://elaspic:elaspic@192.168.6.19/elaspic')



#%%
class TestUniprotDomain:

    domain_query = """
    select
    udmut.uniprot_domain_id,
    udmut.uniprot_id,
    udmut.mutation
    from elaspic.uniprot_domain ud
    join elaspic.uniprot_domain_template using (uniprot_domain_id)
    join elaspic.uniprot_domain_model using (uniprot_domain_id)
    join elaspic.uniprot_domain_mutation udmut using (uniprot_domain_id)
    join uniprot_kb.uniprot_sequence us on (ud.uniprot_id = us.uniprot_id)
    where organism_name = 'Homo sapiens'
    and ddg is not null
    limit {},1;
    """

    delete_provean_command = """
    delete from elaspic.provean
    where uniprot_id = '{}'
    """

    delete_model_command = """
    delete from elaspic.uniprot_domain_model
    where uniprot_domain_id = {}
    """

    provean_validation_query = """
    select *
    from elaspic.provean
    where uniprot_id = '{}'
    """

    model_validation_query = """
    select *
    from elaspic.uniprot_domain_model
    where uniprot_domain_id = {}
    """

    mutation_validation_query = """
    select *
    from elaspic.uniprot_domain_mutation
    where uniprot_domain_id = {}
    """

    def setup_class(self):
        self.pipeline = Pipeline(default_config_file)
        self.domain_for_test = pd.read_sql_query(self.domain_query.format(result_index), engine)
        self.uniprot_domain_id, self.uniprot_id, self.mutation = (
            self.domain_for_test.loc[0, ['uniprot_domain_id', 'uniprot_id', 'mutation']].values
        )
        print(
            'Testing `uniprot_domain`\n',
            'uniprot_domain_id: {}, uniprot_id: {}, mutation: {}'
            .format(self.uniprot_domain_id, self.uniprot_id, self.mutation),
            sep='',
        )

        # Remove precalculated Provean
        engine.execute(self.delete_provean_command.format(self.uniprot_id))
        provean = pd.read_sql_query(
            self.provean_validation_query.format(self.uniprot_id), engine)
        assert(len(provean) == 0)

        # Remove precalculated model and mutation
        engine.execute(self.delete_model_command.format(self.uniprot_domain_id))
        model = pd.read_sql_query(
            self.model_validation_query.format(self.uniprot_domain_id), engine)
        assert(len(model) == 0)
        mutation = pd.read_sql_query(
            self.mutation_validation_query.format(self.uniprot_domain_id), engine)
        assert(len(mutation) == 0)

    def test_t1(self):
        # Calculate new Provean
        self.t1_results = self.pipeline(self.uniprot_id, self.mutation, 1)
        print('t1_results:\n{}'.format(self.t1_results))
        # Make sure that it went successfully
        provean_df = pd.read_sql_query(
            self.provean_validation_query.format(self.uniprot_id), engine)
        print('provean_df:\n{}\n'.format(provean_df))
        assert(len(provean_df) == 1)

    def test_t2(self):
        # Make a new model
        self.t2_results = self.pipeline(self.uniprot_id, self.mutation, 2)
        print('t2_results:\n{}'.format(self.t2_results))
        # Make sure that it went successfully
        model_df = pd.read_sql_query(
            self.model_validation_query.format(self.uniprot_domain_id), engine)
        print('model_df:\n{}\n'.format(model_df))
        assert(len(model_df) == 1)

    def test_t3(self):
        # Calculate mutation results
        self.t3_results = self.pipeline(self.uniprot_id, self.mutation, 3)
        print('t3_results:\n{}'.format(self.t3_results))
        # Make sure that it went successfully
        mutation_df = pd.read_sql_query(
            self.mutation_validation_query.format(self.uniprot_domain_id), engine)
        print('mutation_df:\n{}\n'.format(mutation_df))
        assert(len(mutation_df) == 1)



#%%
class TestUniprotDomainPair:

    domain_pair_query = """
    select
    udpmut.uniprot_domain_pair_id,
    udpmut.uniprot_id,
    udpmut.mutation
    from elaspic.uniprot_domain_pair udp
    join elaspic.uniprot_domain_pair_template using (uniprot_domain_pair_id)
    join elaspic.uniprot_domain_pair_model using (uniprot_domain_pair_id)
    join elaspic.uniprot_domain_pair_mutation udpmut using (uniprot_domain_pair_id)
    join uniprot_kb.uniprot_sequence us on (us.uniprot_id = udpmut.uniprot_id)
    where organism_name = 'Homo sapiens'
    and ddg is not null
    limit {},1;
    """

    delete_model_command = """
    delete from elaspic.uniprot_domain_pair_model
    where uniprot_domain_pair_id = {}
    """

    model_validation_query = """
    select *
    from elaspic.uniprot_domain_pair_model
    where uniprot_domain_pair_id = {}
    """

    mutation_validation_query = """
    select *
    from elaspic.uniprot_domain_pair_mutation
    where uniprot_domain_pair_id = {}
    """

    def setup_class(self):
        # Load data
        self.pipeline = Pipeline(default_config_file)
        self.domain_pair_for_test = pd.read_sql_query(
            self.domain_pair_query.format(result_index), engine)
        self.uniprot_domain_pair_id, self.uniprot_id, self.mutation = (
            self.domain_pair_for_test
                .loc[0, ['uniprot_domain_pair_id', 'uniprot_id', 'mutation']].values
        )
        print(
            'Testing `uniprot_domain pair`\n',
            'uniprot_domain_pair_id: {}, uniprot_id: {}, mutation: {}'
            .format(self.uniprot_domain_pair_id, self.uniprot_id, self.mutation)
        )

        # Remove existing model and mutation
        engine.execute(self.delete_model_command.format(self.uniprot_domain_pair_id))
        model = pd.read_sql_query(
            self.model_validation_query.format(self.uniprot_domain_pair_id), engine)
        assert(len(model) == 0)
        mutation = pd.read_sql_query(
            self.mutation_validation_query.format(self.uniprot_domain_pair_id), engine)
        assert(len(mutation) == 0)

    def test_t2(self):
        # Calculate a new model
        self.t2_results = self.pipeline(self.uniprot_id, self.mutation, 2)
        print('t2_results:\n{}'.format(self.t2_results))
        # Make sure that it went successfully
        model_df = pd.read_sql_query(
            self.model_validation_query.format(self.uniprot_domain_pair_id), engine)
        print('model_df:\n{}\n'.format(model_df))
        assert(len(model_df) == 1)

    def test_t3(self):
        # Calculate new mutation
        self.t3_results = self.pipeline(self.uniprot_id, self.mutation, 3)
        print('t3_results:\n{}'.format(self.t3_results))
        # Make sure that it went successfully
        mutation_df = pd.read_sql_query(
            self.mutation_validation_query.format(self.uniprot_domain_pair_id), engine)
        print('mutation_df:\n{}\n'.format(mutation_df))
        assert(len(mutation_df) == 1)



#%%


