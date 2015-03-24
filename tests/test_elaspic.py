# -*- coding: utf-8 -*-
from __future__ import print_function
from builtins import object
import pandas as pd
from elaspic import pipeline, sql_db


#%%
class TestUniprotDomain(object):

    def setup_class(self, configs):
        """
        Source: http://pytest.org/latest/xunit_setup.html?highlight=class#class-level-setup-teardown
        
        Parameters
        ----------
        configs : dict
            Has the following keys: 
            - db_schema : database schema containing elaspic data
            - db_schema_uniprot : database schema containing uniprot sequences
            - result_index : index of the row which will be deleted and re-analysed
        """
        self.configs = configs.copy()
        self.engine = sql_db.get_engine(configs)
        self.pipeline = pipeline.Pipeline(configs)
        
        domain_query = """
        select
        udmut.uniprot_domain_id,
        udmut.uniprot_id,
        udmut.mutation
        from {db_schema}.uniprot_domain ud
        join {db_schema}.uniprot_domain_template using (uniprot_domain_id)
        join {db_schema}.uniprot_domain_model using (uniprot_domain_id)
        join {db_schema}.uniprot_domain_mutation udmut using (uniprot_domain_id)
        join {db_schema_uniprot}.uniprot_sequence us on (ud.uniprot_id = us.uniprot_id)
        where ddg is not null
        limit 1 offset {result_index};
        """.format(**self.configs)
        
        self.domain_for_test = pd.read_sql_query(domain_query, self.engine)
        self.configs['uniprot_domain_id'], self.configs['uniprot_id'], self.configs['mutation'] = (
            self.domain_for_test.loc[0, ['uniprot_domain_id', 'uniprot_id', 'mutation']].values
        )
        
        self.pipeline.logger.info(
            ('\n' + '#' * 80 +
            '\nResult index: {result_index}'
            '\nTesting `uniprot_domain`'
            '\nuniprot_domain_id: {uniprot_domain_id}, uniprot_id: {uniprot_id}, mutation: {mutation}'
            '\n' + '#' * 80)
            .format(**self.configs)
        )
        
        # SQL commands for `test`
        self.delete_provean_command = """
        delete from {db_schema}.provean
        where uniprot_id = '{uniprot_id}'
        """.format(**self.configs)
        
        self.delete_model_command = """
        delete from {db_schema}.uniprot_domain_model
        where uniprot_domain_id = {uniprot_domain_id}
        """.format(**self.configs)
        
        self.provean_validation_query = """
        select *
        from {db_schema}.provean
        where uniprot_id = '{uniprot_id}'
        """.format(**self.configs)
        
        self.model_validation_query = """
        select *
        from {db_schema}.uniprot_domain_model
        where uniprot_domain_id = {uniprot_domain_id}
        """.format(**self.configs)
        
        self.mutation_validation_query = """
        select *
        from {db_schema}.uniprot_domain_mutation
        where uniprot_domain_id = {uniprot_domain_id}
        """.format(**self.configs)
    
    
    def test(self, make_provean=False):
        self.remove_precalculated(make_provean)
        if make_provean:
            self.make_provean()
        self.make_model()
        self.make_mutation()


    def remove_precalculated(self, make_provean):
        if make_provean:
            # Remove precalculated Provean
            self.engine.execute(self.delete_provean_command)
            provean = pd.read_sql_query(self.provean_validation_query, self.engine)
            assert(len(provean) == 0)

        # Remove precalculated model and mutation
        self.engine.execute(self.delete_model_command)
        model = pd.read_sql_query(self.model_validation_query, self.engine)
        assert(len(model) == 0)
        mutation = pd.read_sql_query(self.mutation_validation_query, self.engine)
        assert(len(mutation) == 0)
        print('Successfully removed provean, mutation and model...')


    def make_provean(self):
        # Calculate new Provean
        self.t1_results = self.pipeline(self.uniprot_id, self.mutation, 1)
        print('t1_results:\n{}'.format(self.t1_results))
        # Make sure that it went successfully
        provean_df = pd.read_sql_query(self.provean_validation_query, self.engine)
        print('provean_df:\n{}\n'.format(provean_df))
        assert(len(provean_df) == 1)


    def make_model(self):
        # Make a new model
        self.t2_results = self.pipeline(self.configs['uniprot_id'], self.configs['mutation'], 2)
        print('t2_results:\n{}'.format(self.t2_results))
        # Make sure that it went successfully
        model_df = pd.read_sql_query(self.model_validation_query, self.engine)
        print('model_df:\n{}\n'.format(model_df))
        assert(len(model_df) == 1)


    def make_mutation(self):
        # Calculate mutation results
        self.t3_results = self.pipeline(self.configs['uniprot_id'], self.configs['mutation'], 3)
        print('t3_results:\n{}'.format(self.t3_results))
        # Make sure that it went successfully
        mutation_df = pd.read_sql_query(self.mutation_validation_query, self.engine)
        print('mutation_df:\n{}\n'.format(mutation_df))
        assert(len(mutation_df) == 1)



#%%
class TestUniprotDomainPair(object):

    def setup_class(self, configs):
        """
        Source: http://pytest.org/latest/xunit_setup.html?highlight=class#class-level-setup-teardown
        
        Parameters
        ----------
        configs : dict
            Has the following keys: 
            - db_schema : database schema containing elaspic data
            - db_schema_uniprot : database schema containing uniprot sequences
            - result_index : index of the row which will be deleted and re-analysed
        """
        # Load data
        self.configs = configs.copy()
        self.engine = sql_db.get_engine(self.configs)
        self.pipeline = pipeline.Pipeline(self.configs)
                
        domain_pair_query = """
        select
        udpmut.uniprot_domain_pair_id,
        udpmut.uniprot_id,
        udpmut.mutation
        from {db_schema}.uniprot_domain_pair udp
        join {db_schema}.uniprot_domain_pair_template using (uniprot_domain_pair_id)
        join {db_schema}.uniprot_domain_pair_model using (uniprot_domain_pair_id)
        join {db_schema}.uniprot_domain_pair_mutation udpmut using (uniprot_domain_pair_id)
        join {db_schema_uniprot}.uniprot_sequence us on (us.uniprot_id = udpmut.uniprot_id)
        where ddg is not null
        limit 1 offset {result_index};
        """.format(**self.configs)
        
        self.domain_pair_for_test = pd.read_sql_query(domain_pair_query, self.engine)
        self.configs['uniprot_domain_pair_id'], self.configs['uniprot_id'], self.configs['mutation'] = (
            self.domain_pair_for_test
                .loc[0, ['uniprot_domain_pair_id', 'uniprot_id', 'mutation']].values
        )
        
        self.pipeline.logger.info(
            ('\n' + '#' * 80 +
            '\nResult index: {result_index}'
            '\nTesting `uniprot_domain pair`'
            '\nuniprot_domain_pair_id: {uniprot_domain_pair_id}, uniprot_id: {uniprot_id}, mutation: {mutation}'
            '\n' + '#' * 80)
            .format(**self.configs)
        )
        
        # SQL commands for `test`
        self.delete_model_command = """
        delete from {db_schema}.uniprot_domain_pair_model
        where uniprot_domain_pair_id = {uniprot_domain_pair_id}
        """.format(**self.configs)
    
        self.model_validation_query = """
        select *
        from {db_schema}.uniprot_domain_pair_model
        where uniprot_domain_pair_id = {uniprot_domain_pair_id}
        """.format(**self.configs)
    
        self.mutation_validation_query = """
        select *
        from {db_schema}.uniprot_domain_pair_mutation
        where uniprot_domain_pair_id = {uniprot_domain_pair_id}
        """.format(**self.configs)
    

    def test(self):
        self.remove_precalculated()
        self.make_model()
        self.make_mutation()

    def remove_precalculated(self):
        # Remove existing model and mutation
        self.engine.execute(self.delete_model_command)
        model = pd.read_sql_query(self.model_validation_query, self.engine)
        assert(len(model) == 0)
        mutation = pd.read_sql_query(self.mutation_validation_query, self.engine)
        assert(len(mutation) == 0)

    def make_model(self):
        # Calculate a new model
        self.t2_results = self.pipeline(self.configs['uniprot_id'], self.configs['mutation'], 2)
        print('t2_results:\n{}'.format(self.t2_results))
        # Make sure that it went successfully
        model_df = pd.read_sql_query(self.model_validation_query, self.engine)
        print('model_df:\n{}\n'.format(model_df))
        assert(len(model_df) == 1)

    def make_mutation(self):
        # Calculate new mutation
        self.t3_results = self.pipeline(self.configs['uniprot_id'], self.configs['mutation'], 3)
        print('t3_results:\n{}'.format(self.t3_results))
        # Make sure that it went successfully
        mutation_df = pd.read_sql_query(self.mutation_validation_query, self.engine)
        print('mutation_df:\n{}\n'.format(mutation_df))
        assert(len(mutation_df) == 1)





