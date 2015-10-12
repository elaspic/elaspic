# -*- coding: utf-8 -*-
from __future__ import print_function
from __future__ import unicode_literals
from __future__ import absolute_import
from future import standard_library
standard_library.install_aliases()
from builtins import object

import os
import os.path as op
import json
import signal
import six
import logging 
from functools import wraps

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from . import conf, errors, helper, structure_tools, structure_analysis, sequence, model, predictor

logger = logging.getLogger(__name__)
configs = conf.Configs()

sql_db = None
domain_alignment = None
domain_model = None
domain_mutation = None

MAX_DISTANCE_BETWEEN_INTERACTING_CHAINS = 6 # Angstrom
ELASPIC_LOGO = """

8888888888 888             d8888  .d8888b.  8888888b. 8888888 .d8888b.  
888        888            d88888 d88P  Y88b 888   Y88b  888  d88P  Y88b 
888        888           d88P888 Y88b.      888    888  888  888    888 
8888888    888          d88P 888  "Y888b.   888   d88P  888  888        
888        888         d88P  888     "Y88b. 8888888P"   888  888        
888        888        d88P   888       "888 888         888  888    888 
888        888       d8888888888 Y88b  d88P 888         888  Y88b  d88P 
8888888888 88888888 d88P     888  "Y8888P"  888       8888888 "Y8888P"  

"""



#%% 
class Pipeline(object):

    def __init__(self, configurations):
        """
        It should be possible to initialize one pipeline and call it in parallel using different
        mutations as input
        """       
        # Read the configuration file and set the variables
        if isinstance(configurations, six.string_types):
            conf.read_configuration_file(configurations)
        elif isinstance(configurations, dict):
            configs.update(**configurations)
        
        # TODO: remove so error message does not appear in a production release.  
        self._validate_temp_path(conf.configs)

        # Initialize a logger
        for line in ELASPIC_LOGO.split('\n'):
            logger.info(line)

        self.PWD = os.getcwd()
        
        self.sequences = {}
        self.models = {}
        self.predictions = {}



    def _validate_temp_path(self, configs):
        """
        Make sure that we are using a job specific temporary folder if we are on a cluster.
        """
        hostname = helper.get_hostname()
        no_job_specific_folder = configs['temp_dir'].startswith(configs['global_temp_dir'])
        on_node_with_manditory_job_specific_folder = (
            any([(x.lower() in hostname) for x in ['node', 'behemoth', 'grendel', 'beagle']])
        )
        if no_job_specific_folder and on_node_with_manditory_job_specific_folder:
            raise Exception('You should be using a temp folder that it specific to the particular job!')
    

    def run(self):
        raise NotImplementedError()
        
        
    def run_sequence(self):
        raise NotImplementedError()


    def _get_model(self):
        raise NotImplementedError()
    
    
    def run_model(self):
        raise NotImplementedError()


    def _get_mutation(self):
        raise NotImplementedError


    def run_mutation(self):
        raise NotImplementedError()

