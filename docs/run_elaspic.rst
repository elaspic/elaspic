.. _run_elaspic:

Running ELASPIC
===============

You can run the ELASPIC pipeline using the ``elaspic.py`` file in the root folder of the repository::

    $ ./elaspic.py -c {config_file} -u {uniprot_id} -m {mutation} -t {run_type}


To get a list of options accepted by ``elaspic.py``, use the ``--help`` command::

   $ ./elaspic.py --help

   usage: elaspic.py [-h] -c CONFIG_FILE [-i INPUT_FILE] [-u UNIPROT_ID]
                     [-m [MUTATIONS]] [-t [RUN_TYPE]] [-n [N_CORES]]

   optional arguments:
     -h, --help            show this help message and exit
     -c CONFIG_FILE, --config_file CONFIG_FILE
     -i INPUT_FILE, --input_file INPUT_FILE
     -u UNIPROT_ID, --uniprot_id UNIPROT_ID
     -m [MUTATIONS], --mutations [MUTATIONS]
     -t [RUN_TYPE], --run_type [RUN_TYPE]
     -n [N_CORES], --n_cores [N_CORES]


The options are described below:

.. option:: -c --config_file
   
  Full path to the ELASPIC configuration file (See :ref:`config_file`).
  
.. option:: -u --uniprot_id
   
  The Uniprot ID of the protein that you wish to analyse.
  
.. option:: -m --mutation
   
  The mutation that you are interested in, in Uniprot coordinates (e.g. 'A20L').
  
.. option:: -t --run_type
   
   1. Run Provean.
   2. Create model.
   3. Evaluate mutation.
   4. Create model and evaluate mutation.
   5. Run Provean, create model, and evaluate mutation.

   Default = 5.
  
.. option:: -n --n_cores
   
  Number of cores to use by programs that support multithreading. Not tested. Default = 1. 
  

