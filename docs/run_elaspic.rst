.. _run_elaspic:


Running ELASPIC 
================

If you have followed all the instructions in the :ref:`install_guide` section, you can run 
the ELASPIC pipeline using the ``elaspic`` command::

    elaspic --help

or by running ``python -m elaspic`` from the root of the elaspic repository folder::

    python -m elaspic --help

.. program-output:: elaspic --help


elaspic options
----------------

.. argparse::
   :module: elaspic.__main__
   :func: get_parser
   :prog: elaspic



Working with the ELASPIC database
=================================

In order to make it easier to initialize and upload data to the ELASPIC database, we provide a script
``elaspic_database.py`` which should streamline the most common database tasks, such as creating 
a new schema, uploading data to the database, and testing to make sure that everything works correctly.
If you have followed all the instructions in the :ref:`install_guide` section, you can access
the ``elaspic_database.py`` script using the ``elaspic_database`` command anywhere on your system::

    elaspic_database --help

or by running the script directly from a local clone of the ELASPIC code repository::

    elaspic ./elaspic/elaspic_database.py --help


.. program-output:: elaspic_database --help

elaspic_database options
-------------------------

.. argparse::
   :module: elaspic.elaspic_database
   :func: get_parser
   :prog: elaspic_database






Command-line options
====================

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



-------------------------------------------------------------------------------------------------

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
  

