.. _run_elaspic:


Running ELASPIC 
================

If you have followed all the instructions in the :ref:`install_guide`, you can run 
the ELASPIC pipeline using the ``elaspic`` command::

    elaspic

or by running ``python -m elaspic`` from the root of the elaspic repository folder::

    python -m elaspic


ELASPIC options
----------------

Description of ELASPIC options can be obtained by running ``elaspic --help``:

.. program-output:: python ./../elaspic/__main__.py --help


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



-------------------------------------------------------------------------------------------------

Working with the ELASPIC database
=================================

In order to make it easier to initialize and upload data to the ELASPIC database, we provide a script
``elaspic_database.py`` which should streamline the most common database tasks, such as:

- :ref:`creating <create_new_schema>` a new schema 
- :ref:`uploading <load_data_to_schema>` data to the new schema 
- :ref:`testing <test_database>` the new schema to make sure it works correctly
- :ref:`deleting <delete_database_schema>` the schema if you don't need it or wish to start with a fresh slate

If you have followed all the instructions in the :ref:`install_guide`, you can access
the ``elaspic_database.py`` script using the ``elaspic_database`` command anywhere on your system::

    elaspic_database

or by running the script directly from a local clone of the ELASPIC code repository::

    python ./elaspic/elaspic_database.py


Description of all availible options can be obtained using the ``--help`` command:

.. program-output:: python ./../elaspic/elaspic_database.py --help


.. _create_new_schema:

Create a new database schema
----------------------------

.. program-output:: python ./../elaspic/elaspic_database.py create --help


.. _load_data_to_schema:

Load data to the database
-------------------------

.. program-output:: python ./../elaspic/elaspic_database.py load_data --help


.. _test_database:

Test the database schema
------------------------

.. program-output:: python ./../elaspic/elaspic_database.py test --help


.. _delete_database_schema:

Delete the database schema
--------------------------

.. program-output:: python ./../elaspic/elaspic_database.py delete --help


