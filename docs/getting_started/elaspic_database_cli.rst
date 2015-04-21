.. _elaspic_database_cli:

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

    elaspic_database -c {config_file}

or by running the script directly from a local clone of the ELASPIC code repository::

    python ./elaspic/elaspic_database.py -c {config_file}


Description of all availible options can be obtained using the ``--help`` command:

.. program-output:: elaspic_database --help


.. _create_new_schema:

Create a new database schema
----------------------------

.. program-output:: elaspic_database --help


.. _load_data_to_schema:

Load data to the database
-------------------------

.. program-output:: elaspic_database load_data --help


.. _test_database:

Test the database schema
------------------------

.. program-output:: elaspic_database test --help


.. _delete_database_schema:

Delete the database schema
--------------------------

.. program-output:: elaspic_database delete --help


