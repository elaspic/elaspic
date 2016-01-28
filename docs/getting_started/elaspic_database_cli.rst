.. _elaspic_database_cli:

Working with the ELASPIC database
=================================

In order to make it easier to initialize and upload data to the ELASPIC database, we provide a script
``elaspic_database.py`` which should streamline the most common database tasks, such as:

- :ref:`creating <create_new_schema>` a new schema
- :ref:`uploading <load_data_to_schema>` data to the new schema
- :ref:`deleting <delete_database_schema>` the schema if you don't need it or wish to start with a fresh slate

If you have followed all the instructions in the :ref:`install_guide`, you can access
the ``elaspic_database.py`` script using the ``elaspic_database`` command anywhere on your system::

    elaspic_database -c {config_file}

or by running the script directly from a local clone of the ELASPIC code repository::

    python ./elaspic/elaspic_database.py -c {config_file}


Description of all availible options can be obtained using the ``--help`` command:

  elaspic database --help


.. _create_new_schema:

Create a new database schema
----------------------------

  elaspic database create --help


.. _load_data_to_schema:

Load data to the database
-------------------------

  elaspic database load_data --help


.. _delete_database_schema:

Delete the database schema
--------------------------

  elaspic database delete --help
