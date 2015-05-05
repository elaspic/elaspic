.. _import_precalculated:

Importing precalculated data
=============================


ELASPIC downloads page
----------------------

The `ELASPIC downloads page`_ contains all precalculated data that is required to run the ELASPIC pipeline on a local machine.

- The ``domain.tar.gz`` file in the root folder contains Profs domain definitions for files in the PDB, and corresponds to the :ref:`domain` table in the :ref:`ELASPIC database <database_schema>`.
- The ``domain_contact.tar.gz`` file in the root folder contains a list of interactions between those domains, and corresponds to the :ref:`domain_contact` tables in :ref:`ELASPIC database <database_schema>`.
- All other tables are split into separate folders according to the organism of origin. The files are named using the ``{table_name}.tsv.gz`` convention, where ``table_name`` is the name of the table in the :ref:`ELASPIC database <database_schema>` from which the file originates. 
- The ``provean``, ``uniprot_domain``, and ``uniprot_domain_pair`` subfolders contain precalculated provean supporting sets, and homology models of protein domains and domain-domain interactions, respectively.



Downloading data
----------------

In order to run up ELASPIC on a local computer, you need to download precalculated data 
for your organism of interest. If your goal is to only test the pipeline, you can download a test dataset from the folder `current_release/Homo_sapiens_test`_.

To download all precalculated data for a given organism, use the ``wget`` command::

    wget -r -l2 --no-parent --no-directories -A .tar.gz,.tsv.gz \
        http://elaspic.kimlab.org/static/download/current_release/Homo_sapiens_test 


You need to extract the provean supporting sets and domain homology models into a folder
specified by the :term:`path_to_archive` variable in your :ref:`configuration_file <config_file>`::
 
    mkdir archive # folder specified by the 'path_to_archive' variable in the config file

    tar xzf provean.tar.gz -C archive/
    tar xzf uniprot_domain.tar.gz -C ./archive/
    tar xzf uniprot_domain_pair.tar.gz -C ./archive/



Importing data into a database
------------------------------

You also need to create a local SQL database and fill it with precalculated data.
Modify the database variables in the ELASPIC configuration file to 
match your local ``MySQL``, ``PostgreSQL``, or ``SQLite`` database, and use the ``elaspic_database`` script to create a new database a fill it with precalculated data, as described in the :ref:`elaspic_database_cli` Section. 


First, you need to create an empty database::

    elaspic_database -c {your_configuration_file}.ini create


Next, you need to load all precalculated data for the organism in question to your database::

    elaspic_database -c {your_configuration_file}.ini load_data


To make sure everything ran successfully, you can run an example core and an example interface mutation::

    elaspic_database -c {your_configuration_file}.ini test


.. centered:: **The ELASPIC pipeline is now set up and ready to be used!**


.. _ELASPIC downloads page: http://elaspic.kimlab.org/static/download/
.. _`current_release/Homo_sapiens_test`: http://elaspic.kimlab.org/static/download/current_release/Homo_sapiens_test/
