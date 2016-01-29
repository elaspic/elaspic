.. _import_precalculated:

Importing precalculated data
=============================

ELASPIC downloads page
----------------------

The `ELASPIC downloads page`_ contains all precalculated data that is required to run the ELASPIC pipeline on a local machine.

The ``*.tsv.gz`` files correspond to different tables of the :ref:`ELASPIC database <database_schema>`:

- The ``domain.tar.gz`` file in the root folder contains Profs domain definitions for files in the PDB, and corresponds to the :ref:`domain` table.
- The ``domain_contact.tar.gz`` file in the root folder contains a list of interactions between those domains, and corresponds to the :ref:`domain_contact` table.
- All other tables are split into separate folders according to the organism of origin. The files are named using the ``{table_name}.tsv.gz`` convention, where ``table_name`` is the name of the table in the database.

The ``*.7z`` files contain precalculated data:

- The *provean*, *uniprot_domain*, and *uniprot_domain_pair* subfolders contain precalculated provean supporting sets, and homology models of protein domains and domain-domain interactions, respectively.

Precalculated mutations:

- The *Homo_sapiens* folder contains an additional subfolder `precalculated_mutations`_, which contains :math:`\Delta \Delta G` scores for mutations in various datasets.

.. note::

  The `configure_test.sh`_ and `run_test.sh`_ scripts in the *./scripts* folder contain examples of how to download and set up a local copy of the database.


Downloading data
----------------

In order to run up ELASPIC on a local computer, you need to download precalculated data
for your organism of interest. If your goal is to only test the pipeline, you can download a test dataset from the folder `current_release/Homo_sapiens_test`_.

To download all precalculated data for a given organism, use the ``wget`` command:

.. code-block:: bash

    # Download external files
    wget -P "${TEST_DIR}/elaspic.kimlab.org" \
        http://elaspic.kimlab.org/static/download/current_release/domain.tsv.gz
    wget -P "${TEST_DIR}/elaspic.kimlab.org" \
        http://elaspic.kimlab.org/static/download/current_release/domain_contact.tsv.gz
    wget -P "${TEST_DIR}" \
        -r --no-parent --reject "index.html*" --cut-dirs=4  \
        http://elaspic.kimlab.org/static/download/current_release/Homo_sapiens_test/

You need to extract the provean supporting sets and domain homology models into a folder
specified by the :term:`archive_dir` variable in your :ref:`configuration_file <config_file>`:

.. code-block:: bash

    mkdir archive  # Set 'archive_dir' variable in the config file to this folder

    7z x "${TEST_DIR}/elaspic.kimlab.org/provean/provean.7z" -o"archive"
    7z x "${TEST_DIR}/elaspic.kimlab.org/uniprot_domain/uniprot_domain.7z" -o"archive"
    7z x "${TEST_DIR}/elaspic.kimlab.org/uniprot_domain_pair/uniprot_domain_pair.7z" -o"archive"



Importing data into a database
------------------------------

You also need to create a local SQL database and fill it with precalculated data.

Modify the database variables in the ELASPIC :ref:`configuration file <config_file>` to
match your local *MySQL*, *PostgreSQL*, or *SQLite* database, and use the :ref:`elaspic database <elaspic_database_cli>` CLI to create a new database and fill it with precalculated data.

First, you need to create an empty database::

    elaspic database -c {your_configuration_file}.ini create

Next, you need to load all precalculated data for the organism in question to your database::

    elaspic database -c {your_configuration_file}.ini load_data

To delete the database that you just created, run::

    elaspic database -c {your_configuration_file}.ini delete


.. _ELASPIC downloads page: http://elaspic.kimlab.org/static/download/current_release/

.. _`configure_test.sh`: https://github.com/ostrokach/elaspic/blob/master/scripts/configure_test.sh
.. _`run_test.sh`: https://github.com/ostrokach/elaspic/blob/master/scripts/run_test.sh
.. _`scripts/`: https://github.com/ostrokach/elaspic/blob/master/scripts/

.. _`current_release/Homo_sapiens_test`: http://elaspic.kimlab.org/static/download/current_release/Homo_sapiens_test/
.. _`precalculated_mutations`: http://elaspic.kimlab.org/static/download/current_release/Homo_sapiens/precalculated_mutations/
