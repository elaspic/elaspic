.. _elaspic_cli:
.. _run_elaspic:

Command Line Interface
======================

After following instructions in the :ref:`installation_guide`, you should be able to run ELASPIC
from the command line using the ``elaspic`` command::

  $ elaspic --help
  usage: elaspic [-h] {run,database,train} ...

  optional arguments:
    -h, --help            show this help message and exit

  command:
    {run,database,train}
      run                 Run ELASPIC.
      database            Perform maintenance tasks on the ELASPIC database.
      train               Train the ELASPIC classifiers.

Type ``--help`` to see the options available for each subcommand:

- ``elaspic run --help``
- ``elaspic database --help``
- ``elaspic database load_data --help``
- etc...


elaspic run
-----------

Run the ELASPIC pipeline.

If you wish to mutate an existing PDB, you should specify the name of the PDB file to be mutated, and the mutation(s)::

  elaspic run \
      --structure_file {structure_file} \
      --mutations {mutations}

If you wish to first create a homology model of a protein, you should provide a fasta file containing the sequence of the protein to be modelled, a PDB file containing the structural template, and the mutation(s)::

  elaspic run \
      --sequence_file {sequence_file} \
      --structure_file {structure_file} \
      --mutations {mutations}

If you wish to perform mutagenesis on a proteome-wide scale, you need to download protein domain definitions from the `elaspic downloads page`_, and optionally a local copy of the PDB database. After saving your database information to a configuration file, you can run specify the uniprot id and mutation(s)::

  elaspic run \
    --config_file {config_file} \
    --uniprot_id {uniprot_id} \
    --mutations {mutations}


elaspic train
-------------

Train the machine learning predictor for the ELASPIC pipeline.

This is automatically done at install time, and you *do not* need to do this again unless you update your ``scikit-learn`` version.


.. _`elaspic_database_cli`:

elaspic database
----------------

Perform maintenance tasks on the ELASPIC database.

You must provide a configuration file containing the details of your database installation for any of these commands to work. For more information about configuration files, see :ref:`config_file`.

elaspic database create
~~~~~~~~~~~~~~~~~~~~~~~

Create a new database schema.

elaspic database load_data
~~~~~~~~~~~~~~~~~~~~~~~~~~

Load data to the database.

elaspic database delete
~~~~~~~~~~~~~~~~~~~~~~~

Delete the database schema.


.. _`elaspic downloads page`: http://elaspic.kimlab.org/static/download/
