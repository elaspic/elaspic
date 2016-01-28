.. _run_elaspic:
.. _elaspic_cli:

Running ELASPIC
================

After following instructions in the :ref:`install_guide`, you should be able to run ELASPIC
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


Options
-------

``elaspic run``
~~~~~~~~~~~~~~~

Run the ELASPIC pipeline.

If you wish to mutate an existing PDB, you should specify the name of the PDB file to be mutated, and the mutation::

    elaspic run \
        --structure_file {structure_file} \
        --mutations {mutations}

If you wish to mutate a homology model of the structure, you should provide a file containing the sequence of the protein that you wish to model, the PDB template of the structure to be used for homology modelling, and the mutation::

  elaspic run \
      --sequence_file {sequence_file} \
      --structure_file {structure_file} \
      --mutations {mutations}

If you wish to perform mutagenesis on a proteome-wide scale, you should first load domain definitions from the

**Use ``elaspic run --help`` to see all available commands.**


Run the ELASPIC pipeline.

To see all the available options, type ``run elaspic --help``.
You need to specify a configuration file with the mutation,


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
