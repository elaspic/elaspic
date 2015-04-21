.. _run_elaspic:
.. _elaspic_cli:

Running ELASPIC 
================

If you have followed all the instructions in the :ref:`install_guide`, you can run 
the ELASPIC pipeline using the ``elaspic`` command::

    elaspic -c {config_file} -u {uniprot_id} -m {mutation(s)}

or by running ``python -m elaspic`` from the root of the elaspic repository folder::

    python -m elaspic -c {config_file} -u {uniprot_id} -m {mutation(s)}


ELASPIC options
----------------

Description of ELASPIC options can be obtained by running ``elaspic --help``:

.. program-output:: elaspic --help


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

