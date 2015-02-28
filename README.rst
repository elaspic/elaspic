Introduction
============

Welcome to the ELASPIC code repository! 

Complete documentation is availible on `ReadTheDocs <http://elaspic.readthedocs.org>`_.
For a small number of mutations, you can try running ELASPIC using our `webserver <http://elaspic.kimlab.org/>`_.

.. 
   Continuous testing runs on drone.io:
   .. image:: https://drone.io/bitbucket.org/ostrokach/elaspic/status.png



Installation instructions
=========================

Installing Python using Anaconda
--------------------------------

Downlaod and install the Anaconda python 2.7 distribution for linux: 
http://continuum.io/downloads. Run ``pip install -r requirements.txt``
in the elaspic directory to install the python packages required by elaspic.


Using Python using virtualenv
-----------------------------

It is recommended that you use python 2.7 for this project. If your system
comes with an earlier version of python, you should download and compile 
python2.7 binaries from source. Use ``make altinstall`` instead of 
``make install`` to prevent any system python binaries from being overwritten.
Also, you must specify the following options when running *./configure* 
in order for the required python packages to work::

    --enable-unicode=ucs4 
    --enable-ipv6 
    --with-dbmliborder=gdbm:bdb 
    --with-threads


The ``--enable-unicode=ucs4`` is particularly important, as it is required
for third-party packags such as numpy / scipy.

Once you have python 2.7 installed in your system, you should set up a virual environment 
for the elaspic project. The instructions on how to do this can be found on the 
`virtualenv`_ and `virtuanenvwrapper`_ websites. 

Run ``pip install -r requirements.txt`` in the elaspic directory to install the required 
python packages.


Installing external software
----------------------------

Install `BLAST`_, `MODELLER`_, and `TCOFFEE`_, and modify the .bashrc file
accordingly. If you install all three programs in your home directory, 
your ``.bashrc`` file should look something like this::

    # BLAST
    export PATH=$HOME/ncbi-blast-2.2.29+/bin:$PATH

    # TCOFFEE
    export DIR_4_TCOFFEE=$HOME/tcoffee
    export MAFFT_BINARIES=$DIR_4_TCOFFEE/plugins/linux/
    export CACHE_4_TCOFFEE=$HOME/.t_coffee/cache/
    export TMP_4_TCOFFEE=$DIR_4_TCOFFEE/tmp/
    export LOCKDIR_4_TCOFFEE=$DIR_4_TCOFFEE/lck/
    export PERL5LIB=$DIR_4_TCOFFEE/perl:$PERL5LIB
    export EMAIL_4_TCOFFEE=alex.strokach@utoronto.ca
    export PATH=$DIR_4_TCOFFEE/bin:$PATH

    # MODELLER
    export PATH=$HOME/modeller9.13/bin/mod9.13:$PATH
    export LD_LIBRARY_PATH=$HOME/modeller9.13/lib/x86_64-intel8:$LD_LIBRARY_PATH
    export PYTHONPATH=$HOME/modeller9.13/modlib:\
                      $HOME/modeller9.13/lib/x86_64-intel8/python2.5:\
                      $PYTHONPATH


All other binaries required to run elaspic come precompiled in the ``./elaspic/bin``
folder. If you are running a strange flavour of linux, you may have to compile
some of those libraries from source. The script ``./elaspic/src/install_required_binaries.sh``
is designed to streamline this process.

.. _virtualenv: http://virtualenv.readthedocs.org/en/latest/
.. _virtuanenvwrapper: http://virtualenvwrapper.readthedocs.org/en/latest/
.. _BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
.. _MODELLER: https://salilab.org/modeller/
.. _TCOFFEE: http://www.tcoffee.org/



Configuration options
=====================

Configuration file
------------------

Modify the ELASPIC configuration file ``./config/config_file.ini`` to match your system. Important options are described below.

[DEFAULT]
~~~~~~~~~

global_temp_path : str, default='/tmp/'
  Location for storing temporary files. It will be used only if the :envvar:`TMPDIR` environmental variable is not set. 
temp_path : str, default='elaspic/'
  A folder in the `global_temp_path` that will contain all the files that are relevant to ELASPIC. A unique folder will be created inside `temp_path` for each job.
debug : bool, default=True
  Whether or not to show detailed debugging information. If True, the logging level will be set to ``logging.DEBUG``. If False, the logging level will be set to ``logging.INFO``.
look_for_interactions : bool, default=True
  Whether or not to compute models of protein-protein interactions.
remake_provean_supset : bool, default=False
  Whether or not to remake the Provean supporting set if one or more sequences cannot be found in the BLAST database.
n_cores : int, default=1
  Number of cores to use by programs that support multithreading.
schema_version : str, default='elaspic'
  Database schema to use for storing and retreiving data.
web_server : bool, default=False
  Whether or not the ELASPIC pipeline is being run as part of a webserver.

[SETTINGS]
~~~~~~~~~~

path_to_archive : str
  Location for storing and retreiving precalculated data.
blast_db_path : str
  Location of the blast **nr** and **pdbaa** databases.
sqlite_db_path : str, optional
  Location of the SQLite database with the precalculated data.
pdb_path : str, optional
  Location of all pdb structures, equivalent to the "data/data/structures/divided/pdb/" folder in the PDB ftp site.
bin_path
  Location of external binary files required by ELASPIC.

[GET_MODEL]
~~~~~~~~~~~

modeller_runs : int, default=1
  Number of models that MODELLER should make before choosing the best one. Not implemented!

[GET_MUTATION]
~~~~~~~~~~~~~~

foldx_water : str, default='-IGNORE'
  If '-CRYSTAL' uses the X-ray waters bridging two protein atoms. If '-PREDICT', waters that make 2 or more hydrogen bonds to the protein are predicted. If '-COMPARE' it compares the predicted water bridges with the X-ray ones. (Source: http://foldx.crg.es/manual3.jsp).  
foldx_num_of_runs : int, default=1
  Number of times that FoldX should evaluate a given mutation.
matrix_type : str, default='blosum80'
  Substitution matrix for calculating the mutation conservation score.
gap_start : int, default=-16
  Penalty for starting a gap when calculating the mutation conservation score.
gap_extend : int, default=-4
  Penalty for extending a gap when calculating the mutation conservation score.



Environmental variables
-----------------------

.. envvar:: TMPDIR

  Location to store all temporary files and folders.
  


