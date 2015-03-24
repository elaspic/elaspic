.. _config_file:

Updating the configuration file
===============================

Modify the ELASPIC configuration file ``./config/default_config_file.ini`` to match your system. 
Important options are described below.


Configuration options
---------------------

[DEFAULT]
~~~~~~~~~
.. option:: global_temp_path
   
  Location for storing temporary files. It will be used only if the :envvar:`TMPDIR` environmental variable is not set. Default = "/tmp/".
  
.. option:: temp_path string
   
  A folder in the `global_temp_path` that will contain all the files that are relevant to ELASPIC. A unique folder will be created inside :option:`temp_path` for each job. Default = "elaspic/".
  
.. option:: debug
   
  Whether or not to show detailed debugging information. If True, the logging level will be set to ``logging.DEBUG``. If False, the logging level will be set to ``logging.INFO``. Default = True.
  
.. option:: look_for_interactions
   
  Whether or not to compute models of protein-protein interactions. Default = True.
  
.. option:: remake_provean_supset
   
  Whether or not to remake the Provean supporting set if one or more sequences cannot be found in the BLAST database. Default = False.
  
.. option:: n_cores
   
  Number of cores to use by programs that support multithreading. Default = 1.
  
.. option:: schema_version
   
  Database schema to use for storing and retreiving data. Default = "elaspic".
  
.. option:: web_server
   
  Whether or not the ELASPIC pipeline is being run as part of a webserver. Default = False.

.. option:: provean_temp_path
  
  Location to store provean temporary files if working on any note other than `beagle` or `banting`.
  For internal use only. Default = ''.


[DATABASE]
~~~~~~~~~~
.. option:: db_type

  The database that you are using. Supported databases are `MySQL`, `PostgreSQL`, and `SQLite`.
  
.. option:: sqlite_db_path

  Location of the SQLite database. Required only if :option:`db_type` is `SQLite`.

.. option:: db_schema

  The name of the schema that holds all elaspic data.

.. option:: db_schema_uniprot

  The name of the database schema that holds uniprot sequences. Defaults to :option:`db_schema`.

.. option:: db_database

  The name of the database that contains :option:`db_schema` and :option:`db_schema_uniprot`.
  Required only if :option:`db_type` is `PostgreSQL`. Defaults to :option:`db_schema`. 

.. option:: db_username

  The username for the database. Required only if :option:`db_type` is `MySQL` or `PostgreSQL`. 

.. option:: db_password

  The password for the database. Required only if :option:`db_type` is `MySQL` or `PostgreSQL`. 

.. option:: db_url

  The IP address of the database. Required only if :option:`db_type` is `MySQL` or `PostgreSQL`. 

.. option:: db_port

  The listening port of the database. Required only if :option:`db_type` is `MySQL` or `PostgreSQL`. 


[SETTINGS]
~~~~~~~~~~
.. option:: path_to_archive

  Location for storing and retreiving precalculated data.
  
.. option:: blast_db_path

  Location of the blast **nr** and **pdbaa** databases.

.. option:: pdb_path 

  Location of all pdb structures, equivalent to the "data/data/structures/divided/pdb/" folder in the PDB ftp site. Optional.
  
.. option:: bin_path

  Location of external binary files required by ELASPIC.


[GET_MODEL]
~~~~~~~~~~~
.. option:: modeller_runs

  Number of models that MODELLER should make before choosing the best one. Not implemented! Default = 1.


[GET_MUTATION]
~~~~~~~~~~~~~~
.. option:: foldx_water

  ``-CRYSTAL``: use water molecules in the crystal structure to bridge two protein atoms. 
  
  ``-PREDICT``: predict water molecules that make 2 or more hydrogen bonds to the protein. 
  
  ``-COMPARE``: compare predicted water bridges with bridges observed in the crystal structure.
  
  ``-IGNORE``: don't predict water molecules. Default.
  
  Source: http://foldx.crg.es/manual3.jsp.
  
.. option:: foldx_num_of_runs
  
  Number of times that FoldX should evaluate a given mutation. Default = 1.
  
.. option:: matrix_type
  
  Substitution matrix for calculating the mutation conservation score. Default = "blosum80".
  
.. option:: gap_start 
   
  Penalty for starting a gap when calculating the mutation conservation score. Default = -16.
  
.. option:: gap_extend
   
  Penalty for extending a gap when calculating the mutation conservation score. Default = -4.



Environmental variables
-----------------------

.. envvar:: PATH

  A colon-separated list of paths where ELASPIC should look for required programs, such as BLAST, T-coffee, Modeller, and cd-hit.

.. envvar:: TMPDIR

  Location to store all temporary files and folders.
  

