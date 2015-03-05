Configuration
=============

Configuration file
------------------

Modify the ELASPIC configuration file ``./config/config_file.ini`` to match your system. Important options are described below.


Configuration options
~~~~~~~~~~~~~~~~~~~~~

[DEFAULT]
`````````
.. option:: global_temp_path
   
  Location for storing temporary files. It will be used only if the :envvar:`TMPDIR` environmental variable is not set. Default = "/tmp/".
  
.. option:: temp_path
   
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


[SETTINGS]
``````````
.. option:: path_to_archive

  Location for storing and retreiving precalculated data.
  
.. option:: blast_db_path

  Location of the blast **nr** and **pdbaa** databases.
  
.. option:: sqlite_db_path

  Location of the SQLite database with the precalculated data. Optional.
  
.. option:: pdb_path 

  Location of all pdb structures, equivalent to the "data/data/structures/divided/pdb/" folder in the PDB ftp site. Optional.
  
.. option:: bin_path

  Location of external binary files required by ELASPIC.


[GET_MODEL]
```````````
.. option:: modeller_runs

  Number of models that MODELLER should make before choosing the best one. Not implemented! Default = 1.


[GET_MUTATION]
``````````````
.. option:: foldx_water

  If "-CRYSTAL" uses the X-ray waters bridging two protein atoms. If "-PREDICT", waters that make 2 or more hydrogen bonds to the protein are predicted. If "-COMPARE" it compares the predicted water bridges with the X-ray ones. (Source: http://foldx.crg.es/manual3.jsp). Default = "-IGNORE"
  
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

.. envvar:: TMPDIR

  Location to store all temporary files and folders.
  

