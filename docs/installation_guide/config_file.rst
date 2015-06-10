.. _config_file:

Updating the configuration file
================================

Edit the ELASPIC configuration file `./config/config_file.ini`_ to match your system:
  #. Settings in the :ref:`[SETTINGS]` section should be modified to match the local installation of ELASPIC and other programs and datasets.

  #. Settings in the :ref:`[DATABASE]` section should be modified to match the local MySQL, PostgreSQL, or SQLite database.

  #. Settings in the :ref:`[DEFAULT]`, :ref:`[GET_MODEL]`, and :ref:`[GET_MUTATION]` sections often may be left unchanged, since the default values are appropriate for most systems.


.. _`./config/config_file.ini`: https://bitbucket.org/ostrokach/elaspic/src/b7f17ea1bcc704dc3bdcaecd7286b5fedb300151/config/config_file.ini



---------------------------------------------------------------------------------------------------

Configuration options
----------------------

.. _`[DEFAULT]`:

[DEFAULT]
~~~~~~~~~

.. glossary::

   global_temp_path
     Location for storing temporary files. It will be used only if the :envvar:`TMPDIR` environmental variable is not set. **Default = '/tmp/'**.
     
   temp_path string
     A folder in the :term:`global_temp_path` that will contain all the files that are relevant to ELASPIC. Inside this folder, every job will create its own unique subfolder. **Default = 'elaspic/'**.
     
   debug
     Whether or not to show detailed debugging information. If True, the logging level will be set to ``logging.DEBUG``. If False, the logging level will be set to ``logging.INFO``. **Default = True**.
     
   look_for_interactions
     Whether or not to compute models of protein-protein interactions. **Default = True**.
     
   remake_provean_supset
     Whether or not to remake the Provean supporting set if one or more sequences cannot be found in the BLAST database. **Default = False**.
     
   n_cores
     Number of cores to use by programs that support multithreading. **Default = 1**.
     
   schema_version
     Database schema to use for storing and retreiving data. **Default = 'elaspic'**.
     
   web_server
     Whether or not the ELASPIC pipeline is being run as part of a webserver. **Default = False**.

   provean_temp_path
     Location to store provean temporary files if working on any note other than `beagle` or `banting`.
     For internal use only. **Default = ''**.

  copy_data
    Whether or not to copy calculated data back to the archive. Set to 'False' if you are planning 
    to copy the data yourself (e.g. from inside a PBS or SGE script). **Default = True**.


.. _`[DATABASE]`:

[DATABASE]
~~~~~~~~~~

.. glossary::

  db_type
    The database that you are using. Supported databases are `MySQL`, `PostgreSQL`, and `SQLite`.

  sqlite_db_path
    Location of the SQLite database. Required only if :term:`db_type` is `SQLite`.

  db_schema
    The name of the schema that holds all elaspic data.

  db_schema_uniprot
    The name of the database schema that holds uniprot sequences. Defaults to :term:`db_schema`.

  db_database
    The name of the database that contains :term:`db_schema` and :term:`db_schema_uniprot`.
    Required only if :term:`db_type` is `PostgreSQL`. Defaults to :term:`db_schema`. 

  db_username
    The username for the database. Required only if :term:`db_type` is `MySQL` or `PostgreSQL`. 

  db_password
    The password for the database. Required only if :term:`db_type` is `MySQL` or `PostgreSQL`. 

  db_url
    The IP address of the database. Required only if :term:`db_type` is `MySQL` or `PostgreSQL`. 

  db_port
    The listening port of the database. Required only if :term:`db_type` is `MySQL` or `PostgreSQL`. 

  db_socket
    Path to the socket file, if it is not in the default location.
    Used only if :term:`db_url` is `localhost`. 
    For example: ``/usr/local/mysql5/mysqld.sock`` for `MySQL` and ``/var/lib/postgresql`` for `PostgreSQL`.



.. _`[SETTINGS]`:

[SETTINGS]
~~~~~~~~~~

.. glossary::

   path_to_archive
     Location for storing and retrieving precalculated data.
    
   blast_db_path
     Location of the blast **nr** and **pdbaa** databases.

   pdb_path 
     Location of all pdb structures, equivalent to the "data/data/structures/divided/pdb/" folder in the PDB ftp site. Optional.
     


.. _`[GET_MODEL]`:

[GET_MODEL]
~~~~~~~~~~~

.. glossary::

   modeller_runs
     Number of models that MODELLER should make before choosing the best one. Not implemented! **Default = 1**.



.. _`[GET_MUTATION]`:

[GET_MUTATION]
~~~~~~~~~~~~~~

.. glossary::

   foldx_water
     - ``-CRYSTAL``: use water molecules in the crystal structure to bridge two protein atoms. 
     - ``-PREDICT``: predict water molecules that make 2 or more hydrogen bonds to the protein. 
     - ``-COMPARE``: compare predicted water bridges with bridges observed in the crystal structure.
     - ``-IGNORE``: don't predict water molecules. **Default**.

     Source: http://foldx.crg.es/manual3.jsp.
     
   foldx_num_of_runs
     Number of times that FoldX should evaluate a given mutation. **Default = 1**.
     
   matrix_type
     Substitution matrix for calculating the mutation conservation score. **Default = 'blosum80'**.
     
   gap_start 
     Penalty for starting a gap when calculating the mutation conservation score. **Default = -16**.
     
   gap_extend
     Penalty for extending a gap when calculating the mutation conservation score. **Default = -4**.



Environmental variables
------------------------

.. envvar:: PATH

  A colon-separated list of paths where ELASPIC should look for required programs, such as BLAST, T-coffee, Modeller, and cd-hit.

.. envvar:: TMPDIR

  Location to store all temporary files and folders.
  

