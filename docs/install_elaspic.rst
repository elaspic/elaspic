Installing ELASPIC
==================

1. Install elaspic in release mode::

    python setup.py install
    
   or in development mode (useful if you are planning to make changes to the code)::
   
    python setup.py develop

2. After you install the required external software (:ref:`install_externals`), 
   download BLAST and PDB databases to a local folder(:ref:`download_data`), 
   and modify the configuration file to match your system and database setting (:ref:`config_file`), 
   you can test the ELASPIC pipeline by entering the ``./tests/`` subdirectory and running::
   
    py.test test_database.py -c {path_to_your_configuration_file}

   **Important!** The test configuration file file should refer to a schema that is reserved for testing,
   as it will be **cleared of all data** at the beginning and in the end of every test.
   
