.. _run_tests:

Running tests
==============

You can test the ELASPIC pipeline by entering the ``./tests/`` subdirectory and running::

    py.test test_database.py -c {path_to_your_configuration_file}

**Important!** The test configuration file file should refer to a schema that is reserved for testing,
as it will be **cleared of all data** at the beginning and in the end of every test.
   
