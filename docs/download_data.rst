.. _download_data:

Downloading external databases
==============================

PDB
---

Download the contents of the ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/ folder,
and change the :option:`pdb_path` variable in your configuration file to point to the directory
containing the downloaded data.


Blast
-----

Download and extract the `nr` and `pdbaa` databases from ftp://ftp.ncbi.nlm.nih.gov/blast/db/, 
and change the :option:`blast_db_path` variable in your configuration file to point to the directory
containing the downloaded data.


Profs
-----

Download protein domain definition and structural template data from the `ELASPIC downloads page`_.
Upload data to a local database using the `xxx.py` script.


.. _ELASPIC downloads page: http://elaspic.kimlab.org/static/download/elaspic/
