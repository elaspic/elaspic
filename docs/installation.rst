Installation
============

Install Python using Anaconda
-----------------------------

Downlaod and install the Anaconda Python distribution for linux: http://continuum.io/downloads. 
Run ``pip install -r requirements.txt`` in the elaspic directory to install the python packages required by elaspic.



Install Python using Virtualenv
-------------------------------

It is recommended that you use Python 2.7+ or Python 3.4+ for this project. 
If your system comes with an earlier version of Python, you should download and a more recent version of Python and compile it from source. 
Use ``make altinstall`` instead of ``make install`` to prevent system Python binaries from being overwritten.
You must specify the following options when running *./configure* in order for the required python packages to work::

    --enable-unicode=ucs4 
    --enable-ipv6 
    --with-dbmliborder=gdbm:bdb 
    --with-threads

The ``--enable-unicode=ucs4`` is particularly important, as it is required for third-party packags such as numpy / scipy.

Once you have Python installed in your system, you should set up a virual environment for the elaspic project. 
The instructions on how to do this can be found on the `virtualenv`_ and `virtuanenvwrapper`_ websites. 

Run ``pip install -r requirements.txt`` in the elaspic directory to install the required python packages.



Install external software
-------------------------

Install `BLAST`_, `MODELLER`_, and `TCOFFEE`_, and modify the ``.bashrc`` file accordingly. 
If you install all three programs in your home directory, your ``.bashrc`` file should look something like this::

    # BLAST
    export PATH=$HOME/ncbi-blast-2.2.29+/bin:$PATH

    # TCOFFEE
    export DIR_4_TCOFFEE=$HOME/tcoffee
    export MAFFT_BINARIES=$DIR_4_TCOFFEE/plugins/linux/
    export CACHE_4_TCOFFEE=$HOME/.t_coffee/cache/
    export TMP_4_TCOFFEE=$DIR_4_TCOFFEE/tmp/
    export LOCKDIR_4_TCOFFEE=$DIR_4_TCOFFEE/lck/
    export PERL5LIB=$DIR_4_TCOFFEE/perl:$PERL5LIB
    export EMAIL_4_TCOFFEE=tcoffee.msa@gmail.com # your email goes here
    export PATH=$DIR_4_TCOFFEE/bin:$PATH

    # MODELLER
    export PATH=$HOME/modeller9.14/bin/mod9.14:$PATH
    export LD_LIBRARY_PATH=$HOME/modeller9.14/lib/x86_64-intel8:$LD_LIBRARY_PATH
    # # For Python 2.7+:
    # export PYTHONPATH=$HOME/modeller9.14/modlib:\
    #                   $HOME/modeller9.14/lib/x86_64-intel8/python2.5:\
    #                   $PYTHONPATH
    # For Python 3.4+:
    export PYTHONPATH=$HOME/modeller9.14/modlib:\
                      $HOME/modeller9.14/lib/x86_64-intel8/python3.3:\
                      $PYTHONPATH

All other binaries required to run ELASPIC come precompiled in the ``./elaspic/bin`` folder. 
If you are running a strange flavour of linux, you may have to compile some of those libraries from source. 
The script ``./elaspic/src/install_required_binaries.sh`` is designed to streamline this process.



Download external databases
----------------------------

PDB
~~~

Download the contents of the ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/ folder,
and change the :option:`pdb_path` variable in your configuration file to point to the directory
containing the downloaded data.



Blast
~~~~~

Download and extract the `nr` and `pdbaa` databases from ftp://ftp.ncbi.nlm.nih.gov/blast/db/, 
and change the :option:`blast_db_path` variable in your configuration file to point to the directory
containing the downloaded data.






.. _virtualenv: http://virtualenv.readthedocs.org/en/latest/
.. _virtuanenvwrapper: http://virtualenvwrapper.readthedocs.org/en/latest/
.. _BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
.. _MODELLER: https://salilab.org/modeller/
.. _TCOFFEE: http://www.tcoffee.org/

