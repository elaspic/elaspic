Installation
============

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

