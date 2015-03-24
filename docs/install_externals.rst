.. _install_externals:

Installing external software
============================

Install `BLAST`_, `Modeller`_, and `Tcoffee`_, and modify the ``.bashrc`` file accordingly. 
If you install all three programs in your home directory, your ``.bashrc`` file should look something like this::

    # BLAST
    export PATH=$HOME/ncbi-blast-2.2.29+/bin:$PATH

    # TCOFFEE
    export DIR_4_Tcoffee=$HOME/tcoffee
    export MAFFT_BINARIES=$DIR_4_Tcoffee/plugins/linux/
    export CACHE_4_Tcoffee=$HOME/.t_coffee/cache/
    export TMP_4_Tcoffee=$DIR_4_Tcoffee/tmp/
    export LOCKDIR_4_Tcoffee=$DIR_4_Tcoffee/lck/
    export PERL5LIB=$DIR_4_Tcoffee/perl:$PERL5LIB
    export EMAIL_4_Tcoffee=tcoffee.msa@gmail.com # your email goes here
    export PATH=$DIR_4_Tcoffee/bin:$PATH

    # MODELLER
    export PATH=$HOME/modeller9.14/bin/mod9.14:$PATH
    export LD_LIBRARY_PATH=$HOME/modeller9.14/lib/x86_64-intel8:$LD_LIBRARY_PATH
    # For Python 2.7+:
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

.. _BLAST: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
.. _Modeller: https://salilab.org/modeller/
.. _Tcoffee: http://www.tcoffee.org/

