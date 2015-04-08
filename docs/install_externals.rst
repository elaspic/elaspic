.. _install_externals:

Installing external software
=============================

Blast, Modeller, Tcoffee
------------------------

Install `Blast`_, `Modeller`_, and `Tcoffee`_, by following the installation instructions provided 
by the developers, and modify your ``.bashrc`` file accordingly.
If you install all three programs in your HOME directory, your ``.bashrc`` file should look something like this::

    # BLAST
    export PATH=$HOME/ncbi-blast-2.2.29+/bin:$PATH

    # TCOFFEE
    export DIR_4_Tcoffee=$HOME/tcoffee
    export PATH=$DIR_4_Tcoffee/bin:$PATH
    export MAFFT_BINARIES=$DIR_4_Tcoffee/plugins/linux/
    export PERL5LIB=$DIR_4_Tcoffee/perl:$PERL5LIB
    export CACHE_4_Tcoffee=$HOME/.t_coffee/cache/
    export TMP_4_Tcoffee=$DIR_4_Tcoffee/tmp/
    export LOCKDIR_4_Tcoffee=$DIR_4_Tcoffee/lck/
    export EMAIL_4_Tcoffee=tcoffee.msa@gmail.com # your email goes here

    # MODELLER
    export PATH=$HOME/modeller9.14/bin:$PATH
    export LD_LIBRARY_PATH=$HOME/modeller9.14/lib/x86_64-intel8:$LD_LIBRARY_PATH
    # For Python 2.7+:
    # export PYTHONPATH=$HOME/modeller9.14/modlib:\
    #                   $HOME/modeller9.14/lib/x86_64-intel8/python2.5:\
    #                   $PYTHONPATH
    # For Python 3.4+:
    export PYTHONPATH=$HOME/modeller9.14/modlib:\
                      $HOME/modeller9.14/lib/x86_64-intel8/python3.3:\
                      $PYTHONPATH


FoldX
-----

Download FoldX 3 from the `FoldX website`_, and place the ``foldx3b6`` and ``rotabase.txt`` 
files in the ``elaspic/bin`` folder, renaming ``foldx3b6`` to ``foldx64Linux``.
If you are using Ubuntu linux, the ``foldx64Linux`` binary included in ``elaspic/bin`` may already
work on your system, in which case you can skip this step.


Other programs
--------------

All other programs required to run ELASPIC come precompiled in the ``elaspic/bin`` folder. 
If you are running a strange flavour of Linux (or Mac), you may have to compile some of those programs from source. 
Running the script :ref:`install_required_binaries_script`
from the root folder of ELASPIC should take care
of this for you, by downloading and compiling the binaries, and moving them to the ``elaspic/bin`` folder.
   

.. _Blast: ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
.. _Modeller: https://salilab.org/modeller/
.. _Tcoffee: http://www.tcoffee.org/
.. _FoldX website: http://foldx.crg.es/
