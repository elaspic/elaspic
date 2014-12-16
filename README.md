## INSTALLATION INSTRUCTIONS

### Configuring python

#### Using Anaconda

Downlaod and install the anaconda python 2.7 distribution for linux: 
<http://continuum.io/downloads>.

After adding the anaconda installation folder to the system path, elaspic 
should run successfully.


#### Using virtualenv

It is recommended that you use python 2.7 for this project. If your system
comes with an earlier version of python, you should download and compile 
python 2.7 binaries from source. Use ``make altinstall`` instead of 
``make install`` to prevent any system python binaries from being overwritten.
Also, you must specify the following options while running *./configure* 
in order for the required python packages to work: 

    --enable-unicode=ucs4 
    --enable-ipv6 
    --with-dbmliborder=gdbm:bdb 
    --with-threads

The ``--enable-unicode=ucs4`` is particularly important, as it is required
for 3rd-party packags such as numpy / scipy.


After you have installed python 2.7 binaries, it is recommended that you 
set up a virual python environment to work with the elaspic project. 
The instructions on how to do this can be found on the 
[virtualenv](http://virtualenv.readthedocs.org/en/latest/) and
[virtuanenvwrapper](http://virtualenvwrapper.readthedocs.org/en/latest/)
websites.

Once you activate the virtual environment, the python packages required
for elaspic can be installed using ``pip install -r requirements.txt``.
If you would like a development environment with a working copy of 
[Spyder](https://bitbucket.org/spyder-ide/spyderlib/downloads),
you also have to run ``pip install -r requirements_dev.txt`` and copy the
*elaspic/src/postmkvirtualenv* file to your *.virtualenvs* path (or 
manually create the symbolic links specified in that file). 


### Installing external software

Install [BLAST](ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/), 
[MODELLER](https://salilab.org/modeller/), and [TCOFFEE](http://www.tcoffee.org/),
and modify the .bashrc file accordingly.
If you install all three programs in your home directory, your bashrc
file should look something like this:

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
    export PYTHONPATH=$HOME/modeller9.13/modlib:$HOME/modeller9.13/lib/x86_64-intel8/python2.5:$PYTHONPATH

All other binaries required to run elaspic come precompiled in the *elaspic/bin*
folder. If you are running a strange flavour of linux, you may have to compile
some of those libraries from source. The script *elaspic/src/install_required_binaries.sh*
is designed to streamline this process.

