Installing Python
=================

Using Anaconda
--------------

1. Downlaod and install the `Anaconda Python distribution <http://continuum.io/downloads>`_ for linux. 

2. Create a new environment containing the required Anaconda binaries::

    conda -n elaspic --file requirements_conda_{py2/py3}.txt

3. Activate the `elaspic` environment::

    source activate elaspic

4. Install the remaining packages that are only availible through pip::

    pip install -r requirements_conda_pip.txt

   Note: If you don't want to install ELASPIC in a separate Anaconda environment, 
   remove ``-n elaspic`` from step 2, and skip step 3.


Using Virtualenv
-----------------

It is recommended that you use Python 2.7+ or Python 3.4+ for this project. If your system comes 
with an earlier version of Python, you should download and a more recent version of Python and 
compile it from source. Use ``make altinstall`` instead of ``make install`` to prevent system 
Python binaries from being overwritten.You must specify the following options when running 
``./configure`` in order for the required python packages to work::

    --enable-unicode=ucs4 
    --enable-ipv6 
    --with-dbmliborder=gdbm:bdb 
    --with-threads

The ``--enable-unicode=ucs4`` is particularly important, as it is required for third-party packags 
such as numpy / scipy.

Once you have Python installed in your system, you should set up a virual environment for the elaspic project. 
The instructions on how to do this can be found on the `virtualenv`_ and `virtuanenvwrapper`_ websites. 

Run ``pip install -r requirements_virtualenv_{py2/py3}.txt`` in the elaspic directory to install 
the required python packages.


.. _virtualenv: http://virtualenv.readthedocs.org/en/latest/
.. _virtuanenvwrapper: http://virtualenvwrapper.readthedocs.org/en/latest/

