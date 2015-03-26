.. _install_python:

Installing Python and ELASPIC
=============================

Using Anaconda (recommended)
----------------------------

#. Downlaod and install the `Anaconda`_ Python distribution for Linux.
   
#. Clone the ELASPIC git repository::

    git clone git@bitbucket.org:ostrokach/elaspic.git

#. Create a separate conda environment for ELASPIC, containing all required Python packages::

    # Enter the ELASPIC repository folder
    cd elaspic

    # Create a conda environment `elaspic`
    conda env create -f environment_{py2/py3}.yml # Use 'py2' for Python2, and 'py3' for Python3

    # Activate the conda environment `elaspic`
    source activate elaspic

#. Install the ELASPC package::

    python setup.py install


.. _Conda: http://conda.pydata.org/
.. _Anaconda: https://store.continuum.io/cshop/anaconda/
.. _Miniconda: http://conda.pydata.org/miniconda.html


Using Virtualenv
-----------------

#. Make sure that you have Python 2.7+ or Python 3.4+ installed on your system. If your system comes with an older version of Python, download a `recent version of Python`_ and compile it from source. Use ``make altinstall`` instead of ``make install`` to prevent system Python binaries from being overwritten. 

   You must specify the following options when running ``./configure`` in order for the required python packages to work::

    --enable-unicode=ucs4 # ucs4 is required for numpy / scipy
    --enable-ipv6 
    --with-dbmliborder=gdbm:bdb 
    --with-threads

.. _recent version of Python: https://www.python.org/downloads/

#. Create a virual environment containing all ELASPIC dependencies. The instructions on how to do this can be found on the `virtualenv`_ and `virtuanenvwrapper`_ websites. 


#. Install required Python packages::

    # Clone the ELASPIC code repository
    git clone git@bitbucket.org:ostrokach/elaspic.git

    # Enter the repository folder
    cd elaspic

    # Install required Python packages
    pip install -r requirements.txt


#. Install the ELASPC package::

    python setup.py install


.. _virtualenv: http://virtualenv.readthedocs.org/en/latest/
.. _virtuanenvwrapper: http://virtualenvwrapper.readthedocs.org/en/latest/

