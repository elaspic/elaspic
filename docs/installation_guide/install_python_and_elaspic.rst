.. _install_python_and_elaspic:

Installing Python and ELASPIC
==============================

.. _install_using_conda:

Using Conda (recommended)
-------------------------

#. Downlaod and install the `Anaconda Python Distribution`_ for Linux.
   
#. Add ``bioconda``, ``salilab``, and ``ostrokach`` channels to your ~/.condarc file::

    conda config --add channels bioconda
    conda config --add channels salilab
    conda config --add channels ostrokach

#. Obtain a `Modeller license`_, and export the license as ``KEY_MODELLER`` in your ~/.bashrc file::

    # ~/.bashrc
    export KEY_MODELLER=XXXXXXX


#. Install ELASPC, including all dependencies::

    conda install elaspic


.. _Conda: http://conda.pydata.org/
.. _Anaconda Python Distribution: https://store.continuum.io/cshop/anaconda/
.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _Modeller license: https://salilab.org/modeller/registration.html


Using Virtualenv
-----------------

.. note:: 

    This method is deprecated. Please use conda instead.
    

#. Make sure that you have Python 2.7+ or Python 3.4+ installed on your system. If your system comes with an older version of Python, download a `recent version of Python`_ and compile it from source. Use ``make altinstall`` instead of ``make install`` to prevent system Python binaries from being overwritten. 

   You must specify the following options when running ``./configure`` in order for the required python packages to work::

    --enable-unicode=ucs4 # ucs4 is required for numpy / scipy
    --enable-ipv6 
    --with-dbmliborder=gdbm:bdb 
    --with-threads

#. Clone the ELASPIC git repository::

    git clone git@bitbucket.org:ostrokach/elaspic.git

#. Create a virual environment containing all ELASPIC dependencies. Detailed instructions on how 
   to do this can be found on the `virtualenv`_ and `virtuanenvwrapper`_ websites. 
   In brief::

    # Install virtualenv via pip
    pip install virtualenv

    # Enter the ELASPIC repository folder
    cd elaspic
    
    # Create a virtual environment called `elaspic`
    # (specify the complete path to your Python 2.7 / Python 3.4 installation)
    virtualenv -p /usr/bin/python2.7 elaspic

    # Activate the `elaspic` virtual environment
    source elaspic/bin/activate

#. Install required Python packages::

    # Install required Python packages
    pip install -r requirements.txt

#. Install the ELASPC package::

    # Install in development mode so that local changes are reflected
    python setup.py develop 


.. _recent version of Python: https://www.python.org/downloads/
.. _virtualenv: http://virtualenv.readthedocs.org/en/latest/
.. _virtuanenvwrapper: http://virtualenvwrapper.readthedocs.org/en/latest/

