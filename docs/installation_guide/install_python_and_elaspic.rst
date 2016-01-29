.. _install_python_and_elaspic:

Installing Python and ELASPIC
==============================

#. Download and install the `Anaconda Python Distribution`_ (Python 3) for Linux.

#. Add ``bioconda``, ``salilab``, and ``ostrokach`` channels to your ~/.condarc file::

    conda config --add channels ostrokach
    conda config --add channels salilab
    conda config --add channels bioconda

#. Obtain a `Modeller license`_, and export the license as ``KEY_MODELLER`` in your ~/.bashrc file::

    # ~/.bashrc
    export KEY_MODELLER=XXXXXXX


#. Install ELASPIC and all its dependencies into a new conda environment::

    conda create -n elaspic elaspic

#. Activate the new environment and use elaspic::

    source activate elaspic
    elaspic --help

.. _Conda: http://conda.pydata.org/
.. _Anaconda Python Distribution: https://store.continuum.io/cshop/anaconda/
.. _Miniconda: http://conda.pydata.org/miniconda.html
.. _Modeller license: https://salilab.org/modeller/registration.html
