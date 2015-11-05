# -*- coding: utf-8 -*-
import os
from setuptools import setup
#~ import distutils.command.bdist_conda


def read(fname):
    """Utility function to read the README file.

    Used for the long_description.  It's nice, because now 
      1) we have a top level README file and 
      2) it's easier to type in the README file than to put a raw string in below ...

    Source: https://pythonhosted.org/an_example_pypi_project/setuptools.html
    """
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
    

setup(
    name='elaspic',
    version='1.0.0', # now in meta.yaml
    
    description=(
		'Ensemble Learning Approach for Stability Prediction of Interface '
		'and Core mutations (ELASPIC).'),
    author='kimlab',
    author_email='alex.strokach@utoronto.ca',
    url='http://elaspic.kimlab.org',
    packages=['elaspic'],
    package_data={'elaspic': ['data/*']},
    long_description=read("README.rst"),
    
    # Conda specific features
    #~ distclass=distutils.command.bdist_conda.CondaDistribution,
    #~ conda_buildnum=1,
    #~ conda_features=['mkl'],
    #~ conda_import_tests=False,
    
    # Specify install requirements in the conda `meta.yaml` file
    # install_requires=[],
    tests_require=[
        'pytest',
    ],
    entry_points={
          'console_scripts': [
              'elaspic = elaspic.__main__:elaspic',
              'elaspic_database = elaspic.__main__:elaspic_database',
          ]
      },
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Topic :: Structural Biology",
        "Topic :: Bioinformatics",
    ],
)

