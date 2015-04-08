import os
from setuptools import setup


# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
# Source: https://pythonhosted.org/an_example_pypi_project/setuptools.html
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()
    

setup(
    name='elaspic',
    version='0.1.2',
    description='Ensemble Learning Approach for Stability Prediction of Interface and Core mutations',
    author='kimlab',
    author_email='elaspic@kimlab.org',
    url='http://elaspic.kimlab.org',
    packages=['elaspic', 'tests'],
    long_description=read("README.rst"),
    install_requires=[
        'fastcache',
        # Scientific Python stack
        'numpy',
        'pandas',
        'scikit-learn>=0.15',
        'sqlalchemy',
        'biopython>=1.65',
        # Documentation
        'sphinx>=1.3',
        # Testing
        'pytest',
        # Py2 / Py3 support
        'future',
        'six',
    ],
    tests_require=[
        'pytest',
    ],
    entry_points={
          'console_scripts': [
              'elaspic = elaspic.__main__:main',
              'elaspic_database = elaspic.elaspic_database:main',
          ]
      },
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Topic :: Bioinformatics",
    ],
)

### Template for future reference
#~ # Utility function to read the README file.
#~ # Used for the long_description.  It's nice, because now 1) we have a top level
#~ # README file and 2) it's easier to type in the README file than to put a raw
#~ # string in below ...
#~ def read(fname):
    #~ return open(os.path.join(os.path.dirname(__file__), fname)).read()
    #~ 
#~ setup(
    #~ name = "elaspic",
    #~ version = "0.98",
    #~ author = "temp",
    #~ author_email = "temp@gmail.com",
    #~ description = (
        #~ "Predict the effect of mutations using protein homology models " 
	#~ "and machine learning"),
    #~ license = "BSD",
    #~ keywords = "example documentation tutorial",
    #~ url = "http://packages.python.org/none_yet",
    #~ packages=['an_example_pypi_project', 'tests'],
    #~ long_description=read('README'),
    #~ classifiers=[
        #~ "Development Status :: 3 - Alpha",
        #~ "Topic :: Bioinformatics",
        #~ "License :: OSI Approved :: BSD License",
    #~ ],
#~ )
