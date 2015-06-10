import os
from setuptools import setup



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
    version='0.1.2',
    description='Ensemble Learning Approach for Stability Prediction of Interface and Core mutations',
    author='kimlab',
    author_email='elaspic@kimlab.org',
    url='http://elaspic.kimlab.org',
    packages=['elaspic'],
    package_data={'elaspic': ['data/*']},
    long_description=read("README.rst"),
    # Specify install requirements in the conda `meta.yaml` file
    # install_requires=[],
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
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Topic :: Structural Biology",
        "Topic :: Bioinformatics",
    ],
)

