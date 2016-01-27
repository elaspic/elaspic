import os.path as op
from setuptools import setup, Command

import yaml


# %%
def read(fname):
    return open(op.join(op.dirname(__file__), fname)).read()


class TrainPredictors(Command):
    user_options = []

    def initialize_options(self):
        """Abstract method that is required to be overwritten"""
        pass

    def run(self):
        from elaspic.__main__ import elaspic_train
        elaspic_train(None)

    def finalize_options(self):
        """Abstract method that is required to be overwritten"""
        pass


# %% Load conda configuration file
with open('conda/elaspic/meta.yaml') as ifh:
    meta = yaml.load(ifh)


setup(
    name=meta['package']['name'],
    version=meta['package']['version'],
    description=meta['about']['summary'],
    author='kimlab',
    author_email='alex.strokach@utoronto.ca',
    url=meta['about']['home'],
    packages=['elaspic'],
    package_data={'elaspic': ['data/*']},
    long_description=read("README.rst"),
    # install_requires=meta['requirements']['run'],
    # tests_require=meta['test']['requires'],
    entry_points={'console_scripts': meta['build']['entry_points']},
    classifiers=[
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 3",
        "Topic :: Structural Biology",
        "Topic :: Bioinformatics",
    ],
    cmdclass={'train': TrainPredictors},
)
