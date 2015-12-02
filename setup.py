import os.path as op
import json
import pickle
from setuptools import setup, Command

import yaml
import pandas as pd
# import distutils.command.bdist_conda


# %%
def read(fname):
    """Utility function to read the README file.

    Used for the long_description.  It's nice, because now
      1) we have a top level README file and
      2) it's easier to type in the README file than to put a raw string in below ...

    Source: https://pythonhosted.org/an_example_pypi_project/setuptools.html
    """
    return open(op.join(op.dirname(__file__), fname)).read()


class TrainPredictors(Command):
    user_options = []

    def initialize_options(self):
        """Abstract method that is required to be overwritten"""
        pass

    def run(self):
        from elaspic import DATA_DIR, machine_learning
        SETUP_DIR = op.abspath(op.dirname(__file__))

        # Load data
        with open(op.join(SETUP_DIR, 'training_data', 'core_options_p0.json')) as ifh:
            core_options_p0 = json.load(ifh)
        with open(op.join(SETUP_DIR, 'training_data', 'core_options_p1.json')) as ifh:
            core_options_p1 = json.load(ifh)
        with open(op.join(SETUP_DIR, 'training_data', 'interface_options_p0.json')) as ifh:
            interface_options_p0 = json.load(ifh)
        with open(op.join(SETUP_DIR, 'training_data', 'interface_options_p1.json')) as ifh:
            interface_options_p1 = json.load(ifh)

        core_training_set = (
            pd.read_csv(
                op.join(SETUP_DIR, 'training_data', 'core_training_set.tsv.gz'), sep='\t')
        )
        interface_training_set = (
            pd.read_csv(
                op.join(SETUP_DIR, 'training_data', 'interface_training_set.tsv.gz'), sep='\t')
        )

        # Train predictors
        core_clf_p0 = (
            machine_learning.get_final_predictor(
                core_training_set, core_options_p0['features'], core_options_p0)
        )
        core_clf_p1 = (
            machine_learning.get_final_predictor(
                core_training_set, core_options_p1['features'], core_options_p1)
        )
        interface_clf_p0 = (
            machine_learning.get_final_predictor(
                interface_training_set, interface_options_p0['features'], interface_options_p0)
        )
        interface_clf_p1 = (
            machine_learning.get_final_predictor(
                interface_training_set, interface_options_p1['features'], interface_options_p1)
        )

        # Save predictors and features
        with open(op.join(DATA_DIR, 'ml_features_core_p0.json'), 'w') as ofh:
            json.dump(core_options_p0['features'], ofh)
        with open(op.join(DATA_DIR, 'ml_features_core_p1.json'), 'w') as ofh:
            json.dump(core_options_p1['features'], ofh)
        with open(op.join(DATA_DIR, 'ml_features_interface_p0.json'), 'w') as ofh:
            json.dump(interface_options_p0['features'], ofh)
        with open(op.join(DATA_DIR, 'ml_features_interface_p1.json'), 'w') as ofh:
            json.dump(interface_options_p1['features'], ofh)

        with open(op.join(DATA_DIR, 'ml_clf_core_p0.pickle'), 'wb') as ofh:
            pickle.dump(core_clf_p0, ofh)
        with open(op.join(DATA_DIR, 'ml_clf_core_p1.pickle'), 'wb') as ofh:
            pickle.dump(core_clf_p1, ofh)
        with open(op.join(DATA_DIR, 'ml_clf_interface_p0.pickle'), 'wb') as ofh:
            pickle.dump(interface_clf_p0, ofh)
        with open(op.join(DATA_DIR, 'ml_clf_interface_p1.pickle'), 'wb') as ofh:
            pickle.dump(interface_clf_p1, ofh)

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
