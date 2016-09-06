import os.path as op
from setuptools import setup, Command

try:
    from pypandoc import convert

    def read_md(file):
        return convert(file, 'rst')

except ImportError:
    print("warning: pypandoc module not found, could not convert Markdown to RST")

    def read_md(file):
        with open(op.join(op.dirname(__file__), file)) as ifh:
            return ifh.read()


class TrainPredictors(Command):
    user_options = []

    def initialize_options(self):
        """Abstract method that is required to be overwritten."""
        pass

    def run(self):
        """Train the ELASPIC classifier."""
        from elaspic.__main__ import elaspic_train
        elaspic_train(None)

    def finalize_options(self):
        """Abstract method that is required to be overwritten."""
        pass


setup(
    name='elaspic',
    version='0.1.37',
    description=(
        "Ensemble Learning Approach for Stability Prediction of "
        "Interface and Core mutations (ELASPIC)."),
    url="https://github.com/kimlaborg/elaspic",
    author='kimlab',
    author_email='alex.strokach@utoronto.ca',
    packages=['elaspic'],
    package_data={'elaspic': ['data/*']},
    long_description=read_md("README.md"),
    entry_points={
        'console_scripts': 'elaspic = elaspic.__main__:main'
    },
    classifiers=[
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    cmdclass={'train': TrainPredictors},
)
