import os.path as op
from setuptools import setup, find_packages, Command


def _read_md_as_rst(file):
    """Read MarkDown file and convert it to ReStructuredText."""
    from pypandoc import convert
    return convert(file, 'rst')


def _read_md_as_md(file):
    """Read MarkDown file."""
    with open(op.join(op.dirname(__file__), file)) as ifh:
        return ifh.read()


def read_md(file):
    """Read MarkDown file and try to convert it to ReStructuredText if you can."""
    try:
        return _read_md_as_rst(file)
    except ImportError:
        print("WARNING: pypandoc module not found, could not convert Markdown to RST!")
        return _read_md_as_md(file)


class TrainPredictors(Command):
    user_options = []

    def initialize_options(self):
        """Abstract method that is required to be overwritten."""
        pass

    def run(self):
        """Train the ELASPIC classifier."""
        from elaspic.cli.elaspic_train import elaspic_train
        elaspic_train(None)

    def finalize_options(self):
        """Abstract method that is required to be overwritten."""
        pass


setup(
    name='elaspic',
    version='0.2.0',
    author='kimlab',
    author_email='alex.strokach@utoronto.ca',
    url="https://github.com/kimlaborg/elaspic",
    description=(
        "Ensemble Learning Approach for Stability Prediction of "
        "Interface and Core mutations (ELASPIC)."),
    long_description=read_md("README.md"),
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3 :: Only",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
    license='MIT',
    packages=find_packages(),
    package_data={
        'elaspic.database': ['sql/*.sql'],
        'elaspic.predictor': ['data/*'],
        'elaspic.tools': ['foldx/*'],
    },
    entry_points={
        'console_scripts': 'elaspic = elaspic.__main__:main'
    },
    cmdclass={'train': TrainPredictors},
)
