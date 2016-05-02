import os.path as op
from setuptools import setup, Command
import yaml


with open('devtools/conda-recipe/meta.yaml', 'rb') as ifh:
    META = yaml.load(ifh)


def read(fname):
    """Read the contents of a file."""
    with open(op.join(op.dirname(__file__), fname)) as ifh:
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
    version='1.0.22',
    description=META.get('about', {}).get('summary', ''),
    url=META.get('about', {}).get('home', ''),
    author='kimlab',
    author_email='alex.strokach@utoronto.ca',
    packages=['elaspic'],
    package_data={'elaspic': ['data/*']},
    long_description=read("README.md"),
    entry_points={'console_scripts': META.get('build', {}).get('entry_points', '')},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Topic :: Structural Biology",
        "Topic :: Bioinformatics",
    ],
    cmdclass={'train': TrainPredictors},
)
