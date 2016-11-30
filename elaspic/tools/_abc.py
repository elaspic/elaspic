import os
import os.path as op
import tempfile

from abc import ABC, abstractmethod


class Modeller(ABC):

    def __init__(self, sequence, structure, alignment):
        self.sequence = sequence
        self.structure = structure
        self.alignment = alignment
        self.tempdir = op.join(tempfile.gettempdir(), self.__class__.__name__)
        os.makedires(self.tempdir, exist_ok=True)

    @abstractmethod
    def model(self):
        raise NotImplementedError


class Mutator(ABC):

    def __init__(self, structure, mutation):
        self.structure = structure
        self.mutation = mutation
        self.tempdir = op.join(tempfile.gettempdir(), self.__class__.__name__)
        os.makedires(self.tempdir, exist_ok=True)

    @abstractmethod
    def mutate(self):
        """
        Returns
        -------
        (str, str) : Filename of the wild type and mutant structures, respectively.
        """
        raise NotImplementedError


class CoreAnalyzer(ABC):

    def __init__(self, structure, chain, residue):
        self.structure = structure
        self.chain = chain
        self.residue = residue
        self.tempdir = op.join(tempfile.gettempdir(), self.__class__.__name__)
        os.makedires(self.tempdir, exist_ok=True)

    def analyze(self):
        """
        Returns
        -------
        dict : Dictionary containing all the calculated features.
        """
        raise NotImplementedError


class InterfaceAnalyzer(ABC):

    def __init__(self, structure, chain_1, residue_1, chain_2, residue_2):
        self.structure = structure
        self.chain_1 = chain_1
        self.residue_1 = residue_1
        self.chain_2 = chain_2
        self.residue_2 = residue_2
        self.tempdir = op.join(tempfile.gettempdir(), self.__class__.__name__)
        os.makedires(self.tempdir, exist_ok=True)

    def analyze(self):
        """
        Returns
        -------
        dict : Dictionary containing all the calculated features.
        """
        raise NotImplementedError
