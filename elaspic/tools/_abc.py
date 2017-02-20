import os
import os.path as op
import logging
import tempfile

from abc import ABC, abstractmethod

from kmtools import structure_tools, py_tools

logger = logging.getLogger(__name__)


class ToolError(Exception):
    pass


class MutatorError(ToolError):
    pass


class _Tool(ABC):

    def __init__(self):
        self.tempdir = op.join(tempfile.gettempdir(), self.__class__.__name__)
        os.makedirs(self.tempdir, exist_ok=True)
        # atexit.register(self.cleanup)


class _Mutator(_Tool):

    @property
    @abstractmethod
    def _result_slots(self):
        raise NotImplementedError

    @abstractmethod
    def __init__(self):
        self.result = py_tools.Struct(self._result_slots)
        self.cache = {}
        self.mutations = {}
        super().__init__()

    @abstractmethod
    def build(self):
        raise NotImplementedError

    @abstractmethod
    def mutate(self, mutation):
        raise NotImplementedError

    @property
    def done(self):
        return (all(key in self.result for key in self._result_slots) and
                all(op.isfile(key) for key in self._result_slots if key.endswith('file')))


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


class Analyzer(ABC):
    """Calculate structural properties for a PDB containing one or more chains.

    Runs the program POPS to calculate the interface size of the complexes
    This is done by calculating the surface of the complex and the seperated parts.
    The interface is then given by the substracting.
    """

    def __init__(self, pdb_file, working_dir, vdw_distance=5.0, min_contact_distance=4.0):
        self.pdb_file = pdb_file
        #: Folder with all the binaries (i.e. ./analyze_structure)
        self.working_dir = working_dir
        self.vdw_distance = vdw_distance
        self.min_contact_distance = min_contact_distance

        self._prepare_temp_folder(self.working_dir)

        self.sp = structure_tools.StructureParser()
        self.sp.extract(pdb_file)
        self.sp.save_structure(output_dir=self.working_dir)

        self.chain_ids = self.sp.chain_ids


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
