import logging
import os
import os.path as op
import tempfile
from abc import ABC, abstractmethod

import Bio

from kmtools import py_tools, structure_tools

logger = logging.getLogger(__name__)


class ToolError(ABC, Exception):
    pass


class _Tool(ABC):

    @property
    @abstractmethod
    def _result_slots(self):
        raise NotImplementedError

    def __init__(self, *args, **kwargs):
        self.tempdir = op.join(tempfile.gettempdir(), self.__class__.__name__)
        os.makedirs(self.tempdir, exist_ok=True)
        self.result = py_tools.Struct(self._result_slots, *args, **kwargs)
        if args or kwargs:
            assert self.done

    def build(self):
        if not self.done:
            self._build()
            assert self.done

    @abstractmethod
    def _build(self):
        raise NotImplementedError

    @property
    def done(self):
        return (all(key in self.result for key in self._result_slots) and
                all(op.isfile(self.result[key]) for key in self._result_slots
                    if key.endswith('file')))


class _SequenceTool(_Tool):

    def __init__(self, sequence, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert len(sequence) > 0, "The sequence must not be empty!"
        self.sequence = sequence
        self._sequence_file = None

    @property
    def sequence_file(self):
        if self._sequence_file is None:
            self._sequence_file = op.join(self.tempdir, self.sequence.id + '.fasta')
            with open(self._sequence_file, 'wt') as ofh:
                Bio.SeqIO.write(self.sequence, ofh, 'fasta')
        return self._sequence_file


class _StructureTool(_Tool):

    def __init__(self, structure, *args, **kwargs):
        super().__init__(*args, **kwargs)
        assert len(structure) == 1, "The structure must have only one model!"
        self.structure = structure
        self._structure_file = None

    @property
    def structure_file(self):
        if self._structure_file is None:
            self._structure_file = op.join(self.tempdir, self.structure.id + '.pdb')
            structure_tools.save_structure(self.structure, self._structure_file)
        return self._structure_file


class Modeller(_StructureTool):

    _result_slots = ['homology_structure_file']

    @abstractmethod
    def model(self, alignment):
        raise NotImplementedError


class Mutator(_StructureTool):

    @abstractmethod
    def mutate(self, mutation):
        raise NotImplementedError


class SequenceAnalyzer(_SequenceTool):

    @abstractmethod
    def analyze(self, mutation):
        raise NotImplementedError


class StructureAnalyzer(_StructureTool):

    def analyze(self, mutation):
        assert self.done
        chain_id, chain_mutation = mutation.split('-')
        residue_id = (' ', int(chain_mutation[1:-1]), ' ')
        aa = chain_mutation[0]
        return self._analyze(chain_id, residue_id, aa)

    @abstractmethod
    def _analyze(self, chain_id, residue_id, aa):
        raise NotImplementedError
