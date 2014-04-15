# -*- coding: utf-8 -*-
"""
Created on Sun Feb  3 15:07:51 2013

@author: niklas
"""
import os
import sys

###############################################################################
# Used to find the location of the files being executed
def we_are_frozen():
    # All of the modules are built-in to the interpreter, e.g., by py2exe
    return hasattr(sys, "frozen")

def path_to_pipeline_code():
    encoding = sys.getfilesystemencoding()
    if we_are_frozen():
        return os.path.dirname(unicode(sys.executable, encoding))
    return os.path.dirname(unicode(__file__, encoding))
    
    
###############################################################################    
# Keep track of and raise different kinds of errors
    
class TcoffeeError(Exception):
    def __init__(self, message, error, alignInFile):
        # Call the base class constructor with the parameters it needs
        Exception.__init__(self, message)
        # Now for your custom code...
        self.error = 'tcoffee error for file: %s, with error message: %s' % (alignInFile, error)


class TcoffeeBlastError(Exception):
    def __init__(self, message, error, alignInFile):
        Exception.__init__(self, message)
        self.error = 'tcoffee blast error for file: %s, with error message: %s' % (alignInFile, error)
        
        
class TcoffeePDBidError(Exception):        
    def __init__(self, message, error, alignInFile):
        Exception.__init__(self, message)
        self.error = 'tcoffee pdbid error for file: %s, with error message: %s' % (alignInFile, error)
        


###############################################################################
# Common

class ProveanError(Exception):
    def __init__(self, message):
        Exception.__init__(self, 'provean exited with an error:\n %s' % message)


###############################################################################
# Finding templates

class pdbError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)
        
        
class EmptyPDBSequenceError(Exception):
    def __init__(self, pdb_id, pdb_chain):
        message = 'Empty pdb sequence file for pdb: %s, chain: %s' % (pdb_id, pdb_chain)
        Exception.__init__(self, message)
        self.pdb_id = pdb_id
        self.pdb_chain = pdb_chain

class PDBDomainDefsError(Exception):
    # domain.pdb_domain_defs not found in pdb file
    def __init__(self, message):
        Exception.__init__(self, message)

###############################################################################
# Making models

class MSMSError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)


class ModellerFailure(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)
        
        
class FoldXError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)


class DataError(Exception):
    def __init__(self, inputFile):
        Exception.__init__(self)
        self.inputFile = inputFile


class TemplateCoreError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)


class TemplateInterfaceError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)


###############################################################################
# Computing mutations

class PDBChainError(Exception):
    def __init__(self, pdb_code, chains):
        message = 'PDBChainError in pdb: %s and chain: %s' % (pdb_code, chains,)
        Exception.__init__(self, message)


class NoStructuralTemplates(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)


class NoSequenceFound(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)
        
        
class ProteinDefinitionError(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)

        
class NoTemplatesFound(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)

        
class NoPrecalculatedAlignmentFound(Exception):
    def __init__(self, save_path, alignment_filename):
        Exception.__init__(self)
        self.save_path = save_path
        self.alignment_filename = alignment_filename


class PopsError(Exception):
    def __init__(self, message, pdb, chains):
        Exception.__init__(self, message)
        self.pdb = pdb
        self.chains = chains
   
     
class NoPDBFound(Exception):
    def __init__(self, pdb_filename):
        message = 'PDB with filename %s not found!' % pdb_filename
        Exception.__init__(self, message)


class NoDomainFound(Exception):
    def __init__(self, pdb_filename):
        message = 'PDB with filename %s not found!' % pdb_filename
        Exception.__init__(self, message)


###############################################################################

class ChainsNotInteracting(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)

     
class MutationOutsideDomain(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)


class MutationOutsideInterface(Exception):
    def __init__(self, message):
        Exception.__init__(self, message)

