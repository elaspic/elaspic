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
        self.error = error
        self.alignInFile = alignInFile

class TcoffeeBlastError(Exception):
    def __init__(self, message, error, alignInFile):
        Exception.__init__(self, message)
        self.error = error
        self.alignInFile = alignInFile
        
class KNOTerror(Exception):
    pass

class ModellError(Exception):
    pass

class FoldXError(Exception):
    def __init__(self, error):
        # Call the base class constructor with the parameters it needs
        Exception.__init__(self)
        self.error = error

class DataError(Exception):
    def __init__(self, inputFile):
        # Call the base class constructor with the parameters it needs
        Exception.__init__(self)
        self.inputFile = inputFile

class ConfigError(Exception):
    def __init__(self, option):
        Exception.__init__(self)
        self.option = option

class TemplateCoreError(Exception):
    def __init__(self, error):
        Exception.__init__(self)
        self.error = error

class TemplateInterfaceError(Exception):
    def __init__(self, error):
        Exception.__init__(self)
        self.error = error

class pdbError(Exception):
    def __init__(self, error):
        Exception.__init__(self)
        self.error = error

class NoStructuralTemplates(Exception):
    def __init__(self, error):
        Exception.__init__(self)
        self.error = error

class NoSequenceFound(Exception):
    def __init__(self, error):
        Exception.__init__(self)
        self.error = error
        
class ProteinDefinitionError(Exception):
    def __init__(self, error):
        Exception.__init__(self)
        self.error = error
        
class NoTemplatesFound(Exception):
    def __init__(self, error):
        Exception.__init__(self)
        self.error = error
        
class EmptyPDBSequenceError(Exception):
    def __init__(self, pdb_id, pdb_chain):
        Exception.__init__(self)
        self.pdb_id = pdb_id
        self.pdb_chain = pdb_chain
        
class NoPrecalculatedAlignmentFound(Exception):
    def __init__(self, save_path, alignment_filename):
        Exception.__init__(self)
        self.save_path = save_path
        self.alignment_filename = alignment_filename
        
class MutationOutsideDomain(Exception):
    def __init__(self):
        Exception.__init__(self)

class NotInteracting(Exception):
    def __init__(self):
        Exception.__init__(self)







class PopsError(Exception):
    def __init__(self, e, pdb, chains):
        Exception.__init__(self)
        self.error = e
        self.pdb = pdb
        self.chains = chains
        


class NoPDBFound(Exception):
    def __init__(self, pdb_filename):
        Exception.__init__(self)
        self.error = 'PDB with filename %s not found!' % pdb_filename
