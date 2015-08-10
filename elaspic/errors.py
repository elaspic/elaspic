# -*- coding: utf-8 -*-

#%% Tcoffee

class TcoffeeError(Exception):
    def __init__(self, result, error, alignInFile, system_command):
        message = (
            'tcoffee blast error for file: {0}, with error message: {1}, '
            'when trying to run command: {2}'
            .format(alignInFile, error, system_command))
        Exception.__init__(self, message)
        self.result = result

class TcoffeeBlastError(Exception):
    def __init__(self, result, error, alignInFile, system_command):
        message = (
            'tcoffee blast error for file: {0}, with error message: {1}, '
            'when trying to run command: {2}'
            .format(alignInFile, error, system_command))
        Exception.__init__(self, message)
        self.result = result

class TcoffeePDBidError(Exception):
    def __init__(self, result, error, alignInFile, system_command):
        message = (
            'tcoffee blast error for file: {0}, with error message: {1}, '
            'when trying to run command: {2}'
            .format(alignInFile, error, system_command))
        Exception.__init__(self, message)
        self.result = result



#%% Provean

class ProveanError(Exception):
    pass

class ProveanResourceError(Exception):
    def __init__(self, message, child_process_group_id):
        Exception.__init__(self, message)
        self.child_process_group_id = child_process_group_id



#%% Finding templates (PDB in uppercase to be consistent with Biopython)

class LowIdentity(Exception):
    pass

class PDBError(Exception):
    pass

class PDBNotFoundError(Exception):
    pass

class PDBEmptySequenceError(Exception):
    """ One of the sequences is missing from the alignment. The most likely cause
    is that the alignment domain definitions were incorrect.
    """
    pass

class PDBDomainDefsError(Exception):
    """ PDB domain definitions not found in the pdb file
    """
    pass

class PDBChainError(Exception):
    pass


#%% Making models

class MSMSError(Exception):
    pass

class ModellerError(Exception):
    pass

class FoldxError(Exception):
    pass

class FoldXAAMismatchError(Exception):
    pass

class DataError(Exception):
    pass

class TemplateCoreError(Exception):
    pass

class TemplateInterfaceError(Exception):
    pass

class ResourceError(Exception):
    pass


#%% Computing mutations

class Archive7zipError(Exception):
    def __init__(self, result, error_message, return_code):
        super(Archive7zipError, self).__init__(result)
        self.error_message = error_message
        self.return_code = return_code

class Archive7zipFileNotFoundError(Archive7zipError):
    pass

class NoSequenceFound(Exception):
    pass

class ProteinDefinitionError(Exception):
    pass

class NoTemplatesFoundError(Exception):
    pass

class AlignmentNotFoundError(Exception):
    def __init__(self, save_path, alignment_filename):
        Exception.__init__(self)
        self.save_path = save_path
        self.alignment_filename = alignment_filename

class PopsError(Exception):
    def __init__(self, message, pdb, chains):
        Exception.__init__(self, message)
        self.pdb = pdb
        self.chains = chains


#%%

class ChainsNotInteractingError(Exception):
    pass

class MutationOutsideDomainError(Exception):
    pass

class MutationOutsideInterfaceError(Exception):
    pass

class NoModelFoundError(Exception):
    pass



#%% Database

class ModelHasMutationsError(Exception):
    """Don't delete a model that has precalculated mutations!
    """
    pass


