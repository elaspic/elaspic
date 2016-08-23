# __main__
class ParameterError(Exception):
    pass


# T-Coffee
class TcoffeeError(Exception):

    def __init__(self, result, error, alignInFile, system_command):
        message = (
            'tcoffee blast error for file: {0}, with error message: {1}, '
            'when trying to run command: {2}'
            .format(alignInFile, error, system_command))
        Exception.__init__(self, message)
        self.result = result


# Provean
class ProveanError(Exception):
    pass


class ProveanResourceError(Exception):
    def __init__(self, message, child_process_group_id):
        Exception.__init__(self, message)
        self.child_process_group_id = child_process_group_id


class MutationMismatchError(Exception):
    pass


# Finding templates (PDB in uppercase to be consistent with Biopython)
class PDBError(Exception):
    pass


class PDBNotFoundError(Exception):
    pass


class PDBEmptySequenceError(Exception):
    """One of the sequences is missing from the alignment.

    The most likely cause is that the alignment domain definitions were incorrect.
    """
    pass


class PDBDomainDefsError(Exception):
    """PDB domain definitions not found in the pdb file."""
    pass


class PDBChainError(Exception):
    pass


# Making models
class MSMSError(Exception):
    pass


class ModellerError(Exception):
    pass


class FoldxError(Exception):
    pass


class FoldXAAMismatchError(Exception):
    pass


class ResourceError(Exception):
    pass


class InterfaceMismatchError(Exception):
    pass


# Computing mutations
class NoSequenceFound(Exception):
    pass


class ChainsNotInteractingError(Exception):
    pass


class MutationOutsideDomainError(Exception):
    pass


class MutationOutsideInterfaceError(Exception):
    pass


# Database
class Archive7zipError(Exception):
    def __init__(self, result, error_message, return_code):
        super(Archive7zipError, self).__init__(result)
        self.error_message = error_message
        self.return_code = return_code


class Archive7zipFileNotFoundError(Archive7zipError):
    pass
