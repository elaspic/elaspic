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


# Making models
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
