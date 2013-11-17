# -*- coding: utf-8 -*-
"""
Created on Sun Feb  3 15:07:51 2013

@author: niklas
"""

class TcoffeeError(Exception):
    def __init__(self, message, errors, alignInFile):
        # Call the base class constructor with the parameters it needs
        Exception.__init__(self, message)

        # Now for your custom code...
        self.errors = errors
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