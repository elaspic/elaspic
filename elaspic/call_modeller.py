# -*- coding: utf-8 -*-
"""
Homology modeling by the automodel class
"""
from __future__ import absolute_import
from __future__ import unicode_literals
from builtins import range
from builtins import object

import six
import logging

from modeller import * # Load standard Modeller classes
from modeller.automodel import *	# Load the automodel class

from . import helper_functions as hf
from . import errors
from .conf import configs

logger = logging.getLogger(__name__)

class Modeller(object):
    """
    Runs MODELLER in order to make a homology model of the given protein.

    Parameters
    ----------
    alignment : list
        Contains filenames with modeller input files in a PIR file format
    seqID : string
        Name of the sequence (target)
    templateID : string
        Name of the template (structure)
    filePath : string
        Path to PDB files
    tmpPath : string
        Path for storing tmp files
    modeller_runs : int
       How many rounds of modelling should be done
    loopRefinement : boolean
        If True, calculate loop refinemnts
    """
    def __init__(self, alignment, seqID, templateID, filePath, loopRefinement=True):

        if not isinstance(alignment, list):
            self.alignment = [alignment]
        else:
            self.alignment = alignment
        self.seqID = seqID
        self.templateID = templateID
        self.loopRefinement = loopRefinement
        self.start = 1                # start model
        self.end = configs['modeller_runs']      # end model
        self.loopStart = 1            # start loop refinement model
        self.loopEnd = 4              # end loop refinement model

        # Some environment settings
        self.filePath = filePath
        self.tmpPath = configs['unique_temp_folder']


    def run(self):
        """
        """
        ranking = dict() # key: assessment score
                         # values: alignment, pdb filename, whether or not using loop refinement
        ranking_knotted = dict()
        knotted = True
        counter = 0
        while knotted and counter < 10:
            counter += 1
            ## NB: You can actually supply many alignments and modeller will give
            # you the alignment with the best model
            for aln in self.alignment:
                # There is a chance that modeller fails to automatically select the
                # regions for the loop modelling. In that case it tries to model
                # without loop refinement. Worst thing that could happen is that
                # the model is bad
                try:
                    result, loop, failures = self.__run_modeller(aln, self.loopRefinement)
                except Exception as e:
                    logger.error('Loop refinement failed with an error: {}'.format(e))
                    try:
                        result, loop, failures = self.__run_modeller(aln, False)
                    except ModellerError as e:
                        raise errors.ModellerError(e)
                    except Exception as e:
                        raise e
                if not result:
                    raise errors.ModellerError(failures[-1])

                for i in range(len(result)):
                    pdbFile, normDOPE = result[i][0], result[i][1]
                    if self.__call_knot(pdbFile):  # i.e. knotted
                        ranking_knotted[normDOPE] = (aln, pdbFile, loop,)
                    else:
                        ranking[normDOPE] = (aln, pdbFile, loop,)
                        knotted = False
        if not knotted:
            return min(ranking), ranking[min(ranking)][1], knotted
        else:
            return min(ranking_knotted), ranking_knotted[min(ranking_knotted)][1], knotted


    def __run_modeller(self, alignFile, loopRefinement):
        """
        Parameters
        ----------
        alignFile : string
            File containing the input data
        result : list
            The successfully calculated models are stored in this list
        loopRefinement : boolean
            If `True`, perform loop refinements

        Returns
        -------
        list
            Successfully calculated models
        """

        log.none() # instructs Modeller to display no log output.
        env = environ() # create a new MODELLER environment to build this model in

        # Directories for input atom files
        env.io.atom_files_directory = [str(self.filePath.rstrip('/')),]
        env.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)

        # Selected atoms do not feel the neighborhood
        #env.edat.nonbonded_sel_atoms = 2
        env.io.hetatm = True # read in HETATM records from template PDBs
        env.io.water = True # read in WATER records (including waters marked as HETATMs)

        logger.debug('Performing loop refinement in addition to regular modelling: {}'.format(loopRefinement))
        if loopRefinement == False:
            a = automodel(
                env,
                alnfile=str(alignFile), # alignment filename
                knowns=(str(self.templateID)), # codes of the templates
                sequence=str(self.seqID), # code of the target
                assess_methods=(assess.DOPE, assess.normalized_dope) # wich method for validation should be calculated
            )
        else:
            a = dope_loopmodel(
                env,
                alnfile=str(alignFile), # alignment filename
                knowns=(str(self.templateID)), # codes of the templates
                sequence=str(self.seqID), # code of the target
                assess_methods=(assess.DOPE, assess.normalized_dope), # wich method for validation should be calculated
                loop_assess_methods=(assess.DOPE, assess.normalized_dope)
            )
            a.loop.starting_model = self.loopStart # index of the first loop model
            a.loop.ending_model = self.loopEnd # index of the last loop model
            a.loop.md_level = refine.slow # loop refinement method; this yields


        a.starting_model = self.start # index of the first model
        a.ending_model = self.end # index of the last model

        # Very thorough VTFM optimization:
        a.library_schedule = autosched.slow
        a.max_var_iterations = 300

        # Thorough MD optimization:
#        a.md_level = refine.slow
        a.md_level = None

        # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
#        a.repeat_optimization = 2

        a.max_molpdf = 2e5

        with hf.log_print_statements(logger):
            a.make() # do the actual homology modeling

        # The output produced by modeller is stored in a.loop.outputs or a.outputs
        # it is a dictionary
        # Check for each model if it was successfully calculated, i.e.
        # for each "normal" model and each loop model and append the
        # assessment score to a list which is used to return the best model
        result = []
        loop = False
        failures = []
        # Add the normal output
        for i in range(len(a.outputs)):
            if not a.outputs[i]['failure']:
                model_filename = a.outputs[i]['name']
                model_dope_score = a.outputs[i]['Normalized DOPE score']
                logger.debug('Success! model_filename: {}, model_dope_score: {}'.format(model_filename, model_dope_score))
                result.append((model_filename, model_dope_score))
            else:
                failure = a.outputs[i]['failure']
                logger.debug('Failure! {}'.format(failure))
                failures.append(a.outputs[i]['failure'])

        # Add the loop refinement output
        if loopRefinement:
            logger.debug('Modeller loop outputs:')
            for i in range(len(a.loop.outputs)):
                if not a.loop.outputs[i]['failure']:
                    model_filename = a.loop.outputs[i]['name']
                    model_dope_score = a.loop.outputs[i]['Normalized DOPE score']
                    logger.debug('Success! model_filename: {}, model_dope_score: {}'.format(model_filename, model_dope_score))
                    result.append((model_filename, model_dope_score))
                    loop = True
                else:
                    failure = a.loop.outputs[i]['failure']
                    logger.debug('Failure! {}'.format(failure))
                    failures.append(failure)

        # Return the successfully calculated models and a loop flag indicating
        # whether the returned models are loop refined or not
        return result, loop, failures


    def __call_knot(self, pdbFile):
        """
        Check a PDB structure for knots using the program KNOTS by Willi Taylor
        Make sure to set the command and PDB path!

        N.B.
        For multiprocessing the KNOT program needs to be run from a unique execution folder.

        Parameters
        ----------
        pdbFile : string
            Filename of the PDB structure

        Returns
        -------
        int
            0: no knot; 1: knot
        """


        system_command = 'knot ' + pdbFile
        logger.debug('Knot system command: {}'.format(system_command))
        child_process = hf.run_subprocess_locally('./', system_command)
        result, error_message = child_process.communicate()
        if six.PY3:
            result = str(result, encoding='utf-8')
            error_message = str(error_message, encoding='utf-8')
        return_code = child_process.returncode
        if not return_code:
            logger.error('Knot result: {}'.format(result))
            logger.error('Knot error message: {}'.format(error_message))

        line = [ x for x in result.split('\n') ]

        # I found two different forms in which the output appears,
        # hence two if statements to catch them
        if line[-4].strip().split(' ')[0] == 'len':
            if int(line[-4].strip().split(' ')[2]) == 2:
                return False # i.e. no knot
            if int(line[-4].strip().split(' ')[2]) > 2:
                return True # i.e. knotted
        elif line[-2].strip().split(' ')[0] == 'len':
            if int(line[-2].strip().split(' ')[2]) == 2:
                return False # i.e. no knot
            if int(line[-2].strip().split(' ')[2]) > 2:
                return True # i.e. knotted
        else:
            # in case the output can't be read, the model is classified as
            # knotted and thus disregarded. This could be improved.
            return True

