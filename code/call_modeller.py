# Homology modeling by the automodel class
from modeller import *			# Load standard Modeller classes
from modeller.automodel import *	# Load the automodel class
import helper_functions as hf
import errors 

class modeller:
    """
    run modeller

    input:  alignment                   type: list      list containing filenames with modeller
                                                        input files in PIR format
            seqID                       type: string    name of the sequence (target)
            templateID                  type: string    name of the template (structure)
            path_to_pdb_for_modeller    type: string    path to PDB files
            tmpPath                     type: string    path for storing tmp files
            modeller_runs               type: int       how many rounds of modelling should be done
            loopRefinement              type: boolean   if True calculate loop refinemnts
    """
    def __init__(
            self, alignment, seqID, templateID, path_to_pdb_for_modeller,
            tmpPath, logger, modeller_runs, loopRefinement=True,):

        if not isinstance(alignment, list):
            self.alignment = list()
            self.alignment.append(alignment)
        else:
            self.alignment = alignment
        self.seqID = seqID
        self.templateID = templateID
        self.loopRefinement = loopRefinement
        self.start = 1                # start model
        self.end = modeller_runs      # end model
        self.loopStart = 1            # start loop refinement model
        self.loopEnd = 4              # end loop refinement model

        # some environment settings
        self.filePath = path_to_pdb_for_modeller
        self.tmpPath = tmpPath
        self.logger = logger


    def run(self):
        """
        """
        ranking = dict() # key: assessment score
                         # values: alignment, pdb of modell, from loop refinement yes or no
        ranking_knotted = dict()

        knotted = True
        counter = 0
        while knotted and counter < 10:
            counter += 1
            ## NB: You can actually supply many alignments and modeller will give
            # you the alignment with the best model, and that model
            for aln in self.alignment:
                # There is a chance that modeller fails to automatically select the
                # regions for the loop modelling. In that case it tries to model
                # without loop refinement. Worst thing that could happen is that
                # the model is bad
                try:
                    result, loop, failures = self.__run_modeller(aln, self.loopRefinement)
                except IndexError:
                    raise
                except:
                    try:
                        result, loop, failures = self.__run_modeller(aln, False)
                    except ModellerError as e:
                        raise errors.ModellerError(str(e))
                        
                if not result:
                    raise errors.ModellerError(str(failures[-1]))

                for i in range(len(result)):
                    pdbFile, normDOPE = result[i][0], result[i][1]
                    if self.__call_knot(pdbFile):  # i.e. knotted
                        ranking_knotted[normDOPE] = (aln, str(pdbFile), loop,)
                    else:
                        ranking[normDOPE] = (aln, str(pdbFile), loop,)
                        knotted = False
        if not knotted:
            return min(ranking), ranking[min(ranking)][1], knotted
        else:
            return min(ranking_knotted), ranking_knotted[min(ranking_knotted)][1], knotted


    def __make_alignment(self):
        """
        Functionality for modeller to make the alignment instead of relying on tcoffee
        """
        pass


    def __run_modeller(self, alignFile, loopRefinement):
        """
        more or less taken from the modeller website

        input:  alignFile   type: string        File containing the input data
                result      type: list          The successfully calculated models
                                                are stored in this list
                loopRefinement  type: boolean   if True: perform loopfeinements

        return: result      type: list          successfully calculated models
        """
        log.none() # instructs Modeller to display no log output.
        env = environ() # create a new MODELLER environment to build this model in

        # directories for input atom files
        print 'atom_files_directory: {}: {}'.format(type(self.filePath), self.filePath)
        env.io.atom_files_directory = [str(self.filePath.rstrip('/'))]
        env.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)
        # selected atoms do not feel the neighborhood
        #env.edat.nonbonded_sel_atoms = 2

        # Read in HETATM records from template PDBs
        env.io.hetatm = True

        if loopRefinement == False:
            a = automodel(env,
                          alnfile = alignFile,      # alignment filename
                          knowns = (self.templateID),                # codes of the templates
                          sequence = self.seqID,	            		  # code of the target
                          assess_methods=(assess.DOPE,
                                          assess.normalized_dope)    # wich method for validation should be calculated
                          )

        if loopRefinement == True:
            a = dope_loopmodel(env,
                               alnfile = alignFile,      # alignment filename
                               knowns = (self.templateID),                # codes of the templates
                               sequence = self.seqID,	            		  # code of the target
                               assess_methods=(assess.DOPE,
                                               assess.normalized_dope),   # wich method for validation should be calculated
                               loop_assess_methods=(assess.DOPE,
                                                    assess.normalized_dope)
                               )

        a.starting_model = self.start		# index of the first model
        a.ending_model = self.end   		# index of the last model
                							# (determines how many models to calculate)
        # Very thorough VTFM optimization:
        a.library_schedule = autosched.slow
        a.max_var_iterations = 300

        # Thorough MD optimization:
#        a.md_level = refine.slow
        a.md_level = None

        # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
#        a.repeat_optimization = 2

        if loopRefinement == True:
            a.loop.starting_model = self.loopStart      # index of the first loop model
            a.loop.ending_model   = self.loopEnd        # index of the last loop model
            a.loop.md_level       = refine.slow         # loop refinement method; this yields

        a.max_molpdf = 2e5

        a.make()						# do the actual homology modeling

        # the output produced by modeller is stored in a.loop.outputs or a.outputs
        # it is a dictionary
        #
        # check for each model if it was successfully calculated, i.e.
        # for each "normal" model and each loop model and append the
        # assessment score to a list which is used to return the best model
        loop = False
        result = []
        failures = []
        # add the loop output
        if loopRefinement == True:
            self.logger.debug('Modeller loop outputs:')
            self.logger.debug(a.loop.outputs)
            for i in range(len(a.loop.outputs)):
                if not a.loop.outputs[i]['failure']:
                    model_filename = a.loop.outputs[i]['name']
                    model_dope_score = a.loop.outputs[i]['Normalized DOPE score']
                    result.append((model_filename, model_dope_score))
                    loop = True
                else:
                    failures.append(a.loop.outputs[i]['failure'])

        # add the normal output
        self.logger.debug('Modeller outputs:')
        self.logger.debug(a.outputs)
        for i in range(len(a.outputs)):
            if not a.outputs[i]['failure']:
                model_filename = a.outputs[i]['name']
                model_dope_score = a.outputs[i]['Normalized DOPE score']
                result.append((model_filename, model_dope_score))
            else:
                failures.append(a.outputs[i]['failure'])

        # return the successfully calculated models and a loop flag indicating
        # whether the returned models are loop refined or not
        return result, loop, failures


    def __call_knot(self, pdbFile):
        """
        check a PDB structure for knots using the program KNOTS by Willi Taylor

        Make sure to set the command and PDB path!

        input:  pdbFile     type: string        Filename of the PDB structure

        return: 0 or 1      type: int           0: no knot
                                                1: knotted
        """
        # for multiprocessing the KNOT program needs to be run from a unique
        # execution folder, this is created in the beginning

        system_command = './topol ' + pdbFile
        self.logger.debug('FoldX system command: {}'.format(system_command))
        child_process = hf.run_subprocess_locally('./', system_command)
        result, error_message = child_process.communicate()
        return_code = child_process.returncode
        if not return_code:
            self.logger.error('FoldX result: {}'.format(result))
            self.logger.error('FoldX error message: {}'.format(error_message))

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

