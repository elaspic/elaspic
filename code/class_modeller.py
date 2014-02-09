# Homology modeling by the automodel class
from modeller import *			# Load standard Modeller classes
from modeller.automodel import *	# Load the automodel class

import subprocess, shlex


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
    def __init__(self, 
                 alignment, 
                 seqID, 
                 templateID, 
                 path_to_pdb_for_modeller,
                 tmpPath, 
                 modeller_runs,
                 loopRefinement=True, 
                 ):
        # input
        # make sure that the alignment is contained in a list
        if not isinstance(alignment, list):
            self.alignment = list()
            self.alignment.append(alignment)
        else:
            self.alignment = alignment
        
        self.seqID = seqID
        self.templateID = templateID
        
        # some behaviour
        self.loopRefinement = loopRefinement

        self.start = 1                # start model
        self.end = modeller_runs      # end model
        self.loopStart = 1            # start loop refinement model
        self.loopEnd = 4              # end loop refinement model
        
        # some environment settings
        self.filePath = path_to_pdb_for_modeller
        self.tmpPath = tmpPath
        self.modeller_path = tmpPath + 'modeller/'
        

    def run(self):
        """
        """
        result = list()
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
                    result, loop = self.__run_modeller(aln, result, self.loopRefinement)
                except IndexError:
                    raise
                except:
                    result, loop = self.__run_modeller(aln, result, False)
                    
                for i in range(len(result)):
                    pdbFile, normDOPE = result[i][0], result[i][1]
                    if self.__call_knot(pdbFile, self.modeller_path):  # i.e. knotted
                        ranking_knotted[normDOPE] = (aln, str(pdbFile), loop,)
                    else: 
                        ranking[normDOPE] = (aln, str(pdbFile), loop,)
                        knotted = False                            

        if not knotted:
            return min(ranking), ranking[min(ranking)][1], knotted
        else:
            return min(ranking_knotted), ranking_knotted[min(ranking_knotted)][1], knotted




    def __run_modeller(self, alignFile, result, loopRefinement):
        """
        more or less taken from the modeller website
        
        input:  alignFile   type: string        File containing the input data
                result      type: list          The successfully calculated models
                                                are stored in this list
                loopRefinement  type: boolean   if True: perform loopfeinements
                
        return: result      type: list          successfully calculated models
        """
        log.none()	# request none output, set to verbose or minimal if needed
#        log.minimal()
#        log.verbose()
        
        env = environ()
        # create a new MODELLER environment to build this model in
        
        # directories for input atom files
        env.io.atom_files_directory = [self.filePath]
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
    
        a.max_molpdf = 1e6
        
        a.make()						# do the actual homology modeling
        

        
        # the output produced by modeller is stored in a.loop.outputs or a.outputs
        # it is a dictionary
        #
        # check for each model if it was successfully calculated, i.e.
        # for each "normal" model and each loop model and append the
        # assessment score to a list which is used to return the best model
        loop = False
        # add the loop output
        if loopRefinement == True:
            for i in range(len(a.loop.outputs)):
                if a.loop.outputs[i]['failure'] == None:
                    result.append((a.loop.outputs[i]['name'], a.loop.outputs[i]['Normalized DOPE score']))
                    loop = True
        # add the normal output
        for i in range(len(a.outputs)):
            result.append((a.outputs[i]['name'], a.outputs[i]['Normalized DOPE score']))
      
        
        # return the successfully calculated models and a loop flag indicating
        # whether the returned models are loop refined or not
        return result, loop



    def __call_knot(self, pdbFile, modeller_path):
        """
        check a PDB structure for knots using the program KNOTS by Willi Taylor
        
        Make sure to set the command and PDB path!
        
        input:  pdbFile     type: string        Filename of the PDB structure
        
        return: 0 or 1      type: int           0: no knot
                                                1: knotted
        """
        # for multiprocessing the KNOT program needs to be run from a unique
        # execution folder, this is created in the beginning

        system_command = modeller_path + 'topol ' + modeller_path + pdbFile

        cmd = shlex.split(system_command)
        childProcess = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=False)
        output, error = childProcess.communicate()
        rc = childProcess.returncode
    
        line = [ x for x in output.split('\n') ]
        
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




