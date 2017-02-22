"""Homology modeling by the automodel class.


Can use either dunamic_sphere or dynamic_lennard contraints
https://salilab.org/modeller/9.17/manual/node108.html


selection.randomize_xyz

"""
import logging
import os.path as op

import Bio
import elaspic
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from elaspic.tools._abc import Modeller, Mutator, ToolError
from modeller import ModellerError, environ, log, physical
from modeller.automodel import (assess, automodel, autosched, dope_loopmodel,
                                refine)

from kmtools import py_tools

logger = logging.getLogger(__name__)
ToolError.register(ModellerError)


class SaliLabModeller(Modeller):

    _result_slots = ['homology_structure_file', 'norm_dope_score']

    def _build(self):
        self.result['homology_structure_file'], self.result['norm_dope_score'] = (
            _build_modeller(self.alignment, self.tempdir)
        )


class SaliLabMutator(Mutator):

    _result_slots = ['wt_structure_file', 'mut_structure_file', 'norm_dope_score']

    def _build(self):
        template_id = self.structure.id + ''.join(chain.id for chain in self.structure[0])
        template_sequence = '/'.join(chain.sequence for chain in self.structure[0])

        target_id = template_id + '-' + self.mutation
        target_sequence = []
        for chain in self.structure[0]:
            if chain.id == self.mutation.split('-')[0]:
                assert chain.sequence[:int(self.mutation[1:-1])] == self.mutation.split('-')[1][0]
                target_sequence.append(
                    chain.sequence[:int(self.mutation[1:-1])] +
                    self.mutation[-1] +
                    chain.sequence[int(self.mutation[1:-1]) + 1:])
            else:
                target_sequence.append(chain.sequence)
        target_sequence = '/'.join(target_sequence)

        self.alignment = MultipleSeqAlignment([
            SeqRecord(id=template_id, seq=Seq(template_sequence)),
            SeqRecord(id=target_id, seq=Seq(target_sequence)),
        ])

        self.result['wt_structure_file'] = self.structure_file
        self.result['mut_structure_file'], self.result['norm_dope_score'] = (
            _build_modeller(self.alignment, self.tempdir)
        )


def _build_modeller(alignment, tempdir):
        results = []
        counter = 0
        max_counter = 3
        while not results or counter < max_counter:
            logger.debug("counter: %s, ranking: '%s'", counter, results)
            counter += 1
            # NB: You can actually supply many alignments and modeller will give
            # you the alignment with the best model

            # There is a chance that modeller fails to automatically select the
            # regions for the loop modelling. In that case it tries to model
            # without loop refinement. Worst thing that could happen is that
            # the model is bad
            try:
                result, loop, failures = _run_modeller(alignment, tempdir, True)
            except Exception as e:
                logger.error('Loop refinement failed with an error: {}'.format(e))
                try:
                    result, loop, failures = _run_modeller(alignment, tempdir, False)
                except ModellerError as e:
                    raise elaspic.exc.ModellerError(e)
                except Exception as e:
                    raise e
            if not result:
                raise elaspic.exc.ModellerError(failures[-1])
            results.extend(result)
        results.sort(key=lambda x: x[1])
        return results[0][0], results[0][1]


def _run_modeller(alignment, tempdir, loop_refinement):
    """.

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
    alignment_file = _write_alignment(alignment, tempdir)

    log.none()  # instructs Modeller to display no log output.
    env = environ()  # create a new MODELLER environment to build this model in

    # Directories for input atom files
    env.io.atom_files_directory = [tempdir, ]
    env.schedule_scale = physical.values(default=1.0, soft_sphere=0.7)

    # Selected atoms do not feel the neighborhood
    # env.edat.nonbonded_sel_atoms = 2
    env.io.hetatm = True  # read in HETATM records from template PDBs
    env.io.water = True  # read in WATER records (including waters marked as HETATMs)

    logger.debug(
        'Performing loop refinement in addition to regular modelling: {}'
        .format(loop_refinement)
    )
    if not loop_refinement:
        a = automodel(
            env,
            # alignment filename
            alnfile=alignment_file,
            # codes of the templates
            knowns=(alignment[0].id, ),
            # code of the target
            sequence=alignment[1].id,
            # wich method for validation should be calculated
            assess_methods=(assess.DOPE, assess.normalized_dope)
        )
    else:
        a = dope_loopmodel(
            env,
            # alignment filename
            alnfile=alignment_file,
            # codes of the templates
            knowns=(alignment[0].id, ),
            # code of the target
            sequence=alignment[1].id,
            # wich method for validation should be calculated
            assess_methods=(assess.DOPE, assess.normalized_dope),
            loop_assess_methods=(assess.DOPE, assess.normalized_dope)
        )
        # index of the first loop model
        a.loop.starting_model = 1  # self.loopStart
        # index of the last loop model
        a.loop.ending_model = 4  # self.loopEnd
        # loop refinement method; this yields
        a.loop.md_level = refine.slow

    a.starting_model = 1  # self.start  # index of the first model
    a.ending_model = elaspic.CONFIGS.get('modeller_runs', 1)  # index of the last model

    # Very thorough VTFM optimization:
    a.library_schedule = autosched.slow
    a.max_var_iterations = 300

    # Thorough MD optimization:
    # a.md_level = refine.slow
    a.md_level = None

    # Repeat the whole cycle 2 times and do not stop unless obj.func. > 1E6
    # a.repeat_optimization = 2

    a.max_molpdf = 2e5

    with py_tools.log_print_statements(logger):
        a.make()  # do the actual homology modeling

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
            logger.debug(
                'Success! model_filename: {}, model_dope_score: {}'
                .format(model_filename, model_dope_score))
            result.append((model_filename, model_dope_score))
        else:
            failure = a.outputs[i]['failure']
            logger.debug('Failure! {}'.format(failure))
            failures.append(a.outputs[i]['failure'])

    # Add the loop refinement output
    if loop_refinement:
        logger.debug('Modeller loop outputs:')
        for i in range(len(a.loop.outputs)):
            if not a.loop.outputs[i]['failure']:
                model_filename = a.loop.outputs[i]['name']
                model_dope_score = a.loop.outputs[i]['Normalized DOPE score']
                logger.debug(
                    'Success! model_filename: {}, model_dope_score: {}'
                    .format(model_filename, model_dope_score))
                result.append((model_filename, model_dope_score))
                loop = True
            else:
                failure = a.loop.outputs[i]['failure']
                logger.debug('Failure! {}'.format(failure))
                failures.append(failure)

    # Return the successfully calculated models and a loop flag indicating
    # whether the returned models are loop refined or not
    return result, loop, failures


def _write_alignment(alignment, tempdir):
    assert len(alignment) == 2
    alignment_file = op.join(tempdir, alignment[0].id + '-' + alignment[1].id + '.pir')
    with open(alignment_file, 'wt') as ofh:
        Bio.SeqIO.write(alignment, ofh, 'pir')
    return alignment_file
