"""Homology modeling by the automodel class.


Can use either dunamic_sphere or dynamic_lennard contraints
https://salilab.org/modeller/9.17/manual/node108.html


selection.randomize_xyz

"""
import os
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

import modeller
from modeller.optimizers import molecular_dynamics, conjugate_gradients

from kmtools import py_tools, structure_tools

logger = logging.getLogger(__name__)
ToolError.register(ModellerError)


class SaliLabModeller(Modeller):

    _result_slots = Modeller._result_slots + ['norm_dope_score']

    def _build(self):
        self.result['homology_structure_file'], self.result['norm_dope_score'] = (
            _build_modeller(self.alignment, self.tempdir)
        )


class SaliLabMutator(Mutator):

    _result_slots = Mutator._result_slots + ['norm_dope_score']

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


def _optimize(atmsel, sched):
    """Run conjugate gradient descent to optimize the model."""
    for step in sched:
        step.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)
    # md
    _refine(atmsel)
    cg = conjugate_gradients()
    cg.optimize(atmsel, max_iterations=200, min_atom_shift=0.001)


def _refine(atmsel):
    """Run molecular dynamics to refine the model."""
    # at T=1000, max_atom_shift for 4fs is cca 0.15 A.
    md = molecular_dynamics(cap_atom_shift=0.39, md_time_step=4.0,
                            md_return='FINAL')
    init_vel = True
    for (its, equil, temps) in ((200, 20, (150.0, 250.0, 400.0, 700.0, 1000.0)),
                                (200, 600, (1000.0, 800.0, 600.0, 500.0, 400.0, 300.0))):
        for temp in temps:
            md.optimize(atmsel, init_velocities=init_vel, temperature=temp,
                        max_iterations=its, equilibrate=equil)
            init_vel = False


def _make_restraints(mdl1, aln):
    """use homologs and dihedral library for dihedral angle restraints
    """
    rsr = mdl1.restraints
    rsr.clear()
    s = modeller.selection(mdl1)
    for typ in ('stereo', 'phi-psi_binormal'):
        rsr.make(s, restraint_type=typ, aln=aln, spline_on_site=True)
    for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
        rsr.make(s, restraint_type=typ + '_dihedral', spline_range=4.0,
                 spline_dx=0.3, spline_min_points=5, aln=aln,
                 spline_on_site=True)


def _modeller_mutate(structure_file, mutation):
    """
    Source: https://salilab.org/modeller/wiki/Mutate%20model

    mutate_model.py

        Usage:   python mutate_model.py modelname respos resname chain > logfile

        Example: python mutate_model.py 1t29 1699 LEU A > 1t29.log

    Creates a single in silico point mutation to sidechain type and at residue position
    input by the user, in the structure whose file is modelname.pdb
    The conformation of the mutant sidechain is optimized by conjugate gradient and
    refined using some MD.

    Note: if the model has no chain identifier, specify "" for the chain argument.
    """

    # first argument
    # modelname, respos, restyp, chain, = sys.argv[1:]
    CWD = os.getcwd()

    tempdir = op.dirname(structure_file)
    os.chdir(tempdir)

    structure_file = op.basename(structure_file)
    pdb_id = op.splitext(structure_file)[0]
    modelname = pdb_id

    chain_id, chain_mutation = mutation.split('-')
    chain = chain_id

    respos = int(mutation[1:-1])
    restyp = structure_tools.A_DICT[mutation[-1]]

    log.verbose()

    # Set a different value for rand_seed to get a different final model
    env = environ(rand_seed=42)

    env.io.hetatm = True
    # soft sphere potential
    env.edat.dynamic_sphere = False
    # lennard-jones potential (more accurate)
    env.edat.dynamic_lennard = True
    env.edat.contact_shell = 4.0
    env.edat.update_dynamic = 0.39

    # Read customized topology file with phosphoserines (or standard one)
    env.libs.topology.read(file='$(LIB)/top_heav.lib')

    # Read customized CHARMM parameter library with phosphoserines (or standard one)
    env.libs.parameters.read(file='$(LIB)/par.lib')

    # Read the original PDB file and copy its sequence to the alignment array:
    mdl1 = modeller.model(env, file=modelname)
    ali = modeller.alignment(env)
    ali.append_model(mdl1, atom_files=modelname, align_codes=modelname)

    # set up the mutate residue selection segment
    s = modeller.selection(mdl1.chains[chain].residues[respos])

    # perform the mutate residue operation
    s.mutate(residue_type=restyp)
    # get two copies of the sequence.  A modeller trick to get things set up
    ali.append_model(mdl1, align_codes=modelname)

    # Generate molecular topology for mutant
    mdl1.clear_topology()
    mdl1.generate_topology(ali[-1])

    # Transfer all the coordinates you can from the template native structure
    # to the mutant (this works even if the order of atoms in the native PDB
    # file is not standard):
    # here we are generating the model by reading the template coordinates
    mdl1.transfer_xyz(ali)

    # Build the remaining unknown coordinates
    mdl1.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')

    # yes model2 is the same file as model1.  It's a modeller trick.
    mdl2 = modeller.model(env, file=modelname)

    # required to do a transfer_res_numb
    # ali.append_model(mdl2, atom_files=modelname, align_codes=modelname)
    # transfers from "model 2" to "model 1"
    mdl1.res_num_from(mdl2, ali)

    # It is usually necessary to write the mutated sequence out and read it in
    # before proceeding, because not all sequence related information about MODEL
    # is changed by this command (e.g., internal coordinates, charges, and atom
    # types and radii are not updated).
    mdl1.write(file=modelname + restyp + respos + '.tmp')
    mdl1.read(file=modelname + restyp + respos + '.tmp')

    # set up restraints before computing energy
    # we do this a second time because the model has been written out and read in,
    # clearing the previously set restraints
    _make_restraints(mdl1, ali)

    # a non-bonded pair has to have at least as many selected atoms
    mdl1.env.edat.nonbonded_sel_atoms = 1

    sched = autosched.loop.make_for_model(mdl1)

    # only optimize the selected residue (in first pass, just atoms in selected
    # residue, in second pass, include nonbonded neighboring atoms)
    # set up the mutate residue selection segment
    s = modeller.selection(mdl1.chains[chain].residues[respos])

    mdl1.restraints.unpick_all()
    mdl1.restraints.pick(s)

    s.energy()

    s.randomize_xyz(deviation=4.0)

    mdl1.env.edat.nonbonded_sel_atoms = 2
    _optimize(s, sched)

    # feels environment (energy computed on pairs that have at least one member
    # in the selected)
    mdl1.env.edat.nonbonded_sel_atoms = 1
    _optimize(s, sched)

    s.energy()

    # give a proper name
    mdl1.write(file=modelname + restyp + respos + '.pdb')

    # delete the temporary file
    os.remove(modelname + restyp + respos + '.tmp')

    os.chdir(CWD)


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
                result, loop, failures = _modeller_model(alignment, tempdir, True)
            except Exception as e:
                logger.error('Loop refinement failed with an error: {}'.format(e))
                try:
                    result, loop, failures = _modeller_model(alignment, tempdir, False)
                except ModellerError as e:
                    raise elaspic.exc.ModellerError(e)
                except Exception as e:
                    raise e
            if not result:
                raise elaspic.exc.ModellerError(failures[-1])
            results.extend(result)
        results.sort(key=lambda x: x[1])
        return results[0][0], results[0][1]


def _modeller_model(alignment, tempdir, loop_refinement):
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
