from elaspic.tools._abc import ToolError


# Generate pose from primary structure
from rosetta import init, pose_from_sequence

init()

seq = 12 * 'A'
pose = pose_from_sequence(seq)


# Mutate a residue
from toolbox import mutate_residue

mutant = mutate_residue(pose, 5, 'K')


# Load Pose into ParmEd
from parmed import load_rosetta

struct = load_rosetta(mutant)


# Start an OpenMM simulation
from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *

ff = ForceField('amber99sbildn.xml', 'tip3p.xml')

# Solvate Structure
mod = Modeller(struct.topology, struct.positions)
mod.addSolvent(ff, model='tip3p', boxSize=Vec3(4, 4, 4)*nanometers,
               positiveIon='Na+', negativeIon='Cl-',
               ionicStrength=0.1*molar)

# Create OpenMM System
system = ff.createSystem(mod.topology, nonbondedMethod=PME,
                         nonbondedCutoff=1*nanometers,
                         constraints=HBonds)
integrator = LangevinIntegrator(300*kelvin, 1/picoseconds,
                                2*femtoseconds)

system.addForce(MonteCarloBarostat(1*bar, 300*kelvin))

simulation = Simulation(mod.topology, system, integrator)
simulation.context.setPositions(mod.positions)

# Minimize System
simulation.minimizeEnergy(maxIterations=1000)

# Equilibrate System
simulation.reporters.append(
    PDBReporter('dodecaalanine.solv.pdb', 50000))

simulation.context.setVelocitiesToTemperature(300)

simulation.step(50000)




class RosettaError(ToolError):
    pass


def mutate():
    ...
