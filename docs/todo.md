# TODO

## Table of contents

- [Thesis](#thesis)
- [mmCIF](#mmcif)
- [Multiple mutation support](#multiple-mutation-support)
- [Other ideas](#other-ideas)


## Thesis

- Get all interface mutations.
- See how well we can predict if those interface mutations are going to be deleterious.
- Should perform better than conservation...?
- *Edgetic* vs *lose all* vs *lose none*.


## mmCIF

The *PDBxReader.py* files from OpenMM, ParmEd, and hhblits are about the same. 
If the writers and parsers are also the same, should use ParmEd because it can through PyPy. 
PDBFixer uses the version from OpenMM. 

- Need to see what kind of conversion happens when you read and write the same file.
- Align sequence to structure.
- Rename chains and residues in a consistent way?
- See what kind of changes happen with PDBfixer.
- No point in cleaning all of PDB, but need a reference implementation of the PDBx reader.


## Multiple mutation support

We can use a simple heuristic to predict the ΔΔG value for multiple mutations:

  - Introduce each point mutation one at a time.
  - Keep the point mutation with the lowest ΔΔG value (i.e.  keep the most stabilising mutation).
  - Use the mutated protein as a new reference, and repeat (recursive :)).

Justifications:

  - We select the most stabilizing mutation in each run in order to emulate *protein design*.
  - We could even introduce an MD run between each mutation in order to emulate *flexible backbone design*.
  - Both Protherm and  have ΔΔG values for proteins with multiple mutations.
  - Skempi has **730 ΔΔG values** (~1/3) for mutations affecting multiple AA.
  - Could also use Phage display data.


## Other ideas

Use multitask learning to combine different datasets to make predictions.
