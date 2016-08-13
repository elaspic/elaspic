# ELASPIC TODO

## Table of contents

- [Thesis](#thesis)
- [mmCIF](#mmcif)
- [Multiple mutation support](#multiple-mutation-support)
- [Improved sequence alignments](#improved-sequence-alignments)
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


## Improved sequence alignments

- Use only SwissProt canonical sequences.
- Map mutations to SwissProt using gene-level identifiers and sequence alignments.

Advantages of using `hhpred`:
: - Conservation score from MSAs should we useful for predictions.
  - Protein-drug homology models would be useful for tieing together ELASPIC and drug synergy projects.
  - Interlinking human proteome and PDB can be useful for other projects.


Use multitask learning to combine different datasets to make predictions.

1. Extract.

2. Create an `hhblits` PDB database with *all* PDB sequences.
  - Extract PDB sequences from UniParc.
  - Query using Gene3D MSAs?

3. Use `hhblits` to create global alignments against the PDB for all Gene3D domains.
  - Save alignments in the database using the `a2b` and `b2a` format.

4. `hhblits` uses UniProt20 to create an MSA for the query sequence.

[hh-suite](https://github.com/soedinglab/hh-suite)


## Other ideas

```
- elaspic
|-- core
|-- dna
|-- peptide
|-- protein
|-- small_molecule
```

License: LGPL


[hh-suite](https://github.com/soedinglab/hh-suite)


[hhpred](http://doi.org/10.1038/nmeth.1818)
: - Creates a MSA for the query protein using the UniProt20 database.
  - Converts the MSA into an HMM.
  - Performs HMM-HMM alignments against the PDB.


[hhblits documentation](https://toolkit.tuebingen.mpg.de/hhblits/help_params)

[hhpred documentation](https://toolkit.tuebingen.mpg.de/hhpred/help_params)

### core

### dna

### peptide

### protein

### small_molecule

- [Platinum: a database of experimentally measured effects of mutations on structurally defined protein–ligand complexes](http://doi.org/10.1093/nar/gku966)
- [Mycobacterium tuberculosis whole genome sequencing and protein structure modelling provides insights into anti-tuberculosis drug resistance](http://doi.org/10.1186/s12916-016-0575-9)


### Additional features

[sdm](http://mordred.bioc.cam.ac.uk/~sdm/sdm_theory.php)

[Cutoff Scanning Matrix (CSM): structural classification and function prediction by protein inter-residue distance patterns](http://doi.org/10.1186/1471-2164-12-S4-S12)
: - mCSM
  - mCSM-AB
  - DUET
  - pkCSM
  - aCSM

## Post-translational modifications

<p align="center">
<img src="todo/current_setup_crop.png" width="600px" >
</p>
<br>
<p align="center">
<img src="todo/proposed_setup_crop.png" width="800px" >
</p>
