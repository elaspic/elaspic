# Alignments

## Improved sequence alignments

- Use only SwissProt canonical sequences.
- Map mutations to SwissProt using gene-level identifiers and sequence alignments.

Advantages of using `hhpred`:
: - Conservation score from MSAs should we useful for predictions.
  - Protein-drug homology models would be useful for tieing together ELASPIC and drug synergy projects.
  - Interlinking human proteome and PDB can be useful for other projects.


1. Extract.

2. Create an `hhblits` PDB database with *all* PDB sequences.
  - Extract PDB sequences from UniParc.
  - Query using Gene3D MSAs?

3. Use `hhblits` to create global alignments against the PDB for all Gene3D domains.
  - Save alignments in the database using the `a2b` and `b2a` format.

4. `hhblits` uses UniProt20 to create an MSA for the query sequence.

[hh-suite](https://github.com/soedinglab/hh-suite)

[hh-suite](https://github.com/soedinglab/hh-suite)

[hhpred](http://doi.org/10.1038/nmeth.1818)
: - Creates a MSA for the query protein using the UniProt20 database.
  - Converts the MSA into an HMM.
  - Performs HMM-HMM alignments against the PDB.

[hhblits documentation](https://toolkit.tuebingen.mpg.de/hhblits/help_params)

[hhpred documentation](https://toolkit.tuebingen.mpg.de/hhpred/help_params)
