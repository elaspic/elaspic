Database schema
===============


.. image:: ./figures/elaspic_schema.*
   
.. only:: html

  Download: 
  :download:`pdf <./figures/elaspic_schema.pdf>`
  :download:`png <./figures/elaspic_schema.png>`
  :download:`mwb <./figures/elaspic.mwb>`


domain
  Description of all domains in the PDB, extracted from `scowlp`.

domain_contact
  Description of all domain_domain contacts in the PDB, extracted from `scowlp`.

uniprot_domain
  All Pfam domains found in a 2014 version of Uniprot / Trembl. This table was obtained by 
  reformatting precalculated Pfam domains from the SIMAP website, and mapping them to Uniprot IDs.

uniprot_domain_template
  Structural templates for domains in the `uniprot_domain` table. Lists PDB crystal structures
  that will be used for making homology models. The column `domain_def` is the PDBfam domain 
  definition, which is a combination of PFam and Gene3D domains.

uniprot_domain_model
  Homology models for templates in the `uniprot_domain_template` table. The column `model_domain_def`
  describes the region of the domain that has structural coverage.

uniprot_domain_mutation
  Characterization of mutations introduced into structures in the `uniprot_domain_model` table.

uniprot_domain_pair
  Pairwise combinations of domains for proteins that are known to interact, according to 
  Hippie, IRefIndex, and the PDB. 

uniprot_domain_pair_template
  Structural templates for pairs of domains in the `uniprot_domain_pair` table.

uniprot_domain_pair_model
  Structural models of interactions between pairs of domains in the `uniprot_domain_pair_model`
  table.

uniprot_domain_pair_mutation
  Characterization of interface mutations introduced into structures in the `uniprot_domain_pair_model`
  table.




