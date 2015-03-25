Database schema
===============

.. image:: ./figures/elaspic_schema.*

.. only:: html

  Download: 
  :download:`pdf <./figures/elaspic_schema.pdf>`
  :download:`png <./figures/elaspic_schema.png>`
  :download:`mwb <./figures/elaspic.mwb>`


----------------

The :term:`domain` and :term:`domain_contact` tables can be downloaded from the `structural_templates`
folder in the `ELASPIC downloads page`_. All other tables, for individual organisms and for all
Uniprot sequences, can be downloaded from the `elaspic` folder in the `ELASPIC downloads page`_.


.. glossary::

   domain
     Profs domain definitions for all proteins in the PDB. 

   domain_contact
     Interactions between domains in the PDB which passed quality control criteria.

   uniprot_sequence
     Protein sequences from the Uniprot KB, obtained by parsing `uniprot_sprot_fasta.gz`, `uniprot_trembl_fasta.gz`, and `homo_sapiens_variation.txt` files from the `Uniprot ftp site`_.

   uniprot_domain
     Pfam domain definitions for proteins in the :term:`uniprot_sequence` table. This table was obtained by downloading precalculated Pfam domain definitions from the `SIMAP`_ website and mapping those proteins to Uniprot IDs.

   uniprot_domain_template
     Structural templates for domains in the :term:`uniprot_domain` table. Lists PDB crystal structures 
     that will be used for making homology models. The column `domain_def` is the PDBfam domain 
     definition, which is a combination of PFam and Gene3D domains.

   uniprot_domain_model
     Homology models for templates in the :term:`uniprot_domain_template` table. The column `model_domain_def` describes the region of the domain that has structural coverage.

   uniprot_domain_mutation
     Characterization of mutations introduced into structures in the :term:`uniprot_domain_model` table.

   uniprot_domain_pair
     Potentially-interacting pairs of domains for proteins that are known to interact, according to 
     `Hippie`_, `IRefIndex`_, and `Rolland et al. 2014`_.

   uniprot_domain_pair_template
     Structural templates for pairs of domains in the :term:`uniprot_domain_pair` table.

   uniprot_domain_pair_model
     Structural models of interactions between pairs of domains in the :term:`uniprot_domain_pair_model`
     table.

   uniprot_domain_pair_mutation
     Characterization of interface mutations introduced into structures in the :term:`uniprot_domain_pair_model` table.


.. _SIMAP: http://liferay.csb.univie.ac.at/portal/web/simap
.. _Hippie: http://cbdm.mdc-berlin.de/tools/hippie/
.. _IRefIndex: http://irefindex.org
.. _Rolland et al. 2014: http://dx.doi.org/10.1016/j.cell.2014.10.050
.. _Profs: https://bitbucket.org/afgiraldofo/pdbfam
.. _ELASPIC downloads page: http://elaspic.kimlab.org/static/download/
.. _Uniprot ftp site: ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/

