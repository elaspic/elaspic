.. _introduction:

Introduction
============

.. figure:: ./introduction/elaspic_flowchart.png
   :target: _downloads/elaspic_flowchart.pdf
   :align: center
   :figclass: align-center

   Flowchart describing the ELASPIC pipeline.

The general overview of ELASPIC is presented in the figure above.

ELASPIC can be run using two different pipelines: the :ref:`local pipeline` and the :ref:`database pipeline`.


.. _`local pipeline`:

Local pipeline
--------------

.. only:: html

   .. sidebar:: Local pipeline flowchart.

      .. image:: ./introduction/elaspic_flowchart_local.*
         :target: _downloads/elaspic_flowchart_local.png
         :width: 100 px
         :align: center

.. only:: latex

   .. image:: ./introduction/elaspic_flowchart_local.pdf
      :scale: 40
      :align: center

The local pipeline works without downloading and installing a local copy of the ELASPIC and PDB databases, but requires a PDB structure or template to be provided for every protein. This pipeline still requires a local copy of the BLAST database.

The general overview of the database pipleine is presented in the :download:`figure <./introduction/elaspic_flowchart_local.png>` to the right. A user runs the ELASPIC pipeline specifying the Uniprot id of the protein being mutated, and one or more mutations affecting that protein. At each decision node, the pipeline queries the database to check whether or not the required information has been previously calculated. If the required data has not been calculated, the pipeline calculates it on the fly and stores the results in the database for later retrieval. The pipeline proceeds until homology models of all domains in the protein, and all domain-domain interactions involving the protein, have been calculated, and the :math:`\Delta \Delta G` has been predicted for every specified mutation.


.. _`database pipeline`:

Database pipeline
-----------------

.. only:: html

   .. sidebar:: Database pipeline flowchart.

      .. image:: introduction/elaspic_flowchart_database.*
         :target: _downloads/elaspic_flowchart_database.png
         :width: 100 px
         :align: center

.. only:: latex

   .. image:: introduction/elaspic_flowchart_database.pdf
      :scale: 40
      :align: center

The database pipeline allows mutations to be performed on a proteome-wide scale, without having to specify a structural template for each protein. This pipeline requires a local copy of Profs domain definitions and templates, as well as a local copy of the BLAST and PDB databases.

The general overview of the database pipleine is presented in the :download:`figure <./introduction/elaspic_flowchart_database.png>` to the right. A user runs the ELASPIC pipeline specifying the Uniprot id of the protein being mutated, and one or more mutations affecting that protein. At each decision node, the pipeline queries the database to check whether or not the required information has been previously calculated. If the required data has not been calculated, the pipeline calculates it on the fly and stores the results in the database for later retrieval. The pipeline proceeds until homology models of all domains in the protein, and all domain-domain interactions involving the protein, have been calculated, and the :math:`\Delta \Delta G` has been predicted for every specified mutation.
