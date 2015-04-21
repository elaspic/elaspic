.. _get_started:

Getting started 
================

.. only:: html

   .. sidebar:: Flowchart describing every step of the ELASPIC pipeline.
     
     .. image:: ./elaspic_flowchart.*
        :target: ./../_downloads/elaspic_flowchart.png
      
.. only:: latex

   .. image:: ./elaspic_flowchart.pdf
      :scale: 40
      :align: center


The general overview of the ELASPIC pipleine is presented in the :download:`figure <./elaspic_flowchart.png>` to the right. A user runs the ELASPIC pipeline specifying the Uniprot id of the protein being mutated, and one or more mutations affecting that protein. At each decision node, the pipeline queries the database to check whether or not the required information has been previously calculated. If the required data has not been calculated, the pipeline calculates it on the fly and stores the results in the database for later retrieval. The pipeline proceeds until homology models of all domains in the protein, and all domain-domain interactions involving the protein, have been calculated, and the :math:`\Delta \Delta G` has been predicted for every specified mutation. 


..
   :width: 400px
   :align: center
   :subtitle: 


.. toctree::
   :maxdepth: 2
   :numbered:
   :hidden:
   
   elaspic_cli
   elaspic_database_cli

