import logging
import os
import os.path as op
import tempfile

import elaspic.structure_analysis

logger = logging.getLogger(__name__)


class TestAnalyseStructure:
    @classmethod
    def setup_class(cls):
        pdb_file = op.abspath(
            op.join(op.splitext(__file__)[0], "4CPA.ENTI_1_PDB4CPA.ENTB_2-4CPAIBIB.pdb")
        )
        tempdir = op.join(op.abspath(op.dirname(__file__)), "tmp")
        os.makedirs(tempdir, exist_ok=True)
        tempfile.tempdir = tempdir
        working_dir = tempfile.mkdtemp()
        cls.analyse_structure = elaspic.structure_analysis.AnalyzeStructure(
            pdb_file=pdb_file,
            working_dir=working_dir,
            vdw_distance=5.0,
            min_contact_distance=6.0,
        )

    @classmethod
    def teardown_class(cls):
        pass

    def test_get_seasa(self):
        (
            seasa_by_chain,
            seasa_by_chain_separately,
            seasa_by_residue,
            seasa_by_residue_separately,
        ) = self.analyse_structure.get_seasa()
