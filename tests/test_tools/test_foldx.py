import os.path as op
import shutil
import tempfile
from pprint import pprint

import pytest

from elaspic.tools import foldx
from kmtools import py_tools, structure_tools

logger = py_tools.get_logger(__name__)

PDB_IDS = [
    ('2wru', 'A-B', 'G1A:I2A:V3A:E4A:Q5A:C6A:C7A:T8A:S9A:I10A:C11A:S12A:N21A'),
]


@pytest.mark.parametrize("pdb_id, chain_ids, mutations", PDB_IDS)
def test_1(pdb_id, chain_ids, mutations):
    tempdir = tempfile.mkdtemp()
    chain_id, _, chain_id_2 = chain_ids.partition('-')
    print(tempdir)

    structure = structure_tools.fetch_structure(pdb_id)
    structure_file = op.join(tempdir, structure.id + '.pdb')
    structure_tools.save_structure(structure, structure_file)
    print(structure_file)

    repaired_structure_file = foldx._foldx_repair_pdb(structure_file)

    for mutation in mutations.split(':'):
        mutation = '{}-{}'.format(chain_id, mutation)
        structure_wt, structure_mut, result = foldx._foldx_build_model(
            repaired_structure_file, mutation)
        assert result['foldx_wt'] != result['foldx_mut']
        if chain_id_2:
            result['foldx_interface_wt'] = foldx._foldx_analyse_complex(
                structure_wt, [chain_id, chain_id_2])
            result['foldx_interface_mut'] = foldx._foldx_analyse_complex(
                structure_mut, [chain_id, chain_id_2])
            assert result['foldx_interface_wt'] != result['foldx_interface_mut']
        pprint(result)

    shutil.rmtree(tempdir)


# @pytest.mark.parametrize("pdb_id, chain_id, mutations", PDB_IDS)
# def test_2(pdb_id, chain_id, mutations):
#     tempdir = tempfile.mkdtemp()
#     print(tempdir)
#
#     structure = structure_tools.fetch_structure(pdb_id)
#     structure_file = op.join(tempdir, structure.id + '.pdb')
#     structure_tools.save_structure(structure, structure_file)
#     print(structure_file)
#
#     repaired_structure_file = op.join(
#         op.dirname(op.abspath(__file__)),
#         'foldx',
#         '{}-foldx.pdb'.format(pdb_id))
#
#     for mutation in mutations.split(':'):
#         mutation = '{}-{}'.format(chain_id, mutation)
#         wt_structure_file, mut_structure_file = foldx._foldx_build_model(repaired_structure_file)
#         result_core_wt = foldx._foldx_stability(wt_structure_file)
#         result_core_mut = foldx._foldx_stability(mut_structure_file)
#         result_interface_wt = foldx._foldx_analyse_complex(wt_structure_file)
#         result_interface_mut = foldx._foldx_analyse_complex(mut_structure_file)
#         result = {
#             'foldx_core_wt': result_core_wt,
#             'foldx_core_mut': result_core_mut,
#             'foldx_interface_wt': result_interface_wt,
#             'foldx_interface_mut': result_interface_mut,
#         }
#         print(result)
#     shutil.rmtree(tempdir)
