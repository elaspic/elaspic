import pytest

from kmtools import py_tools, structure_tools

import elaspic.tools

logger = py_tools.get_logger(__name__)

PDB_IDS = [('4dkl', 'A', 'M130A.Y252A.R263A.P295A.I302A.I352A'), ]


@pytest.mark.parametrize("pdb_id, pdb_chain, pdb_mutations", PDB_IDS)
def test_msms_1(pdb_id, pdb_chain, pdb_mutations):
    s = structure_tools.fetch_structure(pdb_id, pdb_type='cif', biounit=False)
    msms = elaspic.tools.MSMS(s)
    msms.build()
    assert msms.done
    for pdb_mutation in pdb_mutations.split('.'):
        logger.debug('{}-{}', pdb_chain, pdb_mutation)
        result = msms.mutate(pdb_chain, pdb_mutation)
        logger.debug("{}\n", result)


@pytest.mark.parametrize("pdb_id, chain_id, mutations", PDB_IDS)
def test_msms_2(pdb_id, chain_id, mutations):
    s = structure_tools.fetch_structure(pdb_id, pdb_type='cif', biounit=False)
    structure_tools.process_structure(s)
    structure_tools.extract(s, [(0, chain_id, range(len(s[0][chain_id])))])
    logger.info("{}, {}, {}", s, len(s[0]), len(s[0]['A']))
    msms = elaspic.tools.MSMS(s)
    msms.build()
    assert msms.done
    for mutation in mutations.split('.'):
        key = (0, chain_id, (' ', int(mutation[1:-1]), ' '))
        new_model, new_chain_id, new_mutation_pos = s.xtra['residue_mapping_fw'][key]
        new_mutation = mutation[0] + str(new_mutation_pos[1]) + mutation[-1]
        logger.debug('{}-{}: {}-{}', chain_id, mutation, new_chain_id, new_mutation)
        result = msms.mutate(new_chain_id, new_mutation)
        logger.debug("{}\n", result)
