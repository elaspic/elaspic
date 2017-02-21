import pytest

from kmtools import py_tools, structure_tools

import elaspic.tools

logger = py_tools.get_logger(__name__)

Analyzers = [elaspic.tools.MSMS, elaspic.tools.PhysicoChemical]
PDB_IDS = [('4dkl', 'A', 'M65A:M130A:Y252A:R263A:P295A:I302A:I352A'), ]


@pytest.mark.parametrize("pdb_id, chain_id, mutations", PDB_IDS)
@pytest.mark.parametrize("Analyzer", Analyzers)
def test_1(Analyzer, pdb_id, chain_id, mutations):
    # Prepare structure
    s = structure_tools.fetch_structure(pdb_id, pdb_type='cif', biounit=False)
    logger.debug("{}, {}, {}", s, len(s[0]), len(s[0]['A']))
    # Prepare analyzer
    analyzer = Analyzer(s)
    analyzer.build()
    assert analyzer.done
    for mutation in mutations.split(':'):
        logger.debug('{}-{}', chain_id, mutations)
        result = analyzer.analyze(chain_id, int(mutation[1:-1]), mutation[0])
        logger.debug("{}\n", result)


@pytest.mark.parametrize("pdb_id, chain_id, mutations", PDB_IDS)
@pytest.mark.parametrize("Analyzer", Analyzers)
def test_2(Analyzer, pdb_id, chain_id, mutations):
    # Prepare structure
    s = structure_tools.fetch_structure(pdb_id, pdb_type='cif', biounit=False)
    structure_tools.process_structure(s)
    structure_tools.extract(s, [(0, chain_id, range(len(s[0][chain_id])))])
    logger.debug("{}, {}, {}", s, len(s[0]), len(s[0]['A']))
    # Prepare analyzer
    analyzer = Analyzer(s)
    analyzer.build()
    assert analyzer.done
    for mutation in mutations.split(':'):
        _, new_chain_id, new_residue_id = s.xtra['residue_mapping_fw'][
            (0, chain_id, (' ', int(mutation[1:-1]), ' '))]
        logger.debug('{}-{}: {}-{}', chain_id, mutation, new_chain_id, new_residue_id)
        result = analyzer.analyze(new_chain_id, new_residue_id, mutation[0])
        logger.debug("{}\n", result)
