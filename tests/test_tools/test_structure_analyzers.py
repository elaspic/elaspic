from pprint import pprint

import pytest

import elaspic.tools
from kmtools import py_tools, structure_tools

logger = py_tools.get_logger(__name__)

ANALYZERS = [
    elaspic.tools.Stride,
    elaspic.tools.MSMS,
    elaspic.tools.PhysicoChemical,
]
PDB_IDS = [
    ('2wru', 'A-B', 'G1A:I2A:V3A:E4A:Q5A:C6A:C7A:T8A:S9A:I10A:C11A:S12A:N21A'),
    ('4dkl', 'A', 'M65A:M130A:Y252A:R263A:P295A:I302A:I352A'),
]


# @pytest.mark.parametrize("pdb_id, chain_id, mutations", PDB_IDS)
# @pytest.mark.parametrize("Analyzer", Analyzers)
# def test_1(Analyzer, pdb_id, chain_id, mutations):
#     # Prepare structure
#     s = structure_tools.fetch_structure(pdb_id, pdb_type='cif', biounit=False)
#     logger.debug("{}, {}, {}", s, len(s[0]), len(s[0]['A']))
#     # Prepare analyzer
#     analyzer = Analyzer(s)
#     analyzer.build()
#     assert analyzer.done
#     for mutation in mutations.split(':'):
#         logger.debug('{}-{}', chain_id, mutations)
#         result = analyzer.analyze(chain_id, int(mutation[1:-1]), mutation[0])
#         logger.debug("{}\n", result)
#
#
# @pytest.mark.parametrize("pdb_id, chain_id, mutations", PDB_IDS)
# @pytest.mark.parametrize("Analyzer", Analyzers)
# def test_2(Analyzer, pdb_id, chain_id, mutations):
#     # Prepare structure
#     s = structure_tools.fetch_structure(pdb_id, pdb_type='cif', biounit=False)
#     structure_tools.process_structure(s)
#     structure_tools.extract(s, [(0, chain_id, range(len(s[0][chain_id])))])
#     logger.debug("{}, {}, {}", s, len(s[0]), len(s[0]['A']))
#     # Prepare analyzer
#     analyzer = Analyzer(s)
#     analyzer.build()
#     assert analyzer.done
#     for mutation in mutations.split(':'):
#         _, new_chain_id, new_residue_id = s.xtra['residue_mapping_fw'][
#             (0, chain_id, (' ', int(mutation[1:-1]), ' '))]
#         logger.debug('{}-{}: {}-{}', chain_id, mutation, new_chain_id, new_residue_id)
#         result = analyzer.analyze(new_chain_id, new_residue_id, mutation[0])
#         logger.debug("{}\n", result)
#

@pytest.mark.parametrize("pdb_id, chain_ids, mutations", PDB_IDS)
def test_3(pdb_id, chain_ids, mutations):
    # Prepare structure
    chain_id, _, chain_id_2 = chain_ids.partition('-')
    structure = structure_tools.fetch_structure(pdb_id, pdb_type='cif', biounit=False)
    logger.debug("{}, {}, {}", structure, len(structure[0]), len(structure[0]['A']))
    # Mutate
    mutator = elaspic.tools.FoldXMutator(structure)
    mutator.build()
    for chain_mutation in mutations.split(':'):
        result = {}
        mutation = '{}-{}'.format(chain_id, chain_mutation)
        logger.debug('Mutation: {}', mutation)
        structure_wt, structure_mut, _result = mutator.mutate(mutation)
        result.update(_result)
        for Analyzer in ANALYZERS:
            logger.debug('Analyzer: {}'.format(Analyzer))
            # WT
            analyzer_wt = Analyzer(structure_wt)
            analyzer_wt.build()
            assert analyzer_wt.done
            result_wt = analyzer_wt.analyze(mutation)
            result[Analyzer.__class__.__name__ + '_wt'] = result_wt
            # MUT
            analyzer_mut = Analyzer(structure_mut)
            analyzer_mut.build()
            assert analyzer_mut.done
            result_mut = analyzer_mut.analyze(mutation)
            result[Analyzer.__class__.__name__ + '_mut'] = result_mut
        pprint(result)


# @pytest.mark.parametrize("pdb_id, chain_ids, mutations", PDB_IDS)
# def test_4(pdb_id, chain_ids, mutations):
#     # Prepare structure
#     chain_id, _, chain_id_2 = chain_ids.partition('-')
#     structure = structure_tools.fetch_structure(pdb_id, pdb_type='cif', biounit=False)
#     structure_tools.process_structure(structure)
#     if chain_id_2:
#         domain_defs = [
#             (0, chain_id, range(len(structure[0][chain_id]))),
#             (0, chain_id_2, range(len(structure[0][chain_id_2]))),
#         ]
#     else:
#         domain_defs = [
#             (0, chain_id, range(len(structure[0][chain_id]))),
#         ]
#     structure_tools.extract(structure, domain_defs)
#     logger.debug("{}, {}, {}", structure, len(structure[0]), len(structure[0]['A']))
#     # Mutate
#     mutator = elaspic.tools.FoldXMutator(structure)
#     mutator.build()
#     for chain_mutation in mutations.split(':'):
#         result = {}
#         _, new_chain_id, new_residue_id = structure.xtra['residue_mapping_fw'][
#             (0, chain_id, (' ', int(chain_mutation[1:-1]), ' '))]
#         mutation = '{}-{}'.format(
#             new_chain_id, chain_mutation[0] + str(new_residue_id[1]) + chain_mutation[-1])
#         logger.debug('Mutation: {}', mutation)
#         structure_wt, structure_mut, _result = mutator.mutate(mutation)
#         result.update(_result)
#         for Analyzer in ANALYZERS:
#             logger.debug('Analyzer: {}'.format(Analyzer))
#             # WT
#             analyzer_wt = Analyzer(structure_wt)
#             analyzer_wt.build()
#             assert analyzer_wt.done
#             result_wt = analyzer_wt.analyze(mutation)
#             result[Analyzer.__class__.__name__ + '_wt'] = result_wt
#             # MUT
#             analyzer_mut = Analyzer(structure_mut)
#             analyzer_mut.build()
#             assert analyzer_mut.done
#             result_mut = analyzer_mut.analyze(mutation)
#             result[Analyzer.__class__.__name__ + '_mut'] = result_mut
#         pprint(result)
