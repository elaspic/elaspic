import logging
import os.path as op

import numpy as np
import pytest
import yaml

from elaspic.call_foldx import FoldX, read_analyse_complex, read_build_model, read_stability
from elaspic.structure_tools import download_pdb_file

logger = logging.getLogger(__name__)

with open(op.join(op.splitext(__file__)[0], 'test_call_foldx.yaml')) as ifh:
    FOLDX_TEST_DATA = yaml.load(ifh)['foldx_test_data']


def test_read_build_model():
    pdb_id = '3zml'
    foldx_mutation = 'QA93A'
    output_file = op.join(
        op.splitext(__file__)[0], '{}-{}'.format(pdb_id, foldx_mutation),
        'Raw_{}-foldx.fxout'.format(pdb_id))
    stability_values_wt, stability_values_mut = read_build_model(
        output_file, 'WT_3zml-foldx_1.pdb', '3zml-foldx_1.pdb')
    assert np.allclose(stability_values_wt, [
        -151.067, -359.241, -147.835, -560.376, -20.8506, 732.907, -749.99, 19.7999, 291.225,
        653.47, 0, 0, 2.24811, 8.88635, 387.994, -22.3372, 0, 0, -1.73284, 0, 2.75874, 0
    ])
    assert np.allclose(stability_values_mut, [
        -150.757, -358.588, -147.255, -559.122, -20.8506, 731.301, -748.723, 19.7991, 289.654,
        653.214, 0, 0, 2.24811, 8.8786, 387.858, -22.3372, 0, 0, -1.73284, 0, 2.75874, 0
    ])


def test_read_stability():
    pdb_id = '3zml'
    foldx_mutation = 'QA93A'
    output_file = op.join(
        op.splitext(__file__)[0], '{}-{}'.format(pdb_id, foldx_mutation),
        '{}-foldx-{}-wt_0_ST.fxout'.format(pdb_id, foldx_mutation))
    stability_values_wt = read_stability(output_file)
    assert np.allclose(stability_values_wt, [
        -151.067, -359.241, -147.835, -560.376, -20.8506, 732.907, -749.99, 19.7999, 291.225,
        653.47, 0, 0, 2.24811, 8.88635, 387.994, -22.3372, 0, 0, -1.73284, 0, 2.75874, 0, 437
    ])


def test_read_analyse_model():
    pdb_id = '3zml'
    foldx_mutation = 'QA93A'
    output_file = op.join(
        op.splitext(__file__)[0], '{}-{}'.format(pdb_id, foldx_mutation),
        'Interaction_{}-foldx-{}-wt_AC.fxout'.format(pdb_id, foldx_mutation))
    stability_values_wt = read_analyse_complex(output_file)
    assert np.allclose(stability_values_wt, [
        12.0101, 12.9476, -34.0241, -3.99629, -18.0443, -27.6111, -6.64662, 36.8038, -36.4064,
        3.66039, 16.691, 1.18677, 0, 0, 0, 0.0681398, 6.23016, 0.978522, 0, 0, -1.73284, 0,
        1.02483, 2.384, 437, 85, 1, 1, 0
    ])


@pytest.mark.parametrize("test_data", FOLDX_TEST_DATA)
def test_foldx(test_data, tmpdir):
    tmp_dir = tmpdir.mkdtemp()
    pdb_file = download_pdb_file(test_data['pdb_id'], tmp_dir)
    pdb_chains = (test_data['pdb_chain'], test_data['pdb_partner_chain'])
    foldx_mutation = (
        test_data['pdb_mutation'][0] + test_data['pdb_chain'] + test_data['pdb_mutation'][1:])
    foldx = FoldX(tmp_dir)
    structure_file_wt, structure_file_mut = _test_build_model(pdb_file, foldx_mutation, foldx,
                                                              test_data['stability_results_'])
    _test_stability(structure_file_wt, structure_file_mut, foldx, test_data['stability_results_'])
    _test_analyse_complex(structure_file_wt, structure_file_mut, pdb_chains, foldx,
                          test_data['analyze_complex_results_'])


def _test_build_model(pdb_file, foldx_mutation, foldx, stability_results_):
    logger.info('Test BuildModel')
    structure_file_wt, structure_file_mut, stability_values_wt, stability_values_mut = (
        foldx.build_model(pdb_file, foldx_mutation))
    stability_results = [mut - wt for wt, mut in zip(stability_values_wt, stability_values_mut)]
    logger.debug("stability_values_wt: %s", stability_values_wt)
    logger.debug("stability_values_mut: %s", stability_values_mut)
    logger.debug("stability_results: %s", stability_results)
    logger.debug("stability_results_: %s", stability_results_)
    stability_results = [round(f, 3) for f in stability_results]
    assert stability_results == stability_results_
    return structure_file_wt, structure_file_mut


def _test_stability(structure_file_wt, structure_file_mut, foldx, stability_results_):
    logger.info('Test Stability')
    stability_values_wt = foldx.stability(structure_file_wt)
    stability_values_mut = foldx.stability(structure_file_mut)
    stability_results = [mut - wt for wt, mut in zip(stability_values_wt, stability_values_mut)]
    logger.debug("stability_values_wt: %s", stability_values_wt)
    logger.debug("stability_values_mut: %s", stability_values_mut)
    logger.debug("stability_results: %s", stability_results)
    logger.debug("stability_results_: %s", stability_results_)
    stability_results = [round(f, 3) for f in stability_results]
    assert stability_results == stability_results_


def _test_analyse_complex(structure_file_wt, structure_file_mut, pdb_chains, foldx,
                          analyze_complex_results_):
    logger.info('Test AnalyseComplex')
    analyze_complex_values_wt = foldx.analyse_complex(structure_file_wt, pdb_chains)
    analyze_complex_values_mut = foldx.analyse_complex(structure_file_mut, pdb_chains)
    analyze_complex_results = [
        mut - wt for wt, mut in zip(analyze_complex_values_wt, analyze_complex_values_mut)
    ]
    logger.debug("analyze_complex_values_wt: %s", analyze_complex_values_wt)
    logger.debug("analyze_complex_values_mut: %s", analyze_complex_values_mut)
    logger.debug("analyze_complex_results: %s", analyze_complex_results)
    logger.debug("analyze_complex_results: %s", analyze_complex_results_)
    analyze_complex_results = [round(f, 3) for f in analyze_complex_results]
    assert analyze_complex_results == analyze_complex_results_
