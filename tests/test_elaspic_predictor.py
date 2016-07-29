import os.path as op
import pandas as pd
import elaspic.elaspic_predictor


_foldx_core_features = [
    # FoldX
    # (wildtype)
    'dg_wt',
    'backbone_hbond_wt', 'sidechain_hbond_wt', 'van_der_waals_wt',
    'electrostatics_wt', 'solvation_polar_wt', 'solvation_hydrophobic_wt',
    'van_der_waals_clashes_wt', 'entropy_sidechain_wt', 'entropy_mainchain_wt',
    'sloop_entropy_wt', 'mloop_entropy_wt', 'cis_bond_wt', 'torsional_clash_wt',
    'backbone_clash_wt', 'helix_dipole_wt', 'water_bridge_wt', 'disulfide_wt',
    'electrostatic_kon_wt', 'partial_covalent_bonds_wt', 'energy_ionisation_wt',
    'entropy_complex_wt',
    'number_of_residues',
    # (change)
    'dg_change',
    'backbone_hbond_change', 'sidechain_hbond_change', 'van_der_waals_change',
    'electrostatics_change', 'solvation_polar_change', 'solvation_hydrophobic_change',
    'van_der_waals_clashes_change', 'entropy_sidechain_change', 'entropy_mainchain_change',
    'sloop_entropy_change', 'mloop_entropy_change', 'cis_bond_change', 'torsional_clash_change',
    'backbone_clash_change', 'helix_dipole_change', 'water_bridge_change', 'disulfide_change',
    'electrostatic_kon_change', 'partial_covalent_bonds_change', 'energy_ionisation_change',
    'entropy_complex_change',
    # 'number_of_residues_change' <-- does not make sense
]

_foldx_interface_features = (
    ['intraclashes_energy_1_wt', 'intraclashes_energy_2_wt',
     'intraclashes_energy_1_change', 'intraclashes_energy_2_change'] +
    _foldx_core_features
)

_physicochem_features = [
    # PhysicoChemical properties
    'pcv_salt_equal_wt', 'pcv_salt_equal_self_wt',
    'pcv_salt_equal_change', 'pcv_salt_equal_self_change',
    'pcv_salt_opposite_wt', 'pcv_salt_opposite_self_wt',
    'pcv_salt_opposite_change', 'pcv_salt_opposite_self_change',
    'pcv_hbond_wt', 'pcv_hbond_self_wt',
    'pcv_hbond_change', 'pcv_hbond_self_change',
    'pcv_vdw_wt', 'pcv_vdw_self_wt',
    'pcv_vdw_change', 'pcv_vdw_self_change',
]

_remaining_features = [
    # Alignment
    'alignment_identity', 'alignment_coverage', 'alignment_score', 'matrix_score',
    # Model
    'norm_dope',
    # Sequence
    'provean_score', 'secondary_structure_wt', 'secondary_structure_change',
    # Structure
    'solvent_accessibility_wt', 'solvent_accessibility_change',
]


def test__get_foldx_features_core():
    expected = _foldx_core_features
    actual = elaspic.elaspic_predictor._get_foldx_features('core')
    xor = set(expected) ^ set(actual)
    assert not xor, xor


def test__get_foldx_features_interface():
    expected = _foldx_interface_features
    actual = elaspic.elaspic_predictor._get_foldx_features('interface')
    xor = set(expected) ^ set(actual)
    assert not xor, xor


def test__get_physicochem_features():
    expected = _physicochem_features
    actual = elaspic.elaspic_predictor._get_physicochem_features()
    xor = set(expected) ^ set(actual)
    assert not xor, xor


def test_feature_columns_core():
    expected = _foldx_core_features + _physicochem_features + _remaining_features
    actual = elaspic.elaspic_predictor.FEATURE_COLUMNS_CORE
    xor = set(expected) ^ set(actual)
    assert not xor, xor


def test_feature_columns_interface():
    expected = _foldx_interface_features + _physicochem_features + _remaining_features
    actual = elaspic.elaspic_predictor.FEATURE_COLUMNS_INTERFACE
    xor = set(expected) ^ set(actual)
    assert not xor, xor


class TestElaspicPredictor:

    @classmethod
    def setup_class(cls):
        cls.predictor = elaspic.Predictor(data_dir=op.abspath(op.splitext(__file__)[0]))

    def test_train(self):
        self.predictor.train()

    def test_score(self):
        self.predictor()
