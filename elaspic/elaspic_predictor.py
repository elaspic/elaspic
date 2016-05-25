import os.path as op
import pickle
import logging
import json

import pandas as pd

from . import conf, call_foldx

logger = logging.getLogger(__name__)
configs = conf.Configs()


# %%
secondary_structure_to_int = {
    '-': 0,  # Coil (none of the below)
    'C': 0,  # Coil (none of the below)
    'B': 1,  # Isolated bridge
    'b': 1,  # Isolated bridge
    'E': 2,  # Extended conformation
    'G': 3,  # 3-10 helix
    'H': 4,  # Alpha helix
    'I': 5,
    'S': 6,
    'T': 7,  # Turn
}


def _get_foldx_features(core_or_interface):
    """
    """
    feature_columns = []
    # FoldX
    if core_or_interface in [False, 0, 'core']:
        feature_columns += call_foldx.names_stability_wt
        feature_columns += [
            c[:-4] + '_change'
            for c in call_foldx.names_stability_mut
            if c.endswith('_mut')
        ]
    else:
        feature_columns += call_foldx.names_stability_complex_wt
        feature_columns += [
            c[:-4] + '_change'
            for c in call_foldx.names_stability_complex_mut
            if c.endswith('_mut')
        ]
    return feature_columns


def _get_physicochem_features():
    # PhysicoChemical properties
    names_phys_chem = ['pcv_salt_equal', 'pcv_salt_opposite', 'pcv_hbond', 'pcv_vdw']
    feature_columns = []
    feature_columns += [(c + '_wt') for c in names_phys_chem]
    feature_columns += [(c + '_self_wt') for c in names_phys_chem]
    feature_columns += [(c + '_change') for c in names_phys_chem]
    feature_columns += [(c + '_self_change') for c in names_phys_chem]
    return feature_columns


def _get_remaining_features():
    feature_columns = []
    # Sequence
    feature_columns += ['provean_score', 'secondary_structure_wt', 'secondary_structure_change']
    # Alignment
    feature_columns += [
        'alignment_identity', 'alignment_coverage', 'alignment_score', 'matrix_score']
    # Model
    feature_columns += ['norm_dope']
    # Structure
    feature_columns += ['solvent_accessibility_wt', 'solvent_accessibility_change']
    return feature_columns


FEATURE_COLUMNS_CORE = (
    _get_foldx_features('core') + _get_physicochem_features() + _get_remaining_features()
)


FEATURE_COLUMNS_INTERFACE = (
    _get_foldx_features('interface') + _get_physicochem_features() + _get_remaining_features()
)


def format_mutation_features(feature_df, core_or_interface):
    """.

    Converts columns containing comma-separated lists of FoldX features and physicochemical
    features into a DataFrame where each feature has its own column.

    Parameters
    ----------
    feature_df : DataFrame
        A pandas DataFrame containing a subset of rows from the :ref:`uniprot_domain_mutation`
        or the :ref:`uniprot_domain_pair_mutation` tables.
    core_or_interface : int or str
        If 0 or 'core', the `feature_df` DataFrame contains columns from the
        :ref:`uniprot_domain_mutation` table.
        If 1 or 'interface, the feature_df DataFrame contains columns from the
        :ref:`uniprot_domain_pair_mutation` table.

    Returns
    -------
    DataFrame
        Contains the same data as `feature_df`, but with columns containing comma-separated lists
        of features converted to columns containing a single feature each.

    """
    if core_or_interface in [False, 0, 'core']:
        foldx_column_name = 'stability_energy'
        foldx_feature_names_wt = call_foldx.names_stability_wt
        foldx_feature_names_mut = call_foldx.names_stability_mut
    elif core_or_interface in [True, 1, 'interface']:
        foldx_column_name = 'analyse_complex_energy'
        foldx_feature_names_wt = call_foldx.names_stability_complex_wt
        foldx_feature_names_mut = call_foldx.names_stability_complex_mut

    # Drop rows that have missing FoldX information
    # (should not happen when callced from inside the pipeline because we have only one column)
    # feature_df = feature_df.dropna(subset=[foldx_column_name + '_wt', foldx_column_name + '_mut'])

    # FoldX output
    for column_index, column_name in enumerate(foldx_feature_names_wt):
        feature_df[column_name] = feature_df[foldx_column_name + '_wt'].apply(
            lambda x: float(x.split(',')[column_index]) if pd.notnull(x) else x)
    del feature_df[foldx_column_name + '_wt']

    for column_index, column_name in enumerate(foldx_feature_names_mut):
        feature_df[column_name] = feature_df[foldx_column_name + '_mut'].apply(
            lambda x: float(x.split(',')[column_index]) if pd.notnull(x) else x)
    del feature_df[foldx_column_name + '_mut']

    # PhysicoChemical properties
    names_phys_chem = ['pcv_salt_equal', 'pcv_salt_opposite', 'pcv_hbond', 'pcv_vdw']
    for column_index, column_name in enumerate(names_phys_chem):
        feature_df[column_name + '_wt'] = (
            feature_df['physchem_wt']
            .apply(lambda x: int(x.split(',')[column_index]) if pd.notnull(x) else x)
        )
        feature_df[column_name + '_self_wt'] = (
            feature_df['physchem_wt_ownchain']
            .apply(lambda x: int(x.split(',')[column_index]) if pd.notnull(x) else x)
        )
        feature_df[column_name + '_mut'] = (
            feature_df['physchem_mut']
            .apply(lambda x: int(x.split(',')[column_index]) if pd.notnull(x) else x)
        )
        feature_df[column_name + '_self_mut'] = (
            feature_df['physchem_mut_ownchain']
            .apply(lambda x: int(x.split(',')[column_index]) if pd.notnull(x) else x)
        )
    del feature_df['physchem_wt']
    del feature_df['physchem_wt_ownchain']
    del feature_df['physchem_mut']
    del feature_df['physchem_mut_ownchain']

    for col in feature_df.columns:
        if 'secondary_structure' in col:
            feature_df[col] = (
                feature_df[col]
                .apply(lambda x: secondary_structure_to_int[x] if pd.notnull(x) else x)
            )
    return feature_df


def convert_features_to_differences(df, keep_mut=False):
    """
    Creates a new set of features (ending in `_change`) that describe the difference between values
    of the wildtype (features ending in `_wt`) and mutant (features ending in `_mut`) features.
    If `keep_mut` is `False`, removes all mutant features (features ending in `_mut`).
    """
    column_list = []
    for column_name, column in df.iteritems():
        if ('_mut' in column_name and
                column_name.replace('_mut', '_wt') in df.columns and
                df[column_name].dtype != object):
            if keep_mut:
                column_list.append(column)
            new_column = column - df[column_name.replace('_mut', '_wt')]
            if 'secondary_structure' in column_name:
                new_column = new_column.apply(lambda x: 1 if x else 0)
            new_column.name = column_name.replace('_mut', '_change')
            column_list.append(new_column)
        else:
            column_list.append(column)
#    new_df = pd.DataFrame(column_list).T
    new_df = pd.concat(column_list, axis=1)
    return new_df


# %%
class Predictor:

    feature_name_conversion = {
        'normDOPE': 'norm_dope',
        'seq_id_avg': 'alignment_identity'
    }

    def __init__(self):

        def _load_data(filename):
            if op.splitext(filename)[-1] in ['.pkl', '.pickle']:
                with open(op.join(configs['data_dir'], filename), 'rb') as ifh:
                    return pickle.load(ifh)
            elif op.splitext(filename)[-1] in ['.jsn', '.json']:
                with open(op.join(configs['data_dir'], filename), 'r') as ifh:
                    return json.load(ifh)

        self.clf_domain = _load_data('ml_clf_core_p1.pickle')
        self.clf_domain_features = _load_data('ml_features_core_p1.json')
        self.clf_interface = _load_data('ml_clf_interface_p1.pickle')
        self.clf_interface_features = _load_data('ml_features_interface_p1.json')

        self.clf_domain_p1 = _load_data('ml_clf_core_p1.pickle')
        self.clf_domain_features_p1 = _load_data('ml_features_core_p1.json')
        self.clf_interface_p1 = _load_data('ml_clf_interface_p1.pickle')
        self.clf_interface_features_p1 = _load_data('ml_features_interface_p1.json')

    def score(self, df, core_or_interface):
        """
        Parameters
        ----------
        df : DataFrame
            One or more rows with all data required to predict $\Delta \Delta G$ score.
            Like something that you would get when you join the appropriate rows in the database.

        Returns
        -------
        df : Dataframe
            Same as the input dataframe, except with one additional column: `ddg`.
        """
        if core_or_interface in ['core', 0]:
            clf = self.clf_domain
            clf_features = self.clf_domain_features
        elif core_or_interface in ['interface', 1]:
            clf = self.clf_interface
            clf_features = self.clf_interface_features

        feature_name_conversion = {
            'normDOPE': 'norm_dope',
            'seq_id_avg': 'alignment_identity'}
        clf_features = [feature_name_conversion.get(x, x) for x in clf_features]

        df_features = format_mutation_features(df, core_or_interface)
        # keep mut, remove it in next step
        df_features_asdifferences = convert_features_to_differences(df_features, True)
        df_features_asdifferences = df_features_asdifferences[clf_features]

        ddg = clf.predict(df_features_asdifferences)[0]

        return ddg
