import os.path as op
import pickle
import logging
import json

import numpy as np
import pandas as pd
import sklearn.ensemble

from . import call_foldx

logger = logging.getLogger(__name__)

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


def format_mutation_features(df):
    """.

    Converts columns containing comma-separated lists of FoldX features and physicochemical
    features into a DataFrame where each feature has its own column.

    Parameters
    ----------
    feature_df : DataFrame
        A pandas DataFrame containing a subset of rows from the :ref:`uniprot_domain_mutation`
        or the :ref:`uniprot_domain_pair_mutation` tables.

    Returns
    -------
    DataFrame
        Contains the same data as `feature_df`, but with columns containing comma-separated lists
        of features converted to columns containing a single feature each.

    """
    df = df.copy()
    if 'analyse_complex_energy_wt' in df.columns:
        foldx_column_name = 'analyse_complex_energy'
        foldx_feature_names_wt = call_foldx.names_stability_complex_wt
        foldx_feature_names_mut = call_foldx.names_stability_complex_mut
    else:
        foldx_column_name = 'stability_energy'
        foldx_feature_names_wt = call_foldx.names_stability_wt
        foldx_feature_names_mut = call_foldx.names_stability_mut

    # Drop rows that have missing FoldX information
    # (should not happen when callced from inside the pipeline because we have only one column)
    # feature_df = feature_df.dropna(
    #   subset=[foldx_column_name + '_wt', foldx_column_name + '_mut'])

    # FoldX output
    for column_index, column_name in enumerate(foldx_feature_names_wt):
        df[column_name] = df[foldx_column_name + '_wt'].apply(
            lambda x: float(x.split(',')[column_index]) if pd.notnull(x) else x)
    del df[foldx_column_name + '_wt']

    for column_index, column_name in enumerate(foldx_feature_names_mut):
        df[column_name] = df[foldx_column_name + '_mut'].apply(
            lambda x: float(x.split(',')[column_index]) if pd.notnull(x) else x)
    del df[foldx_column_name + '_mut']

    # PhysicoChemical properties
    names_phys_chem = ['pcv_salt_equal', 'pcv_salt_opposite', 'pcv_hbond', 'pcv_vdw']
    for column_index, column_name in enumerate(names_phys_chem):
        df[column_name + '_wt'] = (
            df['physchem_wt']
            .apply(lambda x: int(x.split(',')[column_index]) if pd.notnull(x) else x)
        )
        df[column_name + '_self_wt'] = (
            df['physchem_wt_ownchain']
            .apply(lambda x: int(x.split(',')[column_index]) if pd.notnull(x) else x)
        )
        df[column_name + '_mut'] = (
            df['physchem_mut']
            .apply(lambda x: int(x.split(',')[column_index]) if pd.notnull(x) else x)
        )
        df[column_name + '_self_mut'] = (
            df['physchem_mut_ownchain']
            .apply(lambda x: int(x.split(',')[column_index]) if pd.notnull(x) else x)
        )
    del df['physchem_wt']
    del df['physchem_wt_ownchain']
    del df['physchem_mut']
    del df['physchem_mut_ownchain']

    for col in df.columns:
        if 'secondary_structure' in col:
            df[col] = (
                df[col]
                .apply(lambda x: secondary_structure_to_int[x] if pd.notnull(x) else x)
            )
    return df


def convert_features_to_differences(df, keep_mut=False):
    """Convert `_wt` and `_mut` columns into `_wt` and `_change` columns.

    Create a new set of features (ending in `_change`) that describe the difference between values
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
    new_df = pd.concat(column_list, axis=1)
    return new_df


def get_core_mutations(df, engine, schema_name='elaspic'):
    values = (
        ("('" + df['uniprot_id'] + "', '" + df['uniprot_mutation'] + "')")
    )
    sql_query = """\
SELECT *
FROM {schema_name}.uniprot_domain
JOIN {schema_name}.uniprot_domain_template USING (uniprot_domain_id)
JOIN {schema_name}.uniprot_domain_model USING (uniprot_domain_id)
JOIN {schema_name}.uniprot_domain_mutation mut USING (uniprot_domain_id)
WHERE (mut.uniprot_id, mut.mutation) in ({values})
""".format(schema_name=schema_name, values=', '.join(values))

    results_df = pd.read_sql_query(sql_query, engine)

    # Format predictor features
    results_df = format_mutation_features(results_df)
    # keep `_wt` and `_mut`, remove later:
    results_df = convert_features_to_differences(results_df, True)

    return results_df


def get_interface_mutations(df, engine, schema_name='elaspic'):
    values = (
        ("('" + df['uniprot_id'] + "', '" + df['uniprot_mutation'] + "')")
    )
    sql_query = """\
SELECT *
FROM {schema_name}.uniprot_domain_pair
JOIN {schema_name}.uniprot_domain_pair_template USING (uniprot_domain_pair_id)
JOIN {schema_name}.uniprot_domain_pair_model USING (uniprot_domain_pair_id)
JOIN {schema_name}.uniprot_domain_pair_mutation mut USING (uniprot_domain_pair_id)
WHERE (mut.uniprot_id, mut.mutation) in ({values})
""".format(schema_name=schema_name, values=', '.join(values))

    # Format alignment features
    results_df = pd.read_sql_query(sql_query, engine)
    results_df['alignment_identity'] = (
        np.sqrt(results_df['identical_1'] * results_df['identical_2'])
    )
    results_df['alignment_coverage'] = (
        np.sqrt(results_df['coverage_1'] * results_df['coverage_2'])
    )
    results_df['alignment_score'] = (
        np.sqrt(results_df['score_1'] * results_df['score_2'])
    )

    # Format predictor features
    results_df = format_mutation_features(results_df)
    # keep `_wt` and `_mut`, remove later:
    results_df = convert_features_to_differences(results_df, True)

    return results_df


def get_final_predictor(data, features, options):
    """Train a predictor using the entire dataset."""
    CLF = sklearn.ensemble.GradientBoostingRegressor

    # Keep only recognized options
    import inspect
    accepted_options = inspect.getargspec(CLF)[0]
    clf_options = {k: v for (k, v) in options.items() if k in accepted_options}
    if clf_options != options:
        extra_options = {k: v for (k, v) in options.items() if k not in accepted_options}
        print("Warning, unknown options provided:\n{}".format(extra_options))

    # Remove rows with NULLs
    data_usable = data[features + ['ddg_exp']].dropna()
    if len(data) != len(data_usable):
        print('Warning, {} rows in the provided data contained null!'
              .format(len(data) - len(data_usable)))
    data = data_usable

    # Train predictor
    data_x = data[features].values
    data_y = data['ddg_exp'].values
    clf = CLF(**clf_options)
    clf.fit(data_x, data_y)
    return clf


class _Predictor:

    def __init__(self):
        self.clf = None
        self.features = None

    def _assert_trained(self):
        if self.clf is None:
            raise Exception(
                "You must either load or train a classifier before you can score mutations!")

    def load(self, data_dir):
        """Load predictor data from files in `data_dir`."""
        with open(op.join(data_dir, self.clf_filename), 'rb') as ifh:
            self.clf = pickle.load(ifh)
        with open(op.join(data_dir, self.features_filename), 'rt') as ifh:
            self.features = json.load(ifh)

    def save(self, data_dir):
        """Save predictor data to files in `data_dir`."""
        self._assert_trained()
        with open(op.join(data_dir, self.clf_filename), 'wb') as ofh:
            pickle.dump(self.clf, ofh)
        with open(op.join(data_dir, self.features_filename), 'wt') as ofh:
            json.dump(self.features, ofh)

    def train(self, df, options):
        features = options['features']
        if ',' in features:
            features = features.split(',')
        elif '.' in features:
            features = features.split('.')

        self.clf = get_final_predictor(df, features, options)
        self.features = features

    def score(self, df):
        """Predict ΔΔG score.

        Parameters
        ----------
        df : DataFrame
            One or more rows with all data required to predict ΔΔG score.
            Like something that you would get when you join the appropriate rows in the database.

        Returns
        -------
        df : Dataframe
            Same as the input dataframe, except with one additional column: `ddg`.
        """
        self._assert_trained()
        df = format_mutation_features(df)
        df = convert_features_to_differences(df, True)  # keep mut, remove it in next step
        ddg = self.clf.predict(df[self.features])
        return ddg


class CorePredictor(_Predictor):

    clf_filename = 'ml_clf_core.pickle'
    features_filename = 'ml_features_core.json'


class InterfacePredictor(_Predictor):

    clf_filename = 'ml_clf_interface.pickle'
    features_filename = 'ml_features_interface.json'
