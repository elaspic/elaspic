import os
import tempfile

import pandas as pd

from sklearn import ensemble
from sklearn import cross_validation
from sklearn import metrics

from common import dat


# %%
def write_row_to_file(results, output_filename):
    """
    TODO: Add a datetime column to each written row.
    """
    print('Saving results to "{}"...'.format(output_filename))
    results_df = pd.DataFrame(results, index=[1])
    results_df = results_df.reindex_axis(sorted(results_df.columns), axis=1)
    print('Results:\n{}'.format(results_df))
    try:
        with dat.open_exclusively(output_filename) as ofh:
            results_df.to_csv(ofh, sep='\t', mode='a', index=False, header=False)
    except Exception as e:
        print('Counld not append result to file: {}'.format(output_filename))
        fd, fd_name = tempfile.mkstemp(
            prefix=output_filename.split('/')[-1],
            suffix='.txt',
            dir='/'.join(output_filename.split('/')[:-1]))
        with os.fdopen(fd, 'w') as ofh:
            results_df.to_csv(ofh, sep='\t', mode='a', index=False, header=False)
        print('Saved it to a separate file instead: {}'.format(fd_name))
        print('The error was: {}\n results_df: {}\n\n'.format(e, results_df))


def cross_validate_predictor(data, features, clf_options, output_filename=None):
    print(clf_options)
    data_x = data[features].values
    data_y = data['ddg_exp'].values
    cv = cross_validation.LeaveOneLabelOut(data['label'].values)
    clf = ensemble.GradientBoostingRegressor(**clf_options)
    y_pred_all = []
    y_true_all = []
    for train, test in cv:
        x_train = data_x[train]
        y_train = data_y[train]
        x_test = data_x[test]
        y_test = data_y[test]
        clf.fit(x_train, y_train)
        probas_ = clf.predict(x_test)
        y_pred_all.extend(probas_)
        y_true_all.extend(y_test)
    results = clf_options.copy()
    results['n_features'] = len(features)
    results['features'] = ','.join(features)
    results['explained_variance_score'] = metrics.explained_variance_score(y_true_all, y_pred_all)
    results['mean_absolute_error'] = metrics.mean_absolute_error(y_true_all, y_pred_all)
    results['mean_squared_error'] = metrics.mean_squared_error(y_true_all, y_pred_all)
    results['r2_score'] = metrics.r2_score(y_true_all, y_pred_all)

    if output_filename is not None:
        write_row_to_file(results, output_filename)
    return results, y_true_all, y_pred_all


def get_final_predictor(data, features, options):
    """Train a predictor using the entire dataset.
    """
    # Keep only recognized options
    import inspect
    accepted_options = inspect.getargspec(ensemble.GradientBoostingRegressor)[0]
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
    clf = ensemble.GradientBoostingRegressor(**clf_options)
    clf.fit(data_x, data_y)
    return clf
