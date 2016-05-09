import os.path as op
import pandas as pd
import elaspic.elaspic_predictor


def test_shape_in_is_shape_out():
    """Make sure that formatting features does not change the size of the DataFrame.

    i.e. that you are not removing for example rows with NaNs
    """
    df = pd.read_csv(op.join(op.splitext(__file__)[0], 'df1.tsv'), sep='\t')
    df_out = elaspic.elaspic_predictor.format_mutation_features(df, 'core')
    assert df.shape[0] == df_out.shape[0]
    df_out = elaspic.elaspic_predictor.convert_features_to_differences(df_out)
    assert df.shape[0] == df_out.shape[0]
