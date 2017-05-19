import pytest
from elaspic import pipeline, errors


class TestPipeline:

    def setup_method(self, method):
        self.p = pipeline.Pipeline(None)

    @pytest.mark.parametrize("mutations_in, mutations_out", [
        ('M1C', ['M1C']),
        ('M1C,M1D', ['M1C', 'M1D']),
        ('M1C:M1D', ['M1C', 'M1D']),
        ('A_M1C,A_M1D', ['A_M1C', 'A_M1D']),
        ('0_M1C:0_M1D', ['0_M1C', '0_M1D']),
    ])
    def test_split_mutations(self, mutations_in, mutations_out):
        assert self.p._split_mutations(mutations_in) == mutations_out

    @pytest.mark.parametrize("run_type_in, run_type_out", [
        ('1', 'sequence'),
        ('2', 'model'),
        ('3', 'mutation'),
        ('4', 'model.mutation'),
        ('5', 'sequence.model.mutation'),
        ('6', 'sequence.model'),
        ('all', 'sequence.model.mutation'),
    ])
    def test_validate_run_type(self, run_type_in, run_type_out):
        assert self.p._validate_run_type(run_type_in) == run_type_out

    @pytest.mark.parametrize("run_type", [
        '1a', '2b', 'xxx', 'sequences', 'models', 'mutations'
    ])
    def test_validate_run_type_2(self, run_type):
        with pytest.raises(errors.ParameterError):
            self.p._validate_run_type(run_type)
