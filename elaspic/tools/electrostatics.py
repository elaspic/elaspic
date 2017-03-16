from elaspic.tools._abc import StructureAnalyzer, ToolError


class ElectrostaticsError(ToolError):
    pass


class Electrostatics(StructureAnalyzer):
    """Calculate electrostatic energy."""

    system_command = """\
pdb2pqr_cli --ff=parse --chain --apbs-input --ph-calc-method={ph_calc_method} \
{pdb_file} {pdr_file} \
"""

    def __init__(self, ff='parse', ph_calc_method='propka31'):
        """
        Parameters
        ----------
        ff : str
            Forcefield to use.
            From pdb2pqr docs: "the PARSE, TYL06, and SWANSON force fields are optimized for
            continuum calculations and, as a result, are probably the best to use".
        ph-calc-method : str
            One of ``propka``, ``propka31``, or ``pdb2pka``.
        """
        return None


if __name__ == '__main__':
    import fire
    fire.Fire()
