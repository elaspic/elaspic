import logging
import os.path as op
import random

import pandas as pd
import pytest

from elaspic import conf

logger = logging.getLogger(__name__)

# Constants
QUICK = False
CONFIG_FILE = op.join(op.dirname(__file__), "test_database_pipeline.ini")

if hasattr(pytest, "config"):
    QUICK = pytest.config.getoption("--quick")
    CONFIG_FILE = pytest.config.getoption("--config-file") or CONFIG_FILE

conf.read_configuration_file(CONFIG_FILE)
assert conf.CONFIGS["db_type"]

logger.debug("Running quick: {}".format(QUICK))
logger.debug("Config file: {}".format(CONFIG_FILE))


# Imports that require a parsed config file
import helper_fns  # noqa

from elaspic import elaspic_database  # noqa

db = elaspic_database.MyDatabase()
conf.CONFIGS["engine"] = db.get_engine()
conf.CONFIGS["engine"].execute("SET sql_mode = ''")


test_cases = []


def append_test_cases(df, num=3, num_mutations=3):
    """.

    Parameters
    ----------
    df : DataFrame
        Contains the following columns:
          - `uniprot_id`
          - `uniprot_sequence`
          - `interacting_aa` OR `model_domain_def` OR `domain_def`
    """
    if QUICK:
        num = 1
    for i in range(num):
        row_idx = random.randint(0, len(df) - 1)
        if df.empty:
            raise Exception("empty dataframe supplied: {}".format(df))
        row = df.iloc[row_idx]
        uniprot_id = row["uniprot_id"]
        uniprot_sequence = row["uniprot_sequence"]
        logger.debug("Protein ID: {}".format(uniprot_id))
        for i in range(num_mutations):
            if "interacting_aa_1" in row and pd.notnull(row["interacting_aa_1"]):
                mutation_pos = random.choice([int(x) for x in row["interacting_aa_1"].split(",")])
                logger.debug("Selected interface AA: {}".format(mutation_pos))
            elif "model_domain_def" in row and pd.notnull(row["model_domain_def"]):
                domain_start, domain_end = (int(x) for x in row["model_domain_def"].split(":"))
                mutation_pos = random.randint(domain_start, domain_end)
                logger.debug(
                    "Selected AA: {} falling inside model domain: {}".format(
                        mutation_pos, row["model_domain_def"]
                    )
                )
            elif "domain_def" in row and pd.notnull(row["domain_def"]):
                domain_start, domain_end = (int(x) for x in row["domain_def"].split(":"))
                mutation_pos = random.randint(domain_start, domain_end)
                logger.debug(
                    "Selected AA: {} falling inside domain: {}".format(
                        mutation_pos, row["domain_def"]
                    )
                )
            else:
                mutation_pos = random.randint(1, len(uniprot_sequence))
            mutation_from = uniprot_sequence[mutation_pos - 1]
            mutation_to = random.choice("GVALICMFWPDESTYQNKRH")
            mutation = "{}{}{}".format(mutation_from, mutation_pos, mutation_to)
            test_cases.append(
                (
                    uniprot_id,
                    mutation,
                )
            )


# %% Everything is missing
sql_query = """
select ud.uniprot_id, us.uniprot_sequence, udt.domain_def
from {db_schema}.uniprot_domain ud
join {db_schema}.uniprot_domain_template udt using (uniprot_domain_id)
join {db_schema}.uniprot_domain_pair udp on (udp.uniprot_domain_id_1 = ud.uniprot_domain_id)
join {db_schema}.uniprot_domain_pair_template udpt using (uniprot_domain_pair_id)
join {db_schema}.uniprot_sequence us using (uniprot_id)
where uniprot_id not in
    (select uniprot_id from {db_schema}.provean)
and uniprot_domain_id not in
    (select uniprot_domain_id from {db_schema}.uniprot_domain_model)
and uniprot_domain_pair_id not in
    (select uniprot_domain_pair_id from {db_schema}.uniprot_domain_pair_model)
and CHAR_LENGTH(us.uniprot_sequence) < 1000
and db = 'sp'
limit 1000;
""".format(
    db_schema=conf.CONFIGS["db_schema"]
)
df = pd.read_sql_query(sql_query, conf.CONFIGS["engine"])

logger.debug("Everything is missing: {}".format(len(df)))
if df.empty:
    logger.error("Skipping...")
else:
    append_test_cases(df)


# %% Have provean and domain model but not interface model
sql_query = """
select ud.uniprot_id, us.uniprot_sequence, udm.model_domain_def
from {db_schema}.uniprot_domain ud
join {db_schema}.provean using (uniprot_id)
join {db_schema}.uniprot_domain_template using (uniprot_domain_id)
join {db_schema}.uniprot_domain_model udm using (uniprot_domain_id)
join {db_schema}.uniprot_domain_pair udp on (udp.uniprot_domain_id_1 = ud.uniprot_domain_id)
join {db_schema}.uniprot_domain_pair_template using (uniprot_domain_pair_id)
join {db_schema}.uniprot_sequence us using (uniprot_id)
where uniprot_domain_pair_id not in
    (select uniprot_domain_pair_id from {db_schema}.uniprot_domain_pair_model)
and CHAR_LENGTH(us.uniprot_sequence) < 1000
and db = 'sp'
limit 1000;
""".format(
    db_schema=conf.CONFIGS["db_schema"]
)
df = pd.read_sql_query(sql_query, conf.CONFIGS["engine"])

logger.debug("Have provean and domain model but not interface model: {}".format(len(df)))
if df.empty:
    logger.error("Skipping...")
elif not QUICK:
    append_test_cases(df)


# %% Have provean and all models
sql_query = """
select ud.uniprot_id, us.uniprot_sequence, udpm.interacting_aa_1
from {db_schema}.uniprot_domain ud
join {db_schema}.provean using (uniprot_id)
join {db_schema}.uniprot_domain_template using (uniprot_domain_id)
join {db_schema}.uniprot_domain_model using (uniprot_domain_id)
join {db_schema}.uniprot_domain_pair udp on (udp.uniprot_domain_id_1 = ud.uniprot_domain_id)
join {db_schema}.uniprot_domain_pair_template using (uniprot_domain_pair_id)
join {db_schema}.uniprot_domain_pair_model udpm using (uniprot_domain_pair_id)
join {db_schema}.uniprot_sequence us using (uniprot_id)
where CHAR_LENGTH(us.uniprot_sequence) < 1000
and db = 'sp'
and udpm.model_filename is not null
and udpm.model_errors is null
limit 1000;
""".format(
    db_schema=conf.CONFIGS["db_schema"]
)
df = pd.read_sql_query(sql_query, conf.CONFIGS["engine"])

logger.debug("Have provean and all models: {}".format(len(df)))
if df.empty:
    logger.error("Skipping...")
else:
    append_test_cases(df, 10)


# %% Fixtures
@pytest.fixture(scope="session", params=test_cases)
def uniprot_id_mutation(request):
    return request.param


# %% Tests
def test_database_pipeline(uniprot_id_mutation):
    return helper_fns.run_database_pipeline(uniprot_id_mutation)


# %%
if __name__ == "__main__":
    import pytest

    pytest.main([__file__, "-svx", "--quick"])
