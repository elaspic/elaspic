import sqlalchemy as sa
from sqlalchemy.ext.declarative import declarative_base, declared_attr

import elaspic
from kmtools import df_tools

naming_convention = {
    "ix": 'ix_%(column_0_label)s',
    "uq": "uq_%(table_name)s_%(column_0_name)s",
    "ck": "ck_%(table_name)s_%(constraint_name)s",
    "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
    "pk": "pk_%(table_name)s"
}

#: Some database-specific parameters that SQLAlchemy can't figure out
db_specific_properties = {
    'mysql': {
        'BINARY_COLLATION': 'utf8_bin',
        'STRING_COLLATION': 'utf8_unicode_ci',
    },
    'postgresql': {
        'BINARY_COLLATION': 'en_US.utf8',
        'STRING_COLLATION': 'en_US.utf8',
    },
    'sqlite': {
        'BINARY_COLLATION': 'RTRIM',
        'STRING_COLLATION': 'NOCASE',
    },
}


if elaspic.CONFIGS.get('db_type') is None:
    print('The `DB_TYPE` has not been set. Do not know what database is being used!')


def get_protein_name():
    raise NotImplementedError


def get_binary_collation():
    return db_specific_properties[elaspic.CONFIGS.get('db_type', 'mysql')]['BINARY_COLLATION']


class MyMixin:
    __indexes__ = []

    @declared_attr
    def __tablename__(cls):
        print("Class name: {}".format(cls.__name__))
        if cls.__name__[0] == '_':
            return 'hidden_' + df_tools.format_column(cls.__name__)
        else:
            return df_tools.format_column(cls.__name__)

    @declared_attr
    def __table_args__(cls):
        """Return a tuple of additional table arguments."""
        table_name = cls.__tablename__
        index_columns = cls.__indexes__
        db_specific_params = []
        table_args = []
        # Create indexes over several columns
        for columns in index_columns:
            if type(columns) == tuple:
                column_names, kwargs = columns
            elif type(columns) == list:
                column_names = columns
                kwargs = {}
                kwargs['unique'] = False
            if 'index_name' in kwargs:
                index_name = kwargs.pop('index_name')
            else:
                index_name = (
                    'ix_{table_name}_{column_0_name}'
                    .format(table_name=table_name, column_0_name=column_names[0])[:255]
                )
            table_args.append(sa.Index(index_name, *column_names, **kwargs))
        # Other table parameters, such as schemas, etc.
        for db_specific_param in db_specific_params:
            table_args.append(
                db_specific_properties[elaspic.CONFIGS['db_type']][db_specific_param])
        table_args.append({'mysql_engine': 'InnoDB'})  # this has to be last
        print(table_args)
        return tuple(table_args)


Base = declarative_base(cls=MyMixin)
Base.metadata.naming_conventions = naming_convention
