import logging
import os
import os.path as op
import tempfile
from contextlib import contextmanager

# import abc
import sqlalchemy as sa
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from sqlalchemy.ext.declarative import declarative_base, declared_attr

import elaspic
from kmtools import df_tools, system_tools, py_tools

# Default sizes for creating varchar fields
SHORT = 15
MEDIUM = 255
LONG = 8192

logger = logging.getLogger(__name__)

naming_convention = {
    "ix": 'ix_%(column_0_label)s',
    "uq": "uq_%(table_name)s_%(column_0_name)s",
    "ck": "ck_%(table_name)s_%(constraint_name)s",
    "fk": "fk_%(table_name)s_%(column_0_name)s_%(referred_table_name)s",
    "pk": "pk_%(table_name)s"
}

# Some database-specific parameters that SQLAlchemy can't figure out
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


class Base(object):
    __indexes__ = []

    @declared_attr
    def __tablename__(cls):
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


Base = declarative_base(cls=Base)
Session = sa.orm.sessionmaker(expire_on_commit=False)


def get_protein_subfolder(organism_name, protein_id):
    """Return the name of the subfolder for storing protein information.

    TODO: Screw the splitting of uniprot ids!

    Examples
    --------
    >>> str(get_protein_subfolder('human', 'P27361'))
    'human/P27/36/P27361'
    """
    return op.join(organism_name, protein_id[:3], protein_id[3:5], protein_id)


class _View:
    @property
    def _create_view_sql(self):
        sql_file = op.join(op.abspath(__file__), 'sql', self.__tablename__)
        with open(sql_file, 'rt') as ifh:
            sql = ifh.read().strip()
        return sql

    @property
    def _drop_view_sql(self):
        return "DROP VIEW {};".format(self.__tablename__)


class _Protein(elaspic.Sequence):
    __indexes__ = [(('protein_id', 'organism_name'), {'unique': True}), ]
    protein_id = sa.Column(sa.String(MEDIUM), primary_key=True)
    protein_name = sa.Column(sa.String(MEDIUM), nullable=False)
    protein_description = sa.Column(sa.Text)
    gene_name = sa.Column(sa.String(MEDIUM))
    organism_name = sa.Column(sa.String(MEDIUM))
    protein_sequence = sa.Column(sa.Text, nullable=False)
    protein_sequence_version = sa.Column(sa.Integer)
    _results = sa.Column(sa.JSON)
    # provean_supset_filename = sa.Column(sa.Text)
    # provean_supset_length = sa.Column(sa.Integer)

    # Private
    _root_dir = None
    _sequence_file = None

    def __init__(self):
        results = {
            'protein_sequence_file'
        }
        self._sequence = elaspic.Sequence(self.protein_sequence_file, self.provean_supset_file)

    def __str__(self):
        return '%s (%s)' % (self.protein_id, self.protein_name)

    @property
    def protein_sequence_file(self):
        if self._sequence_file is None:
            secrecord = SeqRecord(id=self.protein_id, seq=Seq(self.protein_sequence))
            sequence_file = op.join(tempfile.tempdir, self.protein_id + '.fasta')
            with open(sequence_file, 'w') as ofh:
                SeqIO.write(secrecord, ofh, 'fasta')
            self._sequence_file = sequence_file
        return self._sequence_file

    @property
    def root_dir(self):
        raise NotImplementedError

    @property
    def results(self):
        results = self._results.copy()
        for method in results:
            for key, value in results[method].items():
                if key.endswith('_file'):
                    results[method][key] = op.join(self.root_dir, value)

    @results.setter
    def results(self, results):
        results = results.copy()
        for method in results:
            for key, value in results[method].items():
                if key.endswith('_file'):
                    if not value.startswith(self._root_dir):
                        system_tools.copyfile(
                            value, op.join(self.root_dir, op.basename(value)), mode=0o666)
                        results[method][key] = op.basename(value)
                    else:
                        system_tools[method][key] = op.relpath(value)

    @property
    def provean_supset_file(self):
        if self.provean_supset_filename:
            provean_supset_file = op.join(self.root_dir, self.provean_supset_filename)
            assert op.isfile(provean_supset_file)
            return provean_supset_file
        else:
            provean_supset_file = op.join(
                self.root_dir, system_tools.slugify(self.protein_id + '_provean_supset'))
            return provean_supset_file

    def build(self):
        if self._sequence.is_done:
            logger.debug("Already done.")
        else:
            self._sequence.run()
            self._results = self._sequence.results
            self.provean_supset_filename = sequence.provean_supset_filename
            self.provean_supset_length = sequence.provean_supset_length
            assert self._sequence.is_done
            with session_scope() as session:
                session.merge(self)

    def mutate(self, mutation):
        ...

class Protein(Base, _Protein, _View):
    @property
    def root_dir(self):
        if self._root_dir is not None:
            return self._root_dir
        organism_name = self.protein_name.split('_')[-1].lower()
        self._root_dir = op.join(
            elaspic.CONFIGS['archive_dir'], get_protein_subfolder(organism_name, self.protein_id))
        os.makedirs(self._root_dir, exist_ok=True)
        return self._root_dir

    def mutate(self):
        mutation = super().mutate()
        mutation = ProteinMutation(**mutation)
        self.mutations.add(mutation)


class LocalProtein(Base, _Protein):
    @property
    def root_dir(self):
        if self._root_dir is not None:
            return self._root_dir
        self._root_dir = elaspic.CONFIGS['archive_dir']
        os.makedirs(self._root_dir, exist_ok=True)
        return self._root_dir

    def mutate(self):
        mutation = super().mutate()
        mutation = LocalProteinMutaton(**mutation)
        self.mutations.add(mutation)


@contextmanager
def session_scope():
    """Provide a transactional scope around a series of operations.

    Enables the following construct: ``with self.session_scope() as session:``.
    """
    session = Session()
    try:
        yield session
        session.commit()
    except:
        session.rollback()
        raise
    finally:
        session.expunge_all()
        session.close()


if __name__ == '__main__':
    logging.basicConfig(level=logging.DEBUG)
    connection_string = (
        'mysql://root@localhost:8888/datapkg?unix_socket=/tmp/strokach/datapkg/mysql.sock'
    )
    engine = sa.create_engine(connection_string)
    Base.metadata.create_all(engine)
    #
    Session.configure(bind=engine)
    elaspic.conf.read_configuration_file(
        EXTERNAL_DIRS={'blast_db_dir': '/home/kimlab1/database_data/blast/db'},
        DATABASE={'db_connection_string': connection_string})

    # with session_scope() as session:
    for i in range(100, 1100):
        p = Protein(
            protein_id='P{}'.format(i),
            protein_name='P{}_human'.format(i),
            protein_sequence='GGGGGGGGG')
        p.run()
        # session.merge(p)

        pl = LocalProtein(
            protein_id='Q{}'.format(i),
            protein_name='Q{}_human'.format(i),
            protein_sequence='AAAAA')
        # session.merge(pl)
