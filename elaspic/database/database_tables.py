import abc
import os.path as op
import tempfile

import sqlalchemy as sa
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

import elaspic
from .database_config import Base, get_binary_collation


# Default sizes for creating varchar fields
SHORT = 15
MEDIUM = 255
LONG = 8192


def decode_domain_def(domains):
    if domains[-1] == ',':
        domains = domains[:-1]
    domain_fragments = [[r.strip() for r in ro.split(':')] for ro in domains.split(',')]
    domain_merged = domain_fragments[0][0], domain_fragments[-1][-1]
    return domain_merged


def get_protein_name():
    raise NotImplementedError


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


class Protein(Base):
    """Protein table.

    Stores information about the protein.
    """
    id = sa.Column(sa.Integer, primary_key=True)
    protein_id = sa.Column(sa.String(MEDIUM), unique=True)
    protein_name = sa.Column(sa.String(MEDIUM), nullable=False)
    protein_description = sa.Column(sa.Text)
    gene_name = sa.Column(sa.String(MEDIUM))
    organism_name = sa.Column(sa.String(MEDIUM))
    protein_sequence = sa.Column(sa.Text, nullable=False)
    protein_sequence_version = sa.Column(sa.Integer)
    provean_supset_filename = sa.Column(sa.Text)
    provean_supset_length = sa.Column(sa.Integer)

    # Private
    _sequence_file = None

    @property
    @abc.abstractmethod
    def root_dir(self):
        raise NotImplementedError

    def __str__(self):
        return '%s (%s)' % (self.id, self.name)

    @property
    def sequence_file(self):
        if self._sequence_file is None:
            secrecord = SeqRecord(id=self.protein_id, seq=Seq(self.protein_sequence))
            sequence_file = op.join(tempfile.tempdir, self.protein_id + '.fasta')
            with open(sequence_file, 'w') as ofh:
                SeqIO.write(secrecord, ofh, 'fasta')
            self._sequence_file = sequence_file
        return self._sequence_file

    @property
    def provean_supset_file(self):
        if self.provean_supset_filename is None:
            return None
        return op.join(self.root_dir, self.provean_supset_filename)

    @provean_supset_file.setter
    def provean_supset_file(self, provean_supset_file):
        provean_supset_filename = op.relpath(provean_supset_file, self.root_dir)
        assert provean_supset_filename == op.basename(provean_supset_filename)
        self.provean_supset_filename = provean_supset_filename

    # === Database ===
    # @property
    # def root_dir(self):
    #     organism_name = self.protein_name.split('_')[-1].lower()
    #     root_dir = op.join(
    #         elaspic.CONFIGS['archive_dir'], get_protein_subfolder(organism_name, self.protein_id))
    #     return root_dir

    # === Local ===
    # class ProteinLocal(_Protein, Base):
    #     @property
    #     def root_dir(self):
    #         return elaspic.CONFIGS['archive_dir']


class DomainModel(Base):
    __indexes = [
        (['protein_id', 'domain_idx'], {'unique': True}),
    ]

    @property
    def local(self):
        raise NotImplementedError

    @property
    def protein(self):
        return self.getprot()

    # Key
    id = sa.Column('domain_id', sa.Integer, primary_key=True)
    protein_id = sa.Column(sa.String(SHORT))
    domain_idx = sa.Column(sa.Integer, index=True)

    # Domain
    clan = sa.Column('pfam_clan', sa.String(MEDIUM))
    name = sa.Column('pdbfam_name', sa.String(MEDIUM), index=True)
    alignment_def = sa.Column('alignment_def', sa.String(MEDIUM))

    data_path = sa.Column('path_to_data', sa.Text)

    def get_protein_name(self, chain=1):
        return get_protein_name(self.protein_id, local=self.local)

    def getclan(self, chain=1):
        return self.clan if self.clan else '-'

    def getname(self, chain=1):
        return self.name

    def getdefs(self, chain=1):
        if self.model_domain_def:
            return self.model_domain_def
        elif self.domain_def:
            return self.domain_def
        return self.alignment_def

    # def __str__(self):
    #     return self.name

    # Template
    # align_file = models.CharField(max_length=255, blank=True, db_column='alignment_filename')
    align_score = sa.Column('alignment_score', sa.Integer)
    align_coverage = sa.Column('alignment_coverage', sa.Float)
    template_errors = sa.Column('template_errors', sa.Text)
    domain_def = sa.Column(sa.String(MEDIUM))
    cath = sa.Column('cath_id', sa.String(MEDIUM, collation=get_binary_collation()))
    seq_id = sa.Column('alignment_identity', sa.Float)

    def getcath(self, chain=1):
        return self.cath[:-2]

    def getSeqId(self, chain=1):
        return self.seq_id

    def getAlnSc(self, chain=1):
        if chain == 1:
            return self.align_score
        elif chain == 2:
            return '-'

    def getsequenceidentity(self, chain=1):
        return '%0.3f' % self.seq_id

    def getalignscore(self, chain=1):
        return '%0.3f' % self.align_score

    @property
    def error(self):
        if self.model_errors:
            return self.model_errors
        if self.template_errors:
            return self.template_errors
        return None

    # def __str__(self):
    #     return '%d' % self.domain_id

    # Model
    model_errors = sa.Column('model_errors', sa.Text)
    dope_score = sa.Column('norm_dope', sa.Float)
    model_filename = sa.Column(sa.String(MEDIUM))
    alignment_filename = sa.Column(sa.String(MEDIUM))
    chain = sa.Column(sa.String(1))
    model_domain_def = sa.Column(sa.String(MEDIUM))

    def getchain(self, chain):
        return self.chain

    def __str__(self):
        return '{}.{}.{}'.format(self.id, self.protein_id, self.domain_idx)

    # === Database ===

    # local = False
    # # interactions = models.ManyToManyField(
    # #     'self', symmetrical=False, through='InterfaceModel', blank=True)
    # # interactions = sa.orm.relationship(CoreModel, secondary=InterfaceModel)
    # # protein = models.ForeignKey(Protein, db_index=True, db_column='protein_id')
    #
    # def getprot(self, chain=1):
    #     return Protein.objects.get(id=self.protein_id)
    #
    # def getpdbtemplate(self, chain=1, link=True):
    #     pdb = self.cath[:-3] + '_' + self.cath[-3]
    #     if link:
    #         return (
    #             '<a class="click2" target="_blank" href="http://www.cathdb.info/pdb/%s">%s</a>'
    #             % (self.cath[:-3], pdb)
    #         )
    #     return pdb

    # === Local ===

    # local = True
    #
    # def getprot(self, chain=1):
    #     return ProteinLocal.objects.get(id=self.protein_id)
    #
    # def getpdbtemplate(self, chain=1, link=True):
    #     return self.cath
    #


class DomainMutation(Base):

    mutation_type = 'core'

    @property
    def protein(self):
        raise NotImplementedError

    # Key
    protein_id = sa.Column(sa.String(SHORT), primary_key=True)
    # domain_id = models.AutoField(primary_key=True, db_column='domain_id', db_index=True)
    domain_idx = sa.Column(sa.Integer, index=True)
    mut = sa.Column('mutation', sa.String(8))

    # Data
    mut_date_modified = sa.Column(sa.Date)

    model_filename_wt = sa.Column(sa.String(MEDIUM))
    model_filename_mut = sa.Column(sa.String(MEDIUM))

    mut_errors = sa.Column('mutation_errors', sa.Text)

    pdb_chain = sa.Column('chain_modeller', sa.String(2, collation=get_binary_collation()))
    pdb_mut = sa.Column('mutation_modeller', sa.String(8))

    stability_energy_wt = sa.Column(sa.Text)
    stability_energy_mut = sa.Column(sa.Text)

    physchem_wt = sa.Column(sa.String(MEDIUM))
    physchem_wt_ownchain = sa.Column(sa.String(MEDIUM))
    physchem_mut = sa.Column(sa.String(MEDIUM))
    physchem_mut_ownchain = sa.Column(sa.String(MEDIUM))

    secondary_structure_wt = sa.Column(sa.String(1))
    secondary_structure_mut = sa.Column(sa.String(1))
    solvent_accessibility_wt = sa.Column(sa.Float)
    solvent_accessibility_mut = sa.Column(sa.Float)

    matrix_score = sa.Column(sa.Float)
    provean_score = sa.Column(sa.Float)

    ddG = sa.Column('ddg', sa.Float)

    elaspic_version = sa.Column(sa.String(SHORT), default=elaspic.__version__)

    def getdomain(self, chain=1):
        return self.model

    def dGwt(self):
        return self.stability_energy_wt.split(',')[0] if self.stability_energy_wt else None

    def dGmut(self):
        return self.stability_energy_mut.split(',')[0] if self.stability_energy_mut else None

    def getddG(self):
        return self.ddG if self.ddG else '-'

    def findChain(self):
        return 1

    def __str__(self):
        return '{}.{}'.format(self.protein_id, self.mut)

    class Meta:
        abstract = True
        # ordering = ['id']
        unique_together = [
            ("protein_id", "model", "mut"),
        ]
        index_together = [
            ("protein_id", "mut"),
        ]

    # === Database ===

    # # model = models.ForeignKey(CoreModel, db_column='domain_id', related_name='muts')
    #
    # @property
    # def protein(self):
    #     return Protein.objects.get(id=self.protein_id)
    #
    # # protein = models.ForeignKey(Protein, db_index=True, db_column='protein_id')

    # === Local ===

    # # model = models.ForeignKey(CoreModelLocal, db_column='domain_id', related_name='muts')
    #
    # @property
    # def protein(self):
    #     return ProteinLocal.objects.get(id=self.protein_id)
    #
    # # protein = models.ForeignKey(ProteinLocal, db_index=True, db_column='protein_id')
    #
    # class Meta(_CoreMutation.Meta):
    #     db_table = 'elaspic_core_mutation_local'


class InterfaceModel(Base):

    @property
    def local(self):
        raise NotImplementedError

    # Key
    id = sa.Column('interface_id', sa.Integer, primary_key=True)
    protein_id_1 = sa.Column(sa.String(SHORT))
    domain_id_1 = sa.Column(sa.Integer, index=True)
    # domain_idx_1 = sa.Column(sa.Integer, index=True)
    protein_id_2 = sa.Column(sa.String(SHORT))
    domain_id_2 = sa.Column(sa.Integer, index=True)
    # domain_idx_2 = sa.Column(sa.Integer, index=True)

    # domain pair
    data_path = sa.Column('path_to_data', sa.Text)

    # def get_protein_name(self, chain=1):
    #     if chain == 1:
    #         return get_protein_name(self.protein_id_1, local=self.local)
    #     elif chain == 2:
    #         return get_protein_name(self.protein_id_2, local=self.local)
    #     else:
    #         raise ValueError

    def getclan(self, chain):
        if chain == 1:
            return self.domain1.getclan()
        elif chain == 2:
            return self.domain2.getclan()

    def getname(self, chain):
        if chain == 1:
            return self.domain1.name
        elif chain == 2:
            return self.domain2.name

    def getdefs(self, chain):
        try:
            if chain == 1:
                defs = self.model_domain_def_1
            elif chain == 2:
                defs = self.model_domain_def_2
            if defs:
                return defs
        except Exception:
            pass
        if chain == 1:
            return self.domain1.getdefs()
        elif chain == 2:
            return self.domain2.getdefs()

    def getdomain(self, chain):
        if chain == 1:
            return self.domain1
        elif chain == 2:
            return self.domain2

    # def __str__(self):
    #     return '%s-%s' % (self.domain1_id, self.domain2_id)

    # domain pair template
    align_score1 = sa.Column('alignment_score_1', sa.Integer)
    align_score2 = sa.Column('alignment_score_2', sa.Integer)

    align_coverage_1 = sa.Column('alignment_identity_1', sa.Integer)
    align_coverage_2 = sa.Column('alignment_coverage_2', sa.Integer)

    cath_id_1 = sa.Column(sa.String(MEDIUM))
    cath_id_2 = sa.Column(sa.String(MEDIUM))

    # seq_id1 = sa.Column('alignment_identity_1', sa.Float)
    # seq_id2 = sa.Column('alignment_identity_2', sa.Float)

    errors = sa.Column('template_errors', sa.Text)

    def getcath(self, chain):
        if chain == 1:
            return self.cath1[:-2]
        elif chain == 2:
            return self.cath2[:-2]

    def getSeqId(self, chain):
        if chain == 1:
            return self.seq_id1
        elif chain == 2:
            return self.seq_id2

    def getAlnSc(self, chain):
        if chain == 1:
            return self.align_score1
        elif chain == 2:
            return self.align_score2

    def getsequenceidentity(self, chain):
        if chain == 1:
            return '%0.3f, %0.3f' % (self.seq_id1, self.seq_id2)
        elif chain == 2:
            return '%0.3f, %0.3f' % (self.seq_id2, self.seq_id1)

    def getalignscore(self, chain):
        if chain == 1:
            return '%0.3f, %0.3f' % (self.align_score1, self.align_score2)
        elif chain == 2:
            return '%0.3f, %0.3f' % (self.align_score2, self.align_score1)

    # def __str__(self):
    #     return '%d' % self.domain_id

    # domain pair model
    model_domain_def_1 = sa.Column(sa.String(MEDIUM))
    model_domain_def_2 = sa.Column(sa.String(MEDIUM))
    model_errors = sa.Column(sa.Text)
    dope_score = sa.Column('norm_dope', sa.Float)
    model_filename = sa.Column(sa.String(MEDIUM))

    alignment_filename_1 = sa.Column(sa.String(MEDIUM))
    alignment_filename_2 = sa.Column(sa.String(MEDIUM))

    aa1 = sa.Column('interacting_aa_1', sa.Text)
    aa2 = sa.Column('interacting_aa_2', sa.Text)

    chain_1 = sa.Column(sa.String(1))
    chain_2 = sa.Column(sa.String(1))

    interface_area_hydrophobic = sa.Column(sa.Float)
    interface_area_hydrophilic = sa.Column(sa.Float)
    interface_area_total = sa.Column(sa.Float)

    def getchain(self, chain):
        if chain == 1:
            return self.chain_1
        elif chain == 2:
            return self.chain_2

    def __str__(self):
        return '{}-{}-{}'.format(self.id, self.protein_id_1, self.protein_id_2)

    class Meta:
        abstract = True
        ordering = ['id']

    # === Database ===

    # local = False
    # # domain1 = models.ForeignKey(
    # #     CoreModel, db_index=True, related_name='p1', db_column='domain_id_1')
    # # domain2 = models.ForeignKey(
    # #     CoreModel, db_index=True, related_name='p2', db_column='domain_id_2')
    #
    # def getprot(self, chain):
    #     if chain == 1:
    #         return Protein.objects.get(id=self.protein_id_1)
    #     elif chain == 2:
    #         return Protein.objects.get(id=self.protein_id_2)
    #
    # def getpdbtemplate(self, chain, link=True):
    #     if link:
    #         a1 = (
    #             '<a class="click2" target="_blank" href="http://www.cathdb.info/pdb/%s">%s_%s</a>'
    #             % (self.cath1[:-3], self.cath1[:-3], self.cath1[-3])
    #         )
    #         a2 = (
    #             '<a class="click2" target="_blank" href="http://www.cathdb.info/pdb/%s">%s_%s</a>'
    #             % (self.cath2[:-3], self.cath2[:-3], self.cath2[-3])
    #         )
    #     else:
    #         a1 = self.cath1[:-3] + '_' + self.cath1[-3]
    #         a2 = self.cath2[:-3] + '_' + self.cath2[-3]
    #     if chain == 1:
    #         return '%s, %s' % (a1, a2)
    #     elif chain == 2:
    #         return '%s, %s' % (a2, a1)

    # === Local ===

    # # local = True
    # # domain1 = models.ForeignKey(
    # #     CoreModelLocal, db_index=True, related_name='p1', db_column='domain_id_1')
    # # domain2 = models.ForeignKey(
    # #     CoreModelLocal, db_index=True, related_name='p2', db_column='domain_id_2')
    # #
    # # domain_1 = sa.orm.relationship(
    # #     CoreModelLocal, uselist=False, cascade='expunge', lazy='joined',
    # #     primaryjoin=(cath_id_1 == Domain.cath_id))  # many to one
    # # domain_2 = sa.orm.relationship(
    # #     CoreModelLocal, uselist=False, cascade='expunge', lazy='joined',
    # #     primaryjoin=(cath_id_2 == Domain.cath_id))  # many to one
    #
    # def getprot(self, chain):
    #     if chain == 1:
    #         return ProteinLocal.objects.get(id=self.protein_id_1)
    #     elif chain == 2:
    #         return ProteinLocal.objects.get(id=self.protein_id_2)
    #
    # def getpdbtemplate(self, chain, link=True):
    #     if chain == 1:
    #         return '{}, {}'.format(self.cath1, self.cath2)
    #     elif chain == 2:
    #         return '{}, {}'.format(self.cath2, self.cath1)
    #
    # class Meta(_InterfaceModel.Meta):
    #     db_table = 'elaspic_interface_model_local'


class InterfaceMutation(Base):

    mutation_type = 'interface'  # interface

    @property
    def protein(self):
        raise NotImplementedError

    # Key
    # interface_id = models.IntegerField(primary_key=True, db_index=True)
    # protein_id_1 = sa.Column(sa.String(SHORT))
    # domain_id_1 = models.IntegerField(db_index=True)
    # protein_id_2 = sa.Column(sa.String(SHORT))
    # domain_id_2 = models.IntegerField(db_index=True)
    id = sa.Column('interface_id', sa.Integer, primary_key=True)
    protein_id = sa.Column(sa.String(SHORT), primary_key=True)
    mut = sa.Column('mutation', sa.String(8))
    chain_idx = sa.Column(sa.Integer)

    # Data
    mut_date_modified = sa.Column(sa.Date)

    model_filename_wt = sa.Column(sa.String(MEDIUM))
    model_filename_mut = sa.Column(sa.String(MEDIUM))

    mut_errors = sa.Column('mutation_errors', sa.Text)

    pdb_chain = sa.Column('chain_modeller', sa.String(2, collation=get_binary_collation()))
    pdb_mut = sa.Column('mutation_modeller', sa.String(8))

    stability_energy_wt = sa.Column(sa.Text)
    stability_energy_mut = sa.Column(sa.Text)
    analyse_complex_energy_wt = sa.Column(sa.Text)
    analyse_complex_energy_mut = sa.Column(sa.Text)

    physchem_wt = sa.Column(sa.String(MEDIUM))
    physchem_wt_ownchain = sa.Column(sa.String(MEDIUM))
    physchem_mut = sa.Column(sa.String(MEDIUM))
    physchem_mut_ownchain = sa.Column(sa.String(MEDIUM))

    secondary_structure_wt = sa.Column(sa.String(1))
    secondary_structure_mut = sa.Column(sa.String(1))
    solvent_accessibility_wt = sa.Column(sa.Float)
    solvent_accessibility_mut = sa.Column(sa.Float)

    contact_distance_wt = sa.Column(sa.Float)
    contact_distance_mut = sa.Column(sa.Float)

    matrix_score = sa.Column(sa.Float)
    provean_score = sa.Column(sa.Float)

    ddG = sa.Column('ddg', sa.Float)
    elaspic_version = sa.Column(sa.String(SHORT), default=elaspic.__version__)

    def getdomain(self, chain):
        return self.model.getdomain(chain)

    def dGwt(self):
        try:
            dGwt = self.analyse_complex_energy_wt.split(',')[2]
        except AttributeError:
            dGwt = None
        return dGwt

    def dGmut(self):
        try:
            dGmut = self.analyse_complex_energy_mut.split(',')[2]
        except AttributeError:
            dGmut = None
        return dGmut

    def getddG(self):
        return self.ddG if self.ddG else '-'

    def findChain(self):
        if self.chain_idx == 0:
            return 1
        elif self.chain_idx == 1:
            return 2
        elif self.protein_id == self.model.protein_id_1:
            return 1
        elif self.protein_id == self.model.protein_id_2:
            return 2
        else:
            raise ValueError(
                'self.chain_idx: {}, self.protein_id: {}'
                .format(self.chain_idx, self.protein_id))

    def getinacprot(self, chain=None):
        c = chain or self.findChain()
        if c == 1:
            return self.model.getprot(2)
        elif c == 2:
            return self.model.getprot(1)

    def __str__(self):
        return '%s.%s' % (self.protein_id, self.mut)

    class Meta:
        abstract = True

    # === Database ===

    # # model = models.ForeignKey(InterfaceModel, db_column='interface_id', related_name='muts')
    # model = sa.orm.relationship(
    #     'interface_id', InterfaceModel, uselist=False, cascade='expunge', lazy='joined',
    #     backref=sa.orm.backref('mutations', cascade='expunge'))  # many to one
    #
    # @property
    # def protein(self):
    #     return Protein.objects.get(id=self.protein_id)
    #
    # # protein = models.ForeignKey(Protein, db_index=True, db_column='protein_id')
    #
    # class Meta(_InterfaceMutation.Meta):
    #     db_table = 'elaspic_interface_mutation'
    #     managed = False

    # === Local ===

    # # model = models.ForeignKey(InterfaceModelLocal, db_column='interface_id', related_name='muts')
    # model = sa.orm.relationship(
    #     'interface_id', InterfaceModelLocal, uselist=False, cascade='expunge', lazy='joined',
    #     backref=sa.orm.backref('mutations', cascade='expunge'))  # many to one
    #
    # # Relationships
    # model = sa.orm.relationship(
    #     InterfaceModelLocal, cascade='expunge', uselist=False, lazy='joined',
    #     backref=sa.orm.backref('mutations', cascade='expunge'))  # many to one
    #
    # protein = sa.orm.relationship(
    #     ProteinLocal, cascade='expunge', uselist=False, lazy='joined',
    #     backfer=sa.orm.backref('mutations', cascade='expunge')
    # )
