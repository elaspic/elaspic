import datetime
import sqlalchemy as sa
import sqlalchemy.ext.declarative as sa_ext_declarative

from . import conf


# Default sizes for creating varchar fields
SHORT = 15
MEDIUM = 255
LONG = 8192

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

if conf.CONFIGS.get('db_type') is None:
    print('The `DB_TYPE` has not been set. Do not know what database is being used!')


def get_table_args(table_name, index_columns=[], db_specific_params=[]):
    """Return a tuple of additional table arguments."""
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
        table_args.append(db_specific_properties[conf.CONFIGS['db_type']][db_specific_param])
    table_args.append({'mysql_engine': 'InnoDB'})  # this has to be last
    return tuple(table_args)


Base = sa_ext_declarative.declarative_base()
Base.metadata.naming_conventions = naming_convention


class Domain(Base):
    """Profs domain definitions for all proteins in the PDB.

    Columns:
      cath_id
        Unique id identifying each domain in the PDB. Constructed by concatenating the pdb_id,
        pdb_chain, and an index specifying the order of the domain in the chain.

      pdb_id
        The PDB id in which the domain is found.

      pdb_chain
        The PDB chain in which the domain is found.

      pdb_domain_def
        Domain definitions of the domain, in PDB RESNUM coordinates.

      pdb_pdbfam_name
        The Profs name of the domain.

      pdb_pdbfam_idx
        An integer specifying the number of times a domain with domain name ``pdb_pdbfam_name`` has
        occurred in this chain up to this point. It is used to make every
        ``(pdb_id, pdb_chain, pdb_pdbfam_name, pdb_pdbfam_idx)`` tuple unique.

      domain_errors
        List of errors that occurred when annotating this domain, or when using this domain
        to make structural homology models.
    """

    __tablename__ = 'domain'
    if conf.CONFIGS.get('db_type') == 'mysql':
        # MySQL can't handle long indexes
        _indexes = [
            ['pdb_id', 'pdb_chain'],
            (['pdb_pdbfam_name'], {'mysql_length': 255}),
        ]
    else:
        _indexes = [
            (['pdb_id', 'pdb_chain', 'pdb_pdbfam_name', 'pdb_pdbfam_idx'], {'unique': True}),
            (['pdb_pdbfam_name'], {'mysql_length': 255})
        ]
    __table_args__ = get_table_args(__tablename__, _indexes, [])
    cath_id = sa.Column(
        sa.String(
            SHORT,
            collation=db_specific_properties[conf.CONFIGS['db_type']]['BINARY_COLLATION']),
        primary_key=True)
    pdb_id = sa.Column(sa.String(SHORT), nullable=False)
    pdb_chain = sa.Column(
        sa.String(
            SHORT,
            collation=db_specific_properties[conf.CONFIGS['db_type']]['BINARY_COLLATION']),
        nullable=False)
    pdb_domain_def = sa.Column(sa.String(MEDIUM), nullable=False)
    pdb_pdbfam_name = sa.Column(sa.String(LONG), nullable=False)
    pdb_pdbfam_idx = sa.Column(sa.Integer)
    domain_errors = sa.Column(sa.Text)


class DomainContact(Base):
    r"""Interactions between Profs domains in the PDB.

    Only interactions that were predicted to be biologically relevant by `NOXclass`_
    are included in this table.

    Columns:
      domain_contact_id
        A unique integer identifying each domain pair.

      cath_id_1
        Unique id identifying the first interacting domain in the :ref:`domain` table.

      cath_id_2
        Unique id identifying the second interacting domain in the :ref:`domain` table.

      min_interchain_distance
        The closest that any residue in domain one comes to any residue in domain two.

      contact_volume
        The volume covered by contacting residues.

      contact_surface_area
        The surface area of the contacting regions of the first and second domains.

      atom_count_1
        The number of atoms in the first domain.

      atom_count_2
        The number of atoms in the second domain.

      number_of_contact_residues_1
        The number of residues in the first domain that come within 5 \u212B of the second domain.

      number_of_contact_residues_2
        The number of residues in the second domain that come withing 5 \u212B of the first domain.

      contact_residues_1
        A list of all residues in the first domain that come within 5 \u212B of the second domain.
        The residue number corresponds to the position of the residue in the domain.

      contact_residues_2
        A list of all residues in the second domain that come within 5 \u212B of the first domain.
        The residue number corresponds to the position of the residue in the domain.

      crystal_packing
        The probability that the interaction is a crystallization artifacts, as defined by
        `NOXclass`_.

      domain_contact_errors
        List of errors that occurred when annotating this domain pair, or when using this domain
        as a template for making structural homology models.

    .. _NOXclass: http://noxclass.bioinf.mpi-inf.mpg.de/
    """

    __tablename__ = 'domain_contact'
    _indexes = [
        (['cath_id_1', 'cath_id_2'], {'unique': True}),
        (['cath_id_2', 'cath_id_1'], {'unique': True}),
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, [])

    domain_contact_id = sa.Column(sa.Integer, primary_key=True)
    cath_id_1 = sa.Column(
        None, sa.ForeignKey(Domain.cath_id, onupdate='cascade', ondelete='cascade'),
        nullable=False)
    cath_id_2 = sa.Column(
        None, sa.ForeignKey(Domain.cath_id, onupdate='cascade', ondelete='cascade'),
        nullable=False)
#    cath_id_2 = sa.Column(
#        sa.String(SHORT, collation=get_db_specific_param('BINARY_COLLATION')),
#        nullable=False)
    min_interchain_distance = sa.Column(sa.Float)
    contact_volume = sa.Column(sa.Float)
    contact_surface_area = sa.Column(sa.Float)
    atom_count_1 = sa.Column(sa.Integer)
    atom_count_2 = sa.Column(sa.Integer)
    number_of_contact_residues_1 = sa.Column(sa.Integer)
    number_of_contact_residues_2 = sa.Column(sa.Integer)
    contact_residues_1 = sa.Column(sa.Text)
    contact_residues_2 = sa.Column(sa.Text)
    crystal_packing = sa.Column(sa.Float)
    domain_contact_errors = sa.Column(sa.Text)

    # Relationships
    domain_1 = sa.orm.relationship(
        Domain, primaryjoin=cath_id_1 == Domain.cath_id, cascade='expunge', lazy='joined')
    # the second domain may be a ligand or a peptide, and so the foreign key constraint
    # does not work
    domain_2 = sa.orm.relationship(
        Domain, primaryjoin=cath_id_2 == Domain.cath_id, cascade='expunge', lazy='joined')


class UniprotSequence(Base):
    """Protein sequences from the Uniprot KB.

    Obtained by parsing `uniprot_sprot_fasta.gz`, `uniprot_trembl_fasta.gz`, and
    `homo_sapiens_variation.txt` files from the `Uniprot ftp site`_.

    Columns:
      db
        The database to which the protein sequence belongs. Possible values are `sp` for SwissProt
        and `tr` for TrEMBL.

      uniprot_id
        The uniprot id of the protein.

      uniprot_name
        The uniprot name of the protein.

      protein_name
        The protein name.

      organism_name
        Name of the organism in which this protein is found.

      gene_name
        Name of the gene that codes for this protein sequence.

      protein_existence
        Evidence for the existence of the protein:

        1. Experimental evidence at protein level
        2. Experimental evidence at transcript level
        3. Protein inferred from homology
        4. Protein predicted
        5. Protein uncertain

      sequence_version
        Version of the protein amino acid sequence.

      uniprot_sequence
        Amino acid sequence of the protein.

    .. _Uniprot ftp site:
       ftp://ftp.uniprot.org/pub/databases/uniprot/current_release/knowledgebase/

    """

    __tablename__ = 'uniprot_sequence'
    __table_args__ = get_table_args(__tablename__, [], [])

    db = sa.Column(sa.String(SHORT), nullable=False)
    uniprot_id = sa.Column(sa.String(SHORT), primary_key=True)
    uniprot_name = sa.Column(sa.String(SHORT), nullable=False)
    protein_name = sa.Column(sa.String(MEDIUM))
    organism_name = sa.Column(sa.String(MEDIUM), index=True)
    gene_name = sa.Column(sa.String(MEDIUM), index=True)
    protein_existence = sa.Column(sa.Integer)
    sequence_version = sa.Column(sa.Integer)
    uniprot_sequence = sa.Column(sa.Text, nullable=False)


class Provean(Base):
    """Description of the `Provean`_ supporting set calculated for a protein sequence.

    The construction of a supporting set is the most lengthy step in running Provean.
    Therefore, the supporting set is precalculated and stored for every protein sequence.

    Columns:
      uniprot_id
        The uniprot id of the protein.

      provean_supset_filename
        The filename of the Provean supporting set. The supporting set contains the ids
        and sequences of all proteins in the NCBI nr database that are used by Provean
        to construct a multiple sequence alignment for the given protein.

      provean_supset_length
        The number of sequences in Provean supporting set.

      provean_errors
        List of errors that occurred while the Provean supporting set was being calculated.

      provean_date_modified
        Date and time that this row was last modified.

    .. _provean: http://provean.jcvi.org/downloads.php
    """

    __tablename__ = 'provean'
    __table_args__ = get_table_args(__tablename__, [], [])

    uniprot_id = sa.Column(
        None, sa.ForeignKey(
            UniprotSequence.uniprot_id,
            onupdate='cascade', ondelete='cascade'),
        primary_key=True)
    provean_supset_filename = sa.Column(sa.String(MEDIUM))
    provean_supset_length = sa.Column(sa.Integer)
    provean_errors = sa.Column(sa.Text)
    provean_date_modified = sa.Column(
        sa.DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow,
        nullable=False)

    # Relationships
    uniprot_sequence = sa.orm.relationship(
        UniprotSequence, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('provean', uselist=False, cascade='expunge', lazy='joined'))


class UniprotDomain(Base):
    """Pfam domain definitions for proteins in the :ref:`uniprot_sequence` table.

    This table was obtained by downloading Pfam domain definitions for all known proteins
    from the `SIMAP`_ website, and mapping the protein sequence to uniprot using the MD5 hash
    of each sequence.

    Columns:
      uniprot_domain_id
        Unique id identifying each domain.

      uniprot_id
        The uniprot id of the protein containing the domain.

      pdbfam_name
        The Profs name of the domain. In most cases this will be equivalent
        to the Pfam name of the domain.

      pdbfam_idx
        The index of the Profs domain. ``pdbfam_idx`` ranges from 1 to the number of domains with
        the name ``pdbfam_name`` in the given protein. The ``(pdbfam_name, pdbfam_idx)`` tuple
        uniquely identifies each domain.

      pfam_clan
        The Pfam clan to which this Profs domain belongs.

      alignment_def
        Alignment domain definitions of the Profs domain. This field is obtained by removing gaps
        in the ``alignment_subdefs`` column.

      pfam_names
        Pfam names of all Pfam domains that were combined to create the given Profs domain.

      alignment_subdefs
        Comma-separated list of domain definitions for all Pfam domains that were merged to create
        the given Profs domain.

      path_to_data
        Location for storing homology models, mutation results, and all other data that
        are relevant to this domain. This path is prefixed by :term:`archive_dir`.

    .. _SIMAP: http://liferay.csb.univie.ac.at/portal/web/simap
    """

    __tablename__ = 'uniprot_domain'

    uniprot_domain_id = sa.Column(sa.Integer, nullable=False, primary_key=True, autoincrement=True)
    uniprot_id = sa.Column(
        None, sa.ForeignKey(
            UniprotSequence.uniprot_id,
            onupdate='cascade', ondelete='cascade'),
        nullable=False)
    pdbfam_name = sa.Column(sa.String(LONG), nullable=False)
    pdbfam_idx = sa.Column(sa.Integer, nullable=False)
    pfam_clan = sa.Column(sa.Text)
    alignment_def = sa.Column(sa.String(MEDIUM))
    pfam_names = sa.Column(sa.String(LONG))
    alignment_subdefs = sa.Column(sa.Text)
    path_to_data = sa.Column(sa.Text)

    if 'training' in conf.CONFIGS.get('db_schema', ''):
        # The database used for storing training data has an extra column `max_seq_identity`,
        # because we want to make homology models at different sequence identities.
        max_seq_identity = sa.Column(sa.Integer, index=True)
        _indexes = [
            (['uniprot_id', 'alignment_def', 'max_seq_identity'],
             {'unique': True, 'index_name': 'ix_uniprot_id_unique'}),
            (['pdbfam_name'], {'mysql_length': 255}),
            (['uniprot_id', 'uniprot_domain_id'],
             {'unique': True, 'index_name': 'ix_uniprot_id_uniprot_domain_id'}),
        ]
    else:
        _indexes = [
            (['pdbfam_name'], {'mysql_length': 255}),
            (['uniprot_id', 'uniprot_domain_id'],
             {'unique': True, 'index_name': 'ix_uniprot_id_uniprot_domain_id'}),
        ]
    __table_args__ = get_table_args(__tablename__, _indexes, [])

    # Relationships
    uniprot_sequence = sa.orm.relationship(
        UniprotSequence, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('uniprot_domain', cascade='expunge'))  # many to one


class UniprotDomainPair(Base):
    """Potentially-interacting pairs of domains for proteins that are known to interact.

    Columns:
      uniprot_domain_pair_id
        Unique id identifying each domain-domain interaction.

      uniprot_domain_id_1
        Unique id of the first domain.

      uniprot_domain_id_2
        Unique id of the second domain.

      rigids
        Phased out.

      domain_contact_ids
        List of unique ids identifying all domain-domain pairs in the PDB, where one domain
        belongs to the protein containing ``uniprot_domain_id_1`` and the other domain
        belongs to the protein containing ``uniprot_domain_id_2``. This was used as
        crystallographic evidence that the two proteins interact.

      path_to_data
        Location for storing homology models, mutation results, and all other data that is relevant
        to this domain pair. This path is prefixed by :term:`archive_dir`.
    """

    __tablename__ = 'uniprot_domain_pair'

    uniprot_domain_pair_id = sa.Column(sa.Integer, primary_key=True, autoincrement=True)
    uniprot_domain_id_1 = sa.Column(
        None, sa.ForeignKey(
            UniprotDomain.uniprot_domain_id,
            onupdate='cascade', ondelete='cascade'),
        nullable=False)
    uniprot_domain_id_2 = sa.Column(
        None, sa.ForeignKey(
            UniprotDomain.uniprot_domain_id,
            onupdate='cascade', ondelete='cascade'),
        nullable=False)
    rigids = sa.Column(sa.Text)  # Interaction references from iRefsa.Index
    domain_contact_ids = sa.Column(sa.Text)  # interaction references from the PDB
    path_to_data = sa.Column(sa.Text)
    # TODO: Move these columns higher up the next time creating a database
    uniprot_id_1 = sa.Column(sa.String(MEDIUM))
    uniprot_id_2 = sa.Column(sa.String(MEDIUM))

    if 'training' in conf.CONFIGS.get('db_schema', ''):
        # The database used for storing training data has an extra column `max_seq_identity`,
        # because we want to make homology models at different sequence identities.
        max_seq_identity = sa.Column(sa.Integer, index=True)
        _indexes = [
            (['uniprot_id_1', 'uniprot_id_2', 'max_seq_identity']),
            (['uniprot_id_2', 'uniprot_id_1', 'max_seq_identity']),
            (['uniprot_domain_id_1', 'uniprot_domain_id_2', 'max_seq_identity']),
            (['uniprot_domain_id_2', 'uniprot_domain_id_1', 'max_seq_identity']),
        ]
    else:
        _indexes = [
            (['uniprot_id_1', 'uniprot_id_2']),
            (['uniprot_id_2', 'uniprot_id_1']),
            (['uniprot_domain_id_1', 'uniprot_domain_id_2'], {'unique': True}),
            (['uniprot_domain_id_2', 'uniprot_domain_id_1'], {'unique': True}),
        ]
    __table_args__ = get_table_args(__tablename__, _indexes, [])

    # Relationships
    uniprot_domain_1 = sa.orm.relationship(
        UniprotDomain,
        primaryjoin=uniprot_domain_id_1 == UniprotDomain.uniprot_domain_id,
        cascade='expunge', lazy='joined')  # many to one
    uniprot_domain_2 = sa.orm.relationship(
        UniprotDomain,
        primaryjoin=uniprot_domain_id_2 == UniprotDomain.uniprot_domain_id,
        cascade='expunge', lazy='joined')  # many to one


class UniprotDomainTemplate(Base):
    r"""Structural templates for domains in the :ref:`uniprot_domain` table.

    Lists PDB crystal structures that will be used for making homology models.

    Columns:
      uniprot_domain_id
        An integer which uniquely identifies each uniprot domain in the
        :ref:`uniprot_domain` table.

      template_errors
        List of errors that occurred during the process for finding the template.

      cath_id
        The unique id identifying the structural template of the domain.

      domain_start
        The Uniprot position of the first amino acid of the Profs domain.

      domain_end
        The Uniprot position of the last amino acid of the Profs domain.

      domain_def
        Profs domain definitions for domains with structural templates. Domain definitions in this
        column are different from domain definitions in the ``alignment_def`` column of the
        :ref:`uniprot_domain` table in that they have been expanded to match domain boundaries
        of the Profs structural template, identified by the ``cath_id``.

      alignment_identity
        Percent identity of the domain to its structural template.

      alignment_coverage
        Percent coverage of the domain to its structural template.

      alignment_score
        A score obtained by combining ``alignment_identity`` (:math:`SeqId`) and
        ``alignment_coverage`` (:math:`Cov`) using the following equation, as described by
        `Mosca et al.`_:

        .. math::
           :label: score_function

           Score = 0.95 \\cdot \\frac{SeqId}{100} \\cdot \\frac{Cov}{100} + \
                   0.05 \\cdot \\frac{Cov}{100}

      t_date_modified
        The date and time when this row was last modified.

    .. _Mosca et al.: http://doi.org/10.1038/nmeth.2289
    """

    __tablename__ = 'uniprot_domain_template'
    __table_args__ = get_table_args(__tablename__, [], [])

    uniprot_domain_id = sa.Column(
        None, sa.ForeignKey(
            UniprotDomain.uniprot_domain_id,
            onupdate='cascade', ondelete='cascade'),
        nullable=False, primary_key=True)
    template_errors = sa.Column(sa.Text)
    cath_id = sa.Column(
        None, sa.ForeignKey(
            Domain.cath_id,
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False)
    domain_start = sa.Column(sa.Integer, index=True)
    domain_end = sa.Column(sa.Integer, index=True)
    domain_def = sa.Column(sa.String(MEDIUM))
    alignment_identity = sa.Column(sa.Float, sa.CheckConstraint('alignment_identity <= 1'))
    alignment_coverage = sa.Column(sa.Float, sa.CheckConstraint('alignment_coverage <= 1'))
    alignment_score = sa.Column(sa.Float, sa.CheckConstraint('alignment_score <= 1'))
    t_date_modified = sa.Column(
        sa.DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow,
        nullable=False)

    # Relationships
    uniprot_domain = sa.orm.relationship(
        UniprotDomain, uselist=False, cascade='expunge', lazy='joined', innerjoin=True,
        backref=sa.orm.backref(
            'template', uselist=False, cascade='expunge', innerjoin=True, lazy='joined'))
    domain = sa.orm.relationship(
        Domain, uselist=False, cascade='expunge', lazy='joined', innerjoin=True,
        backref=sa.orm.backref('uniprot_domain', cascade='expunge'))


class UniprotDomainModel(Base):
    """Homology models for templates in the :ref:`uniprot_domain_template` table.

    Columns:
      uniprot_domain_id
        An integer which uniquely identifies each uniprot domain in the :ref:`uniprot_domain`
        table.

      model_errors
        List of errors that occurred when making the homology model.

      alignment_filename
        The name of the alignment that was given to Modeller when making the homology model.

      model_filename
        The name of the homology model that was produced by Modeller.

      chain
        The chain that contains the domain in question in the homology (this is now set to 'A'
        in all models).

      norm_dope
        Normalized DOPE score of the model (lower is better).

      sasa_score
        Comma-separated list of the percent solvent-accessible surface area for each residue.

      m_date_modified
        The date and time when this row was last modified.

      model_domain_def
        Domain definitions for the region of the domain that is covered by the structural template.

        In most cases, this field is identical to the ``domain_def`` field in the
        :ref:`uniprot_domain_template` table. However, it sometimes happens that the best
        Profs structural template only covers a fraction of the Pfam domain. In that case, the
        ``alignment_def`` column in the :ref:`uniprot_domain` table, and the ``domain_def`` column
        in the :ref:`uniprot_domain_template` table, will contain the original Pfam domain
        definitions, and the ``model_domain_def`` column will contain domain definitions for only
        the region that is covered by the structural template.
    """

    __tablename__ = 'uniprot_domain_model'
    __table_args__ = get_table_args(__tablename__, [], [])

    uniprot_domain_id = sa.Column(
        None, sa.ForeignKey(
            UniprotDomainTemplate.uniprot_domain_id,
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False, primary_key=True)
    model_errors = sa.Column(sa.Text)
    alignment_filename = sa.Column(sa.String(MEDIUM))
    model_filename = sa.Column(sa.String(MEDIUM))
    chain = sa.Column(sa.String(SHORT))
    norm_dope = sa.Column(sa.Float)
    sasa_score = sa.Column(sa.Text)
    m_date_modified = sa.Column(
        sa.DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow,
        nullable=False)
    model_domain_def = sa.Column(sa.String(MEDIUM))

    # Relationships
    # one to one
    template = sa.orm.relationship(
        UniprotDomainTemplate, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('model', uselist=False, cascade='expunge', lazy='joined'))


class UniprotDomainMutation(Base):
    """Characterization of mutations introduced into structures.

    References the :ref:`uniprot_domain_model` table.

    Columns:
      uniprot_id
        Uniprot ID of the protein that was mutated.

      uniprot_domain_id
        Unique id which identifies the Profs domain that was mutated in the :ref:`uniprot_domain`
        table.

      mutation
        Mutation that was introduced into the protein, in Uniprot coordinates.

      mutation_errors
        List of errors that occured while evaluating the mutation.

      model_filename_wt
        The name of the file which contains the homology model of the domain after the model was
        relaxed with FoldX but before the mutation was introduced.

      model_filename_mut
        The name of the file which contains the homology model of the domain after the model was
        relaxed with FoldX and after the mutation was introduced.

      chain_modeller
        The chain which contains the domain that was mutated in the ``model_filename_wt`` and the
        ``model_filename_mut`` structures.

      mutation_modeller
        The mutation that was introduced into the protein, in PDB RESNUM coordinates.
        This identifies the mutated residue in the ``model_filename_wt`` and the
        ``model_filename_mut`` structures.

      stability_energy_wt
        Comma-separated list of scores returned by FoldX for the wildtype protein.
        The comma-separated list can be converted into a DataFrame with each column clearly
        labelled using the :func:`elaspic.predictor.format_mutation_features`.
        The FoldX energy terms are:

            - dg
            - backbone_hbond
            - sidechain_hbond
            - van_der_waals
            - electrostatics
            - solvation_polar
            - solvation_hydrophobic
            - van_der_waals_clashes
            - entropy_sidechain
            - entropy_mainchain
            - sloop_entropy
            - mloop_entropy
            - cis_bond
            - torsional_clash
            - backbone_clash
            - helix_dipole
            - water_bridge
            - disulfide
            - electrostatic_kon
            - partial_covalent_bonds
            - energy_ionisation
            - entropy_complex
            - number_of_residues

      stability_energy_mut
        Comma-separated list of scores returned by FoldX for the mutant protein.
        FoldX energy terms are the same as in `stability_energy_wt`, but for the mutated amino
        acid rather than the wildtype.

      physchem_wt
        Physicochemical properties describing the interaction of the wildtype residue with residues
        on the opposite chain. The terms are:

            - number of atoms in interacting residues that have the same charge.
            - number of atoms in interacting residues that have an opposite charge.
            - number of hydrogen bonds (very rough calculation).
            - number of carbons in interacting residues within 4 A of the mutated residue
              (rough measure of the van der Waals force).

      physchem_wt_ownchain
        Physicochemical properties describing the interaction of the wildtype residue with residues
        on the same chain. The terms are the same as in `physchem_wt`.

      physchem_mut
        Physicochemical properties describing the interaction of the mutant residue with residues
        on the opposite chain. The terms are the same as in `physchem_wt`.

      physchem_mut_ownchain
        Physicochemical properties describing the interaction of the mutant residue with residues
        on the same chain. The terms are the same as in `physchem_wt`.

      matrix_score
        Score assigned to the wt -> mut transition by the BLOSUM substitution matrix.

      secondary_structure_wt
        Secondary structure of the wildtype residue predicted by `stride`_.

      solvent_accessibility_wt
        Percent solvent accessible surface area of the wildtype residue, predicted by `msms`_.

      secondary_structure_mut
        Secondary structure of the mutated residue predicted by `stride`_.

      solvent_accessibility_mut
        Percent solvent accessible surface area of the mutated residue, predicted by `msms`_.

      provean_score
        Score produced by `Provean`_ for this mutation.

      ddg
        Change in the Gibbs free energy of folding that our classifier predicts for this mutation.

      mut_date_modified
        Date and time that this row was last modified.

    .. _stride: http://webclu.bio.wzw.tum.de/stride/
    .. _msms: http://mgltools.scripps.edu/
    """

    __tablename__ = 'uniprot_domain_mutation'
    _indexes = [
        ['uniprot_id', 'mutation'],
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, [])

    uniprot_id = sa.Column(
        None, sa.ForeignKey(
            UniprotSequence.uniprot_id,
            onupdate='cascade', ondelete='cascade'),
        nullable=False, primary_key=True)
    uniprot_domain_id = sa.Column(
        None, sa.ForeignKey(
            UniprotDomainModel.uniprot_domain_id,
            onupdate='cascade', ondelete='cascade'),
        nullable=False, primary_key=True, index=True)
    mutation = sa.Column(sa.String(SHORT), index=True, nullable=False, primary_key=True)
    mutation_errors = sa.Column(sa.Text)
    model_filename_wt = sa.Column(sa.String(MEDIUM))
    model_filename_mut = sa.Column(sa.String(MEDIUM))
    chain_modeller = sa.Column(sa.String(SHORT))
    mutation_modeller = sa.Column(sa.String(SHORT))
    stability_energy_wt = sa.Column(sa.Text)
    stability_energy_mut = sa.Column(sa.Text)
    physchem_wt = sa.Column(sa.Text)
    physchem_wt_ownchain = sa.Column(sa.Text)
    physchem_mut = sa.Column(sa.Text)
    physchem_mut_ownchain = sa.Column(sa.Text)
    matrix_score = sa.Column(sa.Float)
    secondary_structure_wt = sa.Column(sa.Text)
    solvent_accessibility_wt = sa.Column(sa.Float)
    secondary_structure_mut = sa.Column(sa.Text)
    solvent_accessibility_mut = sa.Column(sa.Float)
    provean_score = sa.Column(sa.Float)
    ddg = sa.Column(sa.Float, index=True)
    mut_date_modified = sa.Column(
        sa.DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow,
        nullable=False)

    # Relationships
    model = sa.orm.relationship(
        UniprotDomainModel, cascade='expunge', uselist=False, lazy='joined',
        backref=sa.orm.backref('mutations', cascade='expunge'))  # many to one


class UniprotDomainPairTemplate(Base):
    r"""Structural templates for pairs of domains in the :ref:`uniprot_domain_pair` table.

    Columns:
      uniprot_domain_pair_id
        Unique id identifying each domain-domain interaction.

      domain_contact_id
          Unique id of the domain pair in the :ref:`domain_contact` table
          that was used as a template for the modelled domain pair.

      cath_id_1
        Unique id of the structural template for the first domain.

      cath_id_2
        Unique id of the structural template for the second domain.

      identical_1
        Fraction of residues in the Blast alignment of the first domain to its template that are
        *identical*.

      conserved_1
        Fraction of residues in the Blast alignment of the first domain to its template
        that are *conserved*.

      coverage_1
        Fraction of the first domain that is covered by the blast alignment.

      score_1
        Score obtained by multiplying ``identical_1`` by ``coverage_1``.

      identical_if_1
        Fraction of interface residues [#f1]_ that are *identical* in the Blast alignment
        of the first domain.

      conserved_if_1
        Fraction of interface residues [#f1]_ that are *conserved* in the Blast alignment
        of the first domain.

      coverage_if_1
        Fraction of interface residues [#f1]_ that are *covered* by the Blast alignment
        of the first domain.

      score_if_1
        Score obtained by combining ``identical_if_1`` and ``coverage_if_1`` using
        :eq:`score_function`.

      identical_2
        Fraction of residues in the Blast alignment of the second domain to its template that are
        *identical*.

      conserved_2
        Fraction of residues in the Blast alignment of the second domain to its template
        that are *conserved*.

      coverage_2
        Fraction of the second domain that is covered by the blast alignment.

      score_2
        Score obtained by multiplying ``identical_2`` by ``coverage_2``.

      identical_if_2
        Fraction of interface residues [#f1]_ that are *identical* in the Blast alignment
        of the second domain.

      conserved_if_2
        Fraction of interface residues [#f1]_ that are *conserved* in the Blast alignment
        of the second domain.

      coverage_if_2
        Fraction of interface residues [#f1]_ that are *covered* by the Blast alignment
        of the second domain.

      score_if_2
        Score obtained by combining ``identical_if_2`` and ``coverage_if_2``
        using :eq:`score_function`.

      score_total
        The product of ``score_1`` and ``score_2``.

      score_if_total
        The product of ``score_if_1`` and ``score_if_2``.

      score_overall
        The product of ``score_total`` and ``score_if_total``. This is the score that was used to
        select the best Profs domain pair to be used as a template.

      t_date_modified
        The date and time when this row was last updated.

      template_errors
        List of errors that occured while looking for the structural template.


    .. [#f1] Interface residues are defined as residues that are within 5 \u212B
              of the partner domain.
    """

    __tablename__ = 'uniprot_domain_pair_template'
    _indexes = [
        ['cath_id_1', 'cath_id_2'],
        ['cath_id_2', 'cath_id_1'],
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, [])

    uniprot_domain_pair_id = sa.Column(
        None, sa.ForeignKey(
            UniprotDomainPair.uniprot_domain_pair_id,
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False, primary_key=True)
    # domain_contact_id = sa.Column(
    #     None, sa.ForeignKey(
    #         DomainContact.domain_contact_id,
    #         onupdate='cascade', ondelete='cascade'),
    #     index=True, nullable=False)
    # Had to comment out code above because in case of training
    # we may not have correct `domain_contact_id`
    domain_contact_id = sa.Column(sa.Integer, nullable=True)
    cath_id_1 = sa.Column(
        None, sa.ForeignKey(
            Domain.cath_id,
            onupdate='cascade', ondelete='cascade'),
        nullable=False)
    cath_id_2 = sa.Column(
        None, sa.ForeignKey(
            Domain.cath_id,
            onupdate='cascade', ondelete='cascade'),
        nullable=False)

    identical_1 = sa.Column(sa.Float, sa.CheckConstraint('identical_1 <= 1'))
    conserved_1 = sa.Column(sa.Float, sa.CheckConstraint('conserved_1 <= 1'))
    coverage_1 = sa.Column(sa.Float, sa.CheckConstraint('coverage_1 <= 1'))
    score_1 = sa.Column(sa.Float, sa.CheckConstraint('score_1 <= 1'))

    identical_if_1 = sa.Column(sa.Float, sa.CheckConstraint('identical_if_1 <= 1'))
    conserved_if_1 = sa.Column(sa.Float, sa.CheckConstraint('conserved_if_1 <= 1'))
    coverage_if_1 = sa.Column(sa.Float, sa.CheckConstraint('coverage_if_1 <= 1'))
    score_if_1 = sa.Column(sa.Float, sa.CheckConstraint('score_if_1 <= 1'))

    identical_2 = sa.Column(sa.Float, sa.CheckConstraint('identical_2 <= 1'))
    conserved_2 = sa.Column(sa.Float, sa.CheckConstraint('conserved_2 <= 1'))
    coverage_2 = sa.Column(sa.Float, sa.CheckConstraint('coverage_2 <= 1'))
    score_2 = sa.Column(sa.Float, sa.CheckConstraint('score_2 <= 1'))

    identical_if_2 = sa.Column(sa.Float, sa.CheckConstraint('identical_if_2 <= 1'))
    conserved_if_2 = sa.Column(sa.Float, sa.CheckConstraint('conserved_if_2 <= 1'))
    coverage_if_2 = sa.Column(sa.Float, sa.CheckConstraint('coverage_if_2 <= 1'))
    score_if_2 = sa.Column(sa.Float, sa.CheckConstraint('score_if_2 <= 1'))

    score_total = sa.Column(sa.Float, sa.CheckConstraint('score_total <= 1'))
    score_if_total = sa.Column(sa.Float, sa.CheckConstraint('score_if_total <= 1'))
    score_overall = sa.Column(sa.Float, sa.CheckConstraint('score_overall <= 1'))

    t_date_modified = sa.Column(
        sa.DateTime, default=datetime.datetime.utcnow,
        onupdate=datetime.datetime.utcnow, nullable=False)
    template_errors = sa.Column(sa.Text)

    # Relationships
    # one to one
    domain_pair = sa.orm.relationship(
        UniprotDomainPair, uselist=False, cascade='expunge', lazy='joined', innerjoin=True,
        backref=sa.orm.backref(
            'template', uselist=False, cascade='expunge', lazy='joined', innerjoin=True))
    # domain_contact = sa.orm.relationship(
    #     DomainContact, uselist=False, cascade='expunge', lazy='joined', innerjoin=True,
    #     backref=sa.orm.backref('uniprot', cascade='expunge'))  # one to one
    domain_1 = sa.orm.relationship(
        Domain, uselist=False, cascade='expunge', lazy='joined',
        primaryjoin=(cath_id_1 == Domain.cath_id))  # many to one
    domain_2 = sa.orm.relationship(
        Domain, uselist=False, cascade='expunge', lazy='joined',
        primaryjoin=(cath_id_2 == Domain.cath_id))  # many to one


class UniprotDomainPairModel(Base):
    r"""Structural models of interactions between pairs of domains.

    References the :ref:`uniprot_domain_pair` table.

    Columns:
      uniprot_domain_pair_id
        Unique id identifying each domain-domain interaction.

      model_errors
        List of errors that occured while making the homology model.

      alignment_filename_1
        Name of the file containing the alignment of the first domain
        with its structural template.

      alignment_filename_2
        Name of the file containing the alignment of the second domain
        with its structural template.

      model_filename
        Name of the file containing the homology model of the domain-domain interaction
        created by Modeller.

      chain_1
        Chain containing the first domain in the model specified by ``model_filename``.

      chain_2
        Chain containing the second domain in the model specified by ``model_filename``.

      norm_dope
        The normalized DOPE score of the model.

      interface_area_hydrophobic
        Hydrophobic surface area of the interface, calculated using `POPS`_.

      interface_area_hydrophilic
        Hydrophilic surface area of the interface, calculated using `POPS`_.

      interface_area_total
        Total surface area of the interface, calculated using `POPS`_.

      interface_dg
        Gibbs free energy of binding for this domain-domain interaction, predicted using `FoldX`_.
        Not implemented yet!

      interacting_aa_1
        List of amino acid positions in the first domain that are within 5 \u212B
        of the second domain. Positions are specified using uniprot coordinates.

      interacting_aa_2
        List of amino acids in the second domain that are within 5 \u212B of the first domain.
        Position are specified using uniprot coordinates.

      m_date_modified
        Date and time that this row was last modified.

      model_domain_def_1
        Domain boundaries of the first domain that are covered by the Profs structural template.

      model_domain_def_2
        Domain boundaries of the second domain that are covered by the Profs structural template.

    .. _POPS: http://mathbio.nimr.mrc.ac.uk/wiki/POPS
    .. _FoldX: http://foldx.crg.es/
    """

    __tablename__ = 'uniprot_domain_pair_model'
    __table_args__ = get_table_args(__tablename__, [], [])

    uniprot_domain_pair_id = sa.Column(
        None, sa.ForeignKey(
            UniprotDomainPairTemplate.uniprot_domain_pair_id,
            onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False, primary_key=True)
    model_errors = sa.Column(sa.Text)
    alignment_filename_1 = sa.Column(sa.String(MEDIUM))
    alignment_filename_2 = sa.Column(sa.String(MEDIUM))
    model_filename = sa.Column(sa.String(MEDIUM))
    chain_1 = sa.Column(sa.String(SHORT))
    chain_2 = sa.Column(sa.String(SHORT))
    norm_dope = sa.Column(sa.Float)
    interface_area_hydrophobic = sa.Column(sa.Float)
    interface_area_hydrophilic = sa.Column(sa.Float)
    interface_area_total = sa.Column(sa.Float)
    interface_dg = sa.Column(sa.Float)
    interacting_aa_1 = sa.Column(sa.Text)
    interacting_aa_2 = sa.Column(sa.Text)
    m_date_modified = sa.Column(
        sa.DateTime, default=datetime.datetime.utcnow, onupdate=datetime.datetime.utcnow,
        nullable=False)
    model_domain_def_1 = sa.Column(sa.String(MEDIUM))
    model_domain_def_2 = sa.Column(sa.String(MEDIUM))

    # Relationships
    # one to one
    template = sa.orm.relationship(
        UniprotDomainPairTemplate, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('model', uselist=False, cascade='expunge', lazy='joined'))


class UniprotDomainPairMutation(Base):
    r"""Characterization of interface mutations introduced into structures.

    References the :ref:`uniprot_domain_pair_model` table.

    Columns:
      uniprot_id
        Uniprot ID of the protein that is being mutated.

      uniprot_domain_pair_id
        Unique id identifying each domain-domain interaction.

      mutation
        Mutation for which the :math:`\Delta \Delta G` score is being predicted, specified in
        Uniprot coordinates.

      mutation_errors
        List of errors obtained when evaluating the impact of the mutation.

      model_filename_wt
        Filename of the homology model relaxed by FoldX but containing the wildtype residue.

      model_filename_mut
        Filename of the homology model relaxed by FoldX and containing the mutated residue.

      chain_modeller
        Chain containing the domain that was mutated, in homology models specified by
        ``model_filename_wt`` and ``model_filename_mut``.

      mutation_modeller
        Mutation for which the :math:`\Delta \Delta G` score is being predicted,
        specified in PDB RESNUM coordinates.

      analyse_complex_energy_wt
        Comma-separated list of FoldX scores describing the effect of the wildtype residue on
        the stability of the protein domain.

      stability_energy_wt
        Comma-separated list of FoldX scores describing the effect of the wildtype residue on
        protein-protein interaction interface.

      analyse_complex_energy_mut
        Comma-separated list of FoldX scores describing the effect of the mutated residue on
        the stability of the protein domain.

      stability_energy_mut
        Comma-separated list of FoldX scores describing the effect of the mutated residue on
        protein-protein interaction interface.

      physchem_wt
        Comma-separated list of physicochemical properties describing the interaction between
        the wildtype residue and other residues on the opposite chain.

      physchem_wt_ownchain
        Comma-separated list of physicochemical properties describing the interaction between
        the wildtype residue and other residues on the same chain.

      physchem_mut
        Comma-separated list of physicochemical properties describing the interaction between
        the mutated residue and other residues on the opposite chain.

      physchem_mut_ownchain
        Comma-separated list of physicochemical properties describing the interaction between
        the mutated residue and other residues on the same chain.

      matrix_score
        Score assigned to the wt -> mut transition by the BLOSUM substitution matrix.

      secondary_structure_wt
        Secondary structure of the wildtype residue, predicted by `stride`_.

      solvent_accessibility_wt
        Percent solvent accessible surface area of the wildtype residue, predicted by `msms`_.

      secondary_structure_mut
        Secondary structure of the mutated residue, predicted by `stride`_.

      solvent_accessibility_mut
        Percent solvent accessible surface area of the mutated residue, predicted by `msms`_.

      contact_distance_wt
        Shortest distance between the wildtype residue and a residue on the opposite chain.

      contact_distance_mut
        Shortest distance between the mutated reside and a residue on the opposite chain.

      provean_score
        `Provean`_ score for this mutation.

      ddg
        Predicted change in Gibbs free energy of binding caused by this mutation.

      mut_date_modified
        Date and time when this row was last modified.
    """

    __tablename__ = 'uniprot_domain_pair_mutation'
    _indexes = [
        ['uniprot_id', 'mutation'],
    ]
    __table_args__ = get_table_args(__tablename__, _indexes, [])

    uniprot_id = sa.Column(None, sa.ForeignKey(
        UniprotSequence.uniprot_id, onupdate='cascade', ondelete='cascade'),
        nullable=False, primary_key=True)
    uniprot_domain_pair_id = sa.Column(None, sa.ForeignKey(
        UniprotDomainPairModel.uniprot_domain_pair_id, onupdate='cascade', ondelete='cascade'),
        index=True, nullable=False, primary_key=True)
    mutation = sa.Column(sa.String(SHORT), nullable=False, primary_key=True)
    mutation_errors = sa.Column(sa.Text)
    model_filename_wt = sa.Column(sa.String(MEDIUM))
    model_filename_mut = sa.Column(sa.String(MEDIUM))
    chain_modeller = sa.Column(sa.String(SHORT))
    mutation_modeller = sa.Column(sa.String(SHORT))
    analyse_complex_energy_wt = sa.Column(sa.Text)
    stability_energy_wt = sa.Column(sa.Text)
    analyse_complex_energy_mut = sa.Column(sa.Text)
    stability_energy_mut = sa.Column(sa.Text)
    physchem_wt = sa.Column(sa.Text)
    physchem_wt_ownchain = sa.Column(sa.Text)
    physchem_mut = sa.Column(sa.Text)
    physchem_mut_ownchain = sa.Column(sa.Text)
    matrix_score = sa.Column(sa.Float)
    secondary_structure_wt = sa.Column(sa.Text)
    solvent_accessibility_wt = sa.Column(sa.Float)
    secondary_structure_mut = sa.Column(sa.Text)
    solvent_accessibility_mut = sa.Column(sa.Float)
    contact_distance_wt = sa.Column(sa.Float)
    contact_distance_mut = sa.Column(sa.Float)
    provean_score = sa.Column(sa.Float)
    ddg = sa.Column(sa.Float, index=False)
    mut_date_modified = sa.Column(sa.DateTime, default=datetime.datetime.utcnow,
                                  onupdate=datetime.datetime.utcnow, nullable=False)
    # Relationships
    model = sa.orm.relationship(
        UniprotDomainPairModel, uselist=False, cascade='expunge', lazy='joined',
        backref=sa.orm.backref('mutations', cascade='expunge'))  # many to one
