CREATE OR REPLACE VIEW core_model AS
SELECT
-- domain
ud.uniprot_id protein_id,
ud.uniprot_domain_id domain_id,
ud.pdbfam_idx domain_idx,
ud.pfam_clan,
ud.pdbfam_name,
ud.alignment_def,
ud.path_to_data,

-- template
udt.alignment_score,
udt.alignment_coverage,
udt.template_errors,
udt.domain_def,
udt.cath_id,
udt.alignment_identity,

-- model
udm.model_errors,
udm.norm_dope,
udm.model_filename,
udm.alignment_filename,
udm.chain,
udm.model_domain_def

FROM elaspic.uniprot_domain ud
JOIN elaspic.uniprot_domain_template udt USING (uniprot_domain_id)
LEFT JOIN elaspic.uniprot_domain_model udm USING (uniprot_domain_id);
