CREATE OR REPLACE VIEW elaspic_interface_model AS
SELECT
-- interface
udp.uniprot_domain_pair_id interface_id,
udp.uniprot_id_1 protein_id_1,
udp.uniprot_domain_id_1 domain_id_1,
udp.uniprot_id_2 protein_id_2,
udp.uniprot_domain_id_2 domain_id_2,

udp.path_to_data,

-- template
udpt.score_1 alignment_score_1,
udpt.score_2 alignment_score_2,
udpt.coverage_1 alignment_coverage_1,
udpt.coverage_2 alignment_coverage_2,

udpt.cath_id_1,
udpt.cath_id_2,
udpt.identical_1 alignment_identity_1,
udpt.identical_2 alignment_identity_2,
udpt.template_errors,

-- model
udpm.model_errors,
udpm.norm_dope,
udpm.model_filename,
udpm.alignment_filename_1,
udpm.alignment_filename_2,
udpm.interacting_aa_1,
udpm.interacting_aa_2,
udpm.chain_1,
udpm.chain_2,
udpm.interface_area_hydrophobic,
udpm.interface_area_hydrophilic,
udpm.interface_area_total,
udpm.model_domain_def_1,
udpm.model_domain_def_2

FROM elaspic.uniprot_domain_pair udp
JOIN elaspic.uniprot_domain_pair_template udpt USING (uniprot_domain_pair_id)
LEFT JOIN elaspic.uniprot_domain_pair_model udpm USING (uniprot_domain_pair_id);
