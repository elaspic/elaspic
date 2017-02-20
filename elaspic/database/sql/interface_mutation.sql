CREATE OR REPLACE VIEW elaspic_interface_mutation AS
SELECT
udp.uniprot_domain_pair_id interface_id,
# udp.uniprot_id_1 protein_id_1,
# udp.uniprot_domain_id_1 domain_id_1,
# udp.uniprot_id_2 protein_id_2,
# udp.uniprot_domain_id_2 domain_id_2,
udpmut.uniprot_id protein_id,
udpmut.mutation,
IF(udpmut.chain_modeller = 'A', 0, 1) chain_idx,

-- mutation
udpmut.model_filename_wt,
udpmut.model_filename_mut,
udpmut.mut_date_modified,
udpmut.mutation_errors,
udpmut.chain_modeller,
udpmut.mutation_modeller,
udpmut.stability_energy_wt,
udpmut.stability_energy_mut,
udpmut.analyse_complex_energy_wt,
udpmut.analyse_complex_energy_mut,
udpmut.physchem_wt,
udpmut.physchem_wt_ownchain,
udpmut.physchem_mut,
udpmut.physchem_mut_ownchain,
udpmut.secondary_structure_wt,
udpmut.secondary_structure_mut,
udpmut.solvent_accessibility_wt,
udpmut.solvent_accessibility_mut,
udpmut.contact_distance_wt,
udpmut.contact_distance_mut,
udpmut.matrix_score,
udpmut.provean_score,
udpmut.ddg

FROM elaspic.uniprot_domain_pair udp
JOIN elaspic.uniprot_domain_pair_mutation udpmut USING (uniprot_domain_pair_id);
