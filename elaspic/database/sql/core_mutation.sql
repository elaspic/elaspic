CREATE OR REPLACE VIEW elaspic_core_mutation AS
SELECT
ud.uniprot_domain_id domain_id,
ud.uniprot_id protein_id,
ud.pdbfam_idx domain_idx,
udmut.mutation,

-- mutation
udmut.model_filename_wt,
udmut.model_filename_mut,
udmut.mut_date_modified,
udmut.mutation_errors,
udmut.chain_modeller,
udmut.mutation_modeller,
udmut.stability_energy_wt,
udmut.stability_energy_mut,
udmut.physchem_wt,
udmut.physchem_wt_ownchain,
udmut.physchem_mut,
udmut.physchem_mut_ownchain,
udmut.secondary_structure_wt,
udmut.secondary_structure_mut,
udmut.solvent_accessibility_wt,
udmut.solvent_accessibility_mut,
udmut.matrix_score,
udmut.provean_score,
udmut.ddG

FROM elaspic.uniprot_domain ud
JOIN elaspic.uniprot_domain_mutation udmut USING (uniprot_id, uniprot_domain_id);
