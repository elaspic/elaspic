CREATE OR REPLACE VIEW protein AS
SELECT
uniprot_id protein_id,
uniprot_name protein_name,
protein_name description,
organism_name,
-- gene_name,
-- protein_existence isoforms,
-- sequence_version,
-- db,
uniprot_sequence sequence,
provean_supset_filename provean_supset_file,
provean_supset_length
FROM uniprot_kb.uniprot_sequence us
LEFT JOIN elaspic.provean USING (uniprot_id);
