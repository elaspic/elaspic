# -*- coding: utf-8 -*-
import os
import numpy as np
import pandas as pd

import class_sql as sql

os.chdir('/home/kimlab1/strokach/working/pipeline/code_supl/cosmic')
cosmic_df = pd.read_csv('CosmicV67-Recep.txt', sep='\t')
cosmic_driver_df = cosmic_df[cosmic_df['mType'] == 'Driver']
uniprots_to_precalculate = set(cosmic_driver_df['uniprot'])

db = sql.MyDatabase(path_to_archive='/home/kimlab1/database_data/elaspic/human')


###############################################################################
import json

with open('/home/kimlab1/database_data/elaspic/human/Q9Y/6K/Q9Y6K1/PWWP*291-374/template.json', 'r') as fh: template_data = json.load(fh)

template3 = sql.UniprotDomainTemplate(**template_data)
#template = sql.UniprotDomainTemplate()
template3.uniprot_domain_id = 1000
#template.template_errors = 'testing2'
db.add_uniprot_template(template3)

db.close()

###############################################################################

mutations = cosmic_driver_df[['uniprot', 'location']].values
num_of_mutations = mutations.shape[0]
sasa_score = np.empty((num_of_mutations,1), dtype=float)
errors_log = np.empty((num_of_mutations,1),dtype=object)
#uniprot_definitions = db.get_uniprot_domain('P04637', copy_data=False)
for i in range(num_of_mutations):
    uniprot_id = mutations[i,0]
    mutation_position = mutations[i,1]
    uniprot_definitions = db.get_uniprot_domain(uniprot_id, copy_data=False)
    sasa = None
    log = []
    for uniprot_domain, uniprot_template, uniprot_model in uniprot_definitions:
        if not uniprot_template \
        or not uniprot_template.domain_def \
        or not uniprot_model \
        or (not uniprot_model.model_errors and not uniprot_model.model_filename ):
            domain_start, domain_end = sql.decode_domain(uniprot_domain.alignment_def)
            if mutation_position < domain_start or mutation_position > domain_end:
                log.append(uniprot_domain.pfam_name + ': mutation outside of pfam domain')
            else:
                log.append(uniprot_domain.pfam_name + ': mutation falls inside pfam domain but no structural template was found')
        else:
            domain_start, domain_end = sql.decode_domain(uniprot_template.domain_def)
            if mutation_position < domain_start or mutation_position > domain_end:
                log.append(uniprot_domain.pfam_name + ': mutation outside of pdbfam domain')
            else:
                mutation_idx = domain_start - mutation_position - 1
                if uniprot_model.sasa_score is None:
                    log.append(uniprot_domain.pfam_name + ': mutation inside pdbfam domain but sasa score could not be calculated')
                else:
                    sasa = [score for score in uniprot_model.sasa_score.split(',')][mutation_idx]
                    log = [uniprot_domain.pfam_name + ': affected by mutation']
                    break
    sasa_score[i] = sasa
    errors_log[i] = ','.join(log)
    
cosmic_driver_df['sasa_score'] = sasa_score
cosmic_driver_df['errors_log'] = errors_log
print cosmic_driver_df
cosmic_driver_df.to_csv('/home/kimlab1/strokach/working/pipeline/results/cosmic/driver_data.tsv', sep='\t')
