# -*- coding: utf-8 -*-

#import __main__
#__main__.pymol_argv = [ 'pymol', '-c'] # Quiet and no GUI (-qc)

import sys
import time
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import Tkinter
sys.modules['tkinter'] = Tkinter

import json
import pymol
import class_sql as sql

from Bio import SeqIO

import class_analyze_structure

#pymol.finish_launching()

###############################################################################

A_DICT = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU', \
          'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', \
          'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP', \
          'Y':'TYR', 'V':'VAL', 'U':'SEC', 'O':'PYL', \
          'B':'ASX', 'Z':'GLX', 'J':'XLE', 'X':'XAA', '*':'TER'}

AAA_DICT = dict([(value,key) for key,value in A_DICT.items()])


os.chdir('/home/kimlab1/strokach/working/pipeline/code_supl/cosmic')
cosmic_df = pd.read_csv('CosmicV67-Recep.txt', sep='\t')
cosmic_driver_df = cosmic_df[cosmic_df['mType'] == 'Driver']
uniprots_to_precalculate = set(cosmic_driver_df['uniprot'])

db = sql.MyDatabase(path_to_archive='/home/kimlab1/database_data/elaspic/human')


###############################################################################
### Make a flatfile with sasa score for every mutation, or a reason why it couldn't be calculated

mutations = cosmic_driver_df[['uniprot', 'location']].values
num_of_mutations = mutations.shape[0]
sasa_score = np.empty((num_of_mutations,1), dtype=float)
sasa_old_score = np.empty((num_of_mutations,1), dtype=float)
pfam_name = np.empty((num_of_mutations,1), dtype=object)
errors_log = np.empty((num_of_mutations,1),dtype=object)
uniprot_domain_mutations = dict()
#uniprot_definitions = db.get_uniprot_domain('P04637', copy_data=False)
for i, (idx, mutation_df) in enumerate(cosmic_driver_df.iterrows()):
    sasa = None
    pfam = None
    log = []
    
    # Uniprot domains (calculate sasa)
    uniprot_definitions = db.get_uniprot_domain(mutation_df.uniprot, copy_data=False)
    for idx, (uniprot_domain, uniprot_template, uniprot_model) in enumerate(uniprot_definitions):
        if not uniprot_template \
        or not uniprot_template.domain_def \
        or not uniprot_model \
        or (not uniprot_model.model_errors and not uniprot_model.model_filename):
            domain_start, domain_end = sql.decode_domain(uniprot_domain.alignment_def)
            if mutation_df.location < domain_start or mutation_df.location > domain_end:
                log.append(uniprot_domain.pfam_name + ': mutation outside of pfam domain')
            else:
                log.append(uniprot_domain.pfam_name + ': mutation falls inside pfam domain but no structural template was found')
        else:
            domain_start, domain_end = sql.decode_domain(uniprot_template.domain_def)
            if mutation_df.location < domain_start or mutation_df.location > domain_end:
                log.append(uniprot_domain.pfam_name + ': mutation outside of pdbfam domain')
            else:
                mutation_idx = domain_start - mutation_df.location - 1
                if uniprot_model.sasa_score is None:
                    log.append(uniprot_domain.pfam_name + ': mutation inside pdbfam domain but sasa score could not be calculated')
                else:
                    domain_length = domain_end - domain_start + 1
                    sasa_old = [score for score in uniprot_model.sasa_score.split(',')]
                    pfam = uniprot_domain.pfam_name
                    
                    analyze_structure = class_analyze_structure.AnalyzeStructure(
                        '/home/kimlab1/database_data/elaspic/' + uniprot_domain.path_to_data, 
                        '/tmp/elaspic/NUORFs/analyze_structure/',
                        uniprot_model.model_filename,
                        ['A'], None, None)
                    sasa = analyze_structure.get_sasa()[0]['A']
                    try:
                        assert len(sasa_old) == len(sasa)
                    except:
                        print idx, domain_length, len(sasa_old), len(sasa)
                        print sasa_old[mutation_idx], sasa[mutation_idx]
                        sasa_old = None
                        sasa = None
                        break
                    
                    sasa_old = sasa_old[mutation_idx]
                    sasa = sasa[mutation_idx]
                    
                    # log = [uniprot_domain.pfam_name + ': affected by mutation']
                    mutation = ''.join([str(x) for x in mutation_df[['fromAA','location','toAA']].values])
                    uniprot_domain_mutations.setdefault((uniprot_domain, uniprot_template, uniprot_model,), list())\
                        .append(mutation)
                    break
    sasa_score[i] = sasa
    sasa_old_score[i] = sasa_old
    pfam_name[i] = pfam
    if not pfam:
        errors_log[i] = '; '.join(log)
        for l in log:
            if 'mutation outside of pdbfam domain' not in l:
                pfam_name[i] = 'no structural template'
                break
            pfam_name[i] = 'mutation outside of domain'
    

cosmic_driver_df['sasa_score'] = sasa_score
cosmic_driver_df['pfam_name'] = pfam_name
cosmic_driver_df['errors_log'] = errors_log
cosmic_driver_df['uniprot_name'] = cosmic_driver_df['uniprot'].apply(lambda x: db.get_uniprot_sequence(x).name)
cosmic_driver_df['uniprot_description'] = cosmic_driver_df['uniprot'].apply(lambda x: db.get_uniprot_sequence(x).description.split(' OS=Homo ')[0])
cosmic_driver_df['mutation'] = cosmic_driver_df['fromAA'] + cosmic_driver_df['location'].apply(lambda x: str(x)) + cosmic_driver_df['toAA']

#print cosmic_driver_df

cosmic_driver_df_final = cosmic_driver_df[['uniprot', 'uniprot_name', 'uniprot_description', 'mutid', 'mutation', 'sasa_score', 'pfam_name']]
cosmic_driver_df_final.columns = ['uniprot_id', 'uniprot_name', 'uniprot_description', 'cosmic_mutation_id', 'mutation', 'sasa_score', 'pfam_domain_affected']
cosmic_driver_df_final.to_csv('driver_data_final.tsv', sep='\t', index=False)


###############################################################################



def show_pymol(domain_definition, list_of_mutations=None):
    """
    Colour residues by their sasa score instead of the b factor
    """
    d, t, m = domain_definition
    print d.uniprot_id, d.pfam_name
    
    path_to_pdb = '/home/alexey/elaspic/elaspic/' + d.path_to_data + m.model_filename
    structure_name = m.model_filename.strip('.pdb')
    
    sasa_score = [float(score) for score in m.sasa_score.split(',')]
    
    # Load the pdb
    print path_to_pdb
    pymol.cmd.reinitialize()
    pymol.cmd.load(path_to_pdb, structure_name)
    
    # clear out the old B Factors
    pymol.cmd.alter('%s'%structure_name, 'b=0.0')    
    
    # update the B Factors with new properties
    for idx, sasa in enumerate(sasa_score):
        resnum = idx + 1
        pymol.cmd.alter('%s and resi %i'%(structure_name,resnum), 'b=%f'%sasa)
        
    # color the protein based on the new B Factors
    pymol.cmd.spectrum('b', palette='rainbow', selection='(all)')
    
    pymol.cmd.hide('all')
    pymol.cmd.show('cartoon')
    pymol.cmd.set('surface_quality', 1)
    pymol.cmd.set('transparency', 0.6)
    pymol.cmd.set('sphere_transparency', 0.2)
    
    if list_of_mutations:
        print list_of_mutations
        for mutation in list_of_mutations:
            pymol.stored.mutation_a = None
            pymol.stored.mutation_label = None
            pymol.stored.structure_b = None
            pymol.stored.structure_aaa = None
            
            from_aa = mutation[0]
            to_aa = mutation[-1]
            uniprot_position = int(mutation[1:-1])
            # print uniprot_position
            # print sql.decode_domain(t.domain_def)[0]
            
            resnum = uniprot_position - sql.decode_domain(t.domain_def)[0] + 1
            # print resnum
            
            pymol.stored.mutation_a = from_aa.upper()
            pymol.cmd.iterate('%s and resi %i and n. CA'%(structure_name,resnum), \
                'print resn; stored.structure_aaa = resn; stored.structure_b = b')
            assert AAA_DICT[pymol.stored.structure_aaa] == pymol.stored.mutation_a
            
            pymol.stored.mutation_label = '%s_%.2f'%(mutation,pymol.stored.structure_b)
            pymol.cmd.label('%s and resi %i and n. CA'%(structure_name,resnum), "stored.mutation_label")
            pymol.cmd.show('sticks', '%s and resi %i'%(structure_name,resnum))
            pymol.cmd.show('spheres', '%s and resi %i'%(structure_name,resnum))

    
    pymol.cmd.show('labels')    
    pymol.cmd.set('label_size', -1)
    pymol.cmd.set('label_position',(3,1,5))
    
    pymol.cmd.refresh()
#    pymol.cmd.save(path_to_pdb.replace('.pdb','-sasa_bfactors.pdb'))
#    return (d.uniprot_id + '_' + d.pfam_name + '.png')

if False:
    idx = 6
    show_pymol(uniprot_domain_mutations.items()[idx][0], uniprot_domain_mutations.items()[idx][1])
    filename = (save_path +
                 uniprot_domain_mutations.items()[idx][0][0].uniprot_id + '*' +
                 uniprot_domain_mutations.items()[idx][0][0].pfam_name +
                 '.png')

if False:
    pymol.cmd.png(filename, ray=1, dpi=300)


################################################################################
#import json
#
#with open('/home/kimlab1/database_data/elaspic/human/Q9Y/6K/Q9Y6K1/PWWP*291-374/template.json', 'r') as fh: template_data = json.load(fh)
#
#template3 = sql.UniprotDomainTemplate(**template_data)
##template = sql.UniprotDomainTemplate()
#template3.uniprot_domain_id = 1000
##template.template_errors = 'testing2'
#db.add_uniprot_template(template3)
#
#db.close()



###############################################################################
###


def make_plot(mutation_stats, title_string=''):
    bar_labels = (" No template\n or errors ",
                  " Mutation outside\n domain ", 
    #              " Interface\n templates ", 
                  " Mutation in\n core ", 
                  " Mutation on\n surface ")
    n_groups = len(bar_labels)
    
    # Plot driver mutation statistics                       
    fig, ax = plt.subplots()
    
    index = np.arange(n_groups) + 0.35
    bar_width = 0.7
    opacity = 1 # 0.4
    # error_config = {'ecolor': '0.3'}
    
    rects1 = plt.bar(index, mutation_stats, bar_width, alpha=opacity, color='b')
    
    plt.ylabel('Number of mutations', size='large')
    plt.title(title_string, size='large')
    plt.xticks(index + bar_width/2, bar_labels, size='large', rotation=0, ha='center')
    plt.xlim(0, 4.35)
    
    plt.tight_layout()
    plt.show()
#    plt.savefig('/home/alexey/documents/presentations/subgroup-meetings/splicing/131204/%s-by-%s.png' % (mutation_type, binning_by), format='png', dpi=600)
#    plt.close()
    

cosmic_driver_df = pd.read_csv('driver_data.tsv', sep='\t')
#h = cosmic_driver_df['sasa_score'].hist(bins=22.8,range=(0.63,30))
#h.xaxis.label.set_size(16)
#h.yaxis.label.set_size(16)


def mutation_outside_domain(error_log):
    errors = [e.strip() for e in error_log.split(';')]
    outside_domain = True
    for error in errors:
        if 'mutation outside of pdbfam domain' not in error:
            outside_domain = False
    return outside_domain



###############################################################################
### Core mutations

number_of_missing = len(cosmic_driver_df[pd.isnull(cosmic_driver_df['sasa_score'])]['errors_log'].apply(mutation_outside_domain))
outside_domain = sum(cosmic_driver_df[pd.isnull(cosmic_driver_df['sasa_score'])]['errors_log'].apply(mutation_outside_domain))
template_error = number_of_missing - outside_domain
core = len(cosmic_driver_df[cosmic_driver_df['sasa_score'] <= 5])
surface = len(cosmic_driver_df[cosmic_driver_df['sasa_score'] > 5])

#make_plot([template_error, outside_domain, core, surface])

                    
###############################################################################
### Interface mutations

interface_mutations = dict()
uniprot_domain_pair_mutations = dict()
cosmic_driver_df_groups = cosmic_driver_df.groupby(['uniprot', 'pfam_name'])
for key, mutation_ids in cosmic_driver_df_groups.groups.iteritems():
    uniprot_id, pfam_name = key
    print uniprot_id, pfam_name
    
    domain_definitions = db.get_uniprot_domain_pair(uniprot_id, False)
    interacting_aa_set = set()
    for d, t, m in domain_definitions:
        if d.uniprot_domain_1.pfam_name != pfam_name and \
        d.uniprot_domain_2.pfam_name != pfam_name:
            continue
        elif d.uniprot_domain_1.uniprot_id == uniprot_id and \
        d.uniprot_domain_1.pfam_name == pfam_name and \
        m and not m.model_errors:
            domain_start = sql.decode_domain(t.domain_def_1)[0]
            interacting_aa = m.interacting_aa_1
        elif d.uniprot_domain_2.uniprot_id == uniprot_id and \
        d.uniprot_domain_2.pfam_name == pfam_name and \
        m and not m.model_errors:
            domain_start = sql.decode_domain(t.domain_def_2)[0]
            interacting_aa = m.interacting_aa_2
            
        if not interacting_aa:
            continue
        interacting_aa_list = [int(i) + domain_start - 1 for i in interacting_aa.split(',') if i and i != '']
    
        # print interacting_aa_list
            
        # print uniprot_id, pfam_name
        for mutation_id in mutation_ids:
            mutation_position = cosmic_driver_df.iloc[mutation_id]['location']
            if mutation_position in interacting_aa_list:
                
                interface_mutations.setdefault((uniprot_id, pfam_name,), list()).append((mutation_id, domain_definitions,))
                
                # (protein_definition) -> list of mutations affecting it
                mutation = ''.join([str(x) for x in cosmic_driver_df.iloc[mutation_idx]\
                    [['fromAA','location','toAA']].values])
                uniprot_domain_pair_mutations.setdefault((d,t,m,),list()).append(mutation)
                    
        print mutation_id

    print '------------------------------------------------------------\n'


num_of_mutations = 0
list_of_sasa = []
for v in interface_mutations.itervalues():    
    mutation_and_sasa = set([(w[0], cosmic_driver_df.iloc[w[0]]['sasa_score'],) for w in v])
    num_of_mutations += len(mutation_and_sasa)
    list_of_sasa.extend([m[1] for m in mutation_and_sasa])
    
h = plt.hist(list_of_sasa, bins=22.8, range=(0.63,30))
h.xaxis.label.set_size(16)
h.yaxis.label.set_size(16)
###############################################################################






  
    
#        else:
#            print 'Error!!!'
#            print (d.uniprot_domain_1.uniprot_id,
#                   d.uniprot_domain_1.pfam_name,
#                   d.uniprot_domain_2.uniprot_id,
#                   d.uniprot_domain_2.pfam_name)
