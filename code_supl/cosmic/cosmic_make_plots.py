import pandas as pd
import cPickle as pickle

import numpy as np
import matplotlib.pyplot as plt
import SOAPpy
from SOAPpy import WSDL

biodbnet = WSDL.Proxy('http://biodbnet.abcc.ncifcrf.gov/webServices/bioDBnet.wsdl')

###############################################################################
# Use bioDBnet to get the uniprotID for all of the cosmic IDs

#in_fh = open('/home/alexey/working/databases/cosmic/accessions.tsv', 'r')
#accession_list = [line.strip() for line in in_fh.readlines()]
#in_fh.close()
#
#dbFindResultFull = ''
#i_start = 0
#for i_end in range(5000, len(accession_list), 5000):
#    website_query = ', '.join(accession_list[i_start:i_end])
#    output = 'UniProt Accession'
#    taxonId = '9606'
#    dbFindParams = SOAPpy.structType({'output': output, 'inputValues': website_query, 'taxonId': taxonId})
#    dbFindResult = biodbnet.dbFind(dbFindParams)
#    dbFindResultFull = dbFindResultFull + dbFindResult
#    i_start = i_end
#
#i_end = len(accession_list)
#website_query = ', '.join(accession_list[i_start:i_end])
#output = 'UniProt Accession'
#taxonId = '9606'
#dbFindParams = SOAPpy.structType({'output': output, 'inputValues': website_query, 'taxonId': taxonId})
#dbFindResult = biodbnet.dbFind(dbFindParams)
#dbFindResultFull = dbFindResultFull + dbFindResult
#    
#accession_db = dict()
#for line in dbFindResultFull.split('\n'):
#    rows = line.split('\t')
#    if rows[0] == 'Input Value':
#        continue
#    accession_db[rows[0]] = rows[-1].split(';')    
#    
#out_fh = open('/home/alexey/working/databases/cosmic/accessions.pickle', 'wb')
#pickle.dump(accession_db, out_fh)
#out_fh.close()

# Accessions as ENST or Refseq (cosmic)
with open('/home/alexey/working/databases/cosmic/accessions.pickle', 'rb') as fh:
    accession_db = pickle.load(fh)

###############################################################################

# Load Niklas' databases
handle = open('/home/alexey/working/pipeline-splicing/code/pipeline_3DID_database.pickle', 'rb')
pipeline_3DID_database = pickle.load(handle)
handle.close()

handle = open('/home/alexey/working/pipeline-splicing/code/pipeline_core_templates.pickle', 'rb')
pipeline_core_templates = pickle.load(handle)
handle.close()

handle = open('/home/alexey/working/pipeline-splicing/code/pipeline_interaction_database.pickle', 'rb')
pipeline_interaction_database = pickle.load(handle)
handle.close()

handle = open('/home/alexey/working/pipeline-splicing/code/pipeline_uniprot_sequences.pickle', 'rb')
pipeline_uniprot_sequences = pickle.load(handle)
handle.close()

###############################################################################

class build_cosmic_db():
    
    def __init__(self):
        self.db = dict()
    
    
    def add_mutation(self, cosmic_accession, cosmic_mutation):
        assert cosmic_mutation[:2] == 'p.'
        mutation = cosmic_mutation[2:]
        
        if self.db.has_key(cosmic_accession):
            self.db[cosmic_accession]['mutations'].append(mutation)
        else:
            self.__create_db_entry(cosmic_accession)
            self.db[cosmic_accession]['mutations'].append(mutation)
            
    
    def __create_db_entry(self, cosmic_accession):
        uniprotIDs = accession_db[cosmic_accession]
        if uniprotIDs[0] == '-':
            stripped_accession = cosmic_accession.split('.')[0]
            for version in range(5):
                uniprotIDs = self.__find_uniprotIDs(stripped_accession + '.' + str(version))
                if not uniprotIDs[0] == '-':
                    break
        
        self.db[cosmic_accession] = {'uniprotIDs': uniprotIDs, 'mutations': []}
    
    
    def __find_uniprotIDs(self, cosmic_accession):
        inputValues = cosmic_accession
        output = 'UniProt Accession'
        taxonId = '9606'
        
        dbFindParams = SOAPpy.structType({'output': output, 'inputValues': inputValues, 
        				    'taxonId': taxonId})
        dbFindResult = biodbnet.dbFind(dbFindParams)
        
        return dbFindResult.split('\n')[1].split('\t')[-1].split(';')
        
    def save_db(self, filename):
        out_fh = open(filename, 'wb')
        pickle.dump(self.db, out_fh)
        out_fh.close()

###############################################################################
# Make (or open) lists of driver and passenger ids

#drivers = build_cosmic_db()
#with open('/home/alexey/working/databases/cosmic/drivers.tsv', 'r') as in_fh:
#    for counter, line in enumerate(in_fh):
#        if counter % 1000 == 0:
#            print counter
#        rows = line.strip().split('\t')
#        drivers.add_mutation(rows[1],rows[2])
#drivers.save_db('/home/alexey/working/databases/cosmic/drivers.pickle')
#
#
#passengers = build_cosmic_db()
#with open('/home/alexey/working/databases/cosmic/passengers.tsv', 'r') as in_fh:
#    for counter, line in enumerate(in_fh):
#        if counter % 1000 == 0:
#            print counter
#        rows = line.strip().split('\t')
#        passengers.add_mutation(rows[1],rows[2])
#passengers.save_db('/home/alexey/working/databases/cosmic/passengers.pickle')


drivers = list()
with open('/home/alexey/working/databases/cosmic/drivers.tsv', 'r') as fh:
    for counter, line in enumerate(fh):
        if counter % 1000 == 0:
            print counter
        rows = line.strip().split('\t')
        drivers.append(rows[1])

passengers = list()
with open('/home/alexey/working/databases/cosmic/passengers.tsv', 'r') as fh:
    for counter, line in enumerate(fh):
        if counter % 100000 == 0:
            print counter
        rows = line.strip().split('\t')
        passengers.append(rows[1])

###############################################################################

# Set parameters to be used for 
# binning_by options: 'proteins', 'mutations'
do_positions = True
do_proteins = False
binning_by = 'mutations'
if binning_by == 'proteins':
    drivers = list(set(drivers))
    passengers = list(set(passengers))

###############################################################################

def categorise_mutations(list_of_ids):  
    # Go over all mutations
    #mut_in_interface = 0 # If interacting, must have both
    mut_no_uniprot = 0
    mut_in_core = 0
    mut_in_both = 0
    mut_in_neither = 0
    
    for key in list_of_ids:
        
        no_uniprot = False
        is_interacting = False
        is_core = False
        
        for uniprotID in accession_db[key]:
            if uniprotID == '-':
                no_uniprot = True
            if pipeline_interaction_database.has_key(uniprotID):
                is_interacting = True
            if pipeline_core_templates.has_key(uniprotID):
                is_core = True
        
        if no_uniprot:
            mut_no_uniprot += 1
        elif is_interacting and not is_core:
    #        driver_interacting += 1
            mut_in_both += 1
        elif not is_interacting and is_core:
            mut_in_core += 1
        elif is_interacting and is_core:
            mut_in_both += 1
        elif not is_interacting and not is_core:
            mut_in_neither += 1
        else:
            raise Exception("Something went wrong")
            
    return mut_no_uniprot, mut_in_neither, mut_in_core, mut_in_both


def make_plot(mutation_type, mutation_stats, title_string):
    bar_labels = (" Can't map\n to uniprot ",
                  " No structural\n templates ", 
    #              " Interface\n templates ", 
                  " Only core\n templates ", 
                  " Interface\n and core\n templates ")
    n_groups = len(bar_labels)
    
    # Plot driver mutation statistics                       
    fig, ax = plt.subplots()
    
    index = np.arange(n_groups) + 0.35
    bar_width = 0.7
    opacity = 0.4
    error_config = {'ecolor': '0.3'}
    
    rects1 = plt.bar(index, mutation_stats, bar_width, alpha=opacity, color='b')
    
    plt.ylabel('Number of %s' % binning_by, size='large')
    plt.title(title_string, size='large')
    plt.xticks(index + bar_width/2, bar_labels, size='large', rotation=0, ha='center')
    plt.xlim(0, 4.35)
    
    plt.tight_layout()
    plt.show()
    plt.savefig('/home/alexey/documents/presentations/subgroup-meetings/splicing/131204/%s-by-%s.png' % (mutation_type, binning_by), format='png', dpi=600)
    plt.close()

if do_proteins:
    driver_stats = categorise_mutations(drivers)
    make_plot('Driver', driver_stats, 'Driver mutations in COSMIC (> 5 records, > 3 pubmeds)')
    
    passenger_stats = categorise_mutations(passengers)
    make_plot('Passenger', passenger_stats, 'Passenger mutations in COSMIC (< 5 records, < 3 pubmeds)')


###############################################################################
# Accessions as Uniprot (Recep)
###############################################################################

with open('/home/alexey/working/pipeline/code_supl/cosmic/CosmicV67.txt', 'r') as fh:
    mutations_df = pd.read_csv(fh, sep='\t')

def categorise_mutations(mutation_type):  
    # Go over all mutations
    #mut_in_interface = 0 # If interacting, must have both
    templates_core = 0
    templates_both = 0
    templates_neither = 0
    mut_in_interface = 0
    
    for mType, mutid, fromAA, toAA, location, uniprot in mutations_df[mutations_df['mType'] == mutation_type].itertuples(index=False):
        
        is_interface = False
        is_core = False
        is_interacting = False
        
        if pipeline_interaction_database.has_key(uniprot):
            is_interface = True
            for interaction in pipeline_interaction_database[uniprot]:
                if int(location) in [int(aa.strip()) for aa in interaction[4].split(',') if interaction[4] != '']:
                    is_interacting = True
        if pipeline_core_templates.has_key(uniprot):
            is_core = True
            
        # Core, interface, both, neither...
        if is_interacting:
            mut_in_interface += 1
        if is_interface and not is_core:
    #        driver_interacting += 1
            templates_both += 1
        elif not is_interface and is_core:
            templates_core += 1
        elif is_interface and is_core:
            templates_both += 1
        elif not is_interface and not is_core:
            templates_neither += 1
        else:
            raise Exception("Something went wrong")
            
    return templates_neither, templates_core, templates_both, mut_in_interface


def make_plot(mutation_type, mutation_stats, title_string):
    bar_labels = (" No structural\n templates ", 
    #              " Interface\n templates ", 
                  " Only core\n templates ", 
                  " Interface\n and core\n templates ",
                  " Mutation in\n interface ")
    n_groups = len(bar_labels)
    
    # Plot driver mutation statistics                       
    fig, ax = plt.subplots()
    
    index = np.arange(n_groups) + 0.35
    bar_width = 0.7
    opacity = 0.4
    error_config = {'ecolor': '0.3'}
    
    rects1 = plt.bar(index, mutation_stats, bar_width, alpha=opacity, color='b')
    
    plt.ylabel('Number of mutations', size='large')
    plt.title(title_string, size='large')
    plt.xticks(index + bar_width/2, bar_labels, size='large', rotation=0, ha='center')
    plt.xlim(0, 4.35)
    
    plt.tight_layout()
    plt.show()
    plt.savefig('/home/alexey/documents/presentations/subgroup-meetings/splicing/131204/%s-by-mutations-positions.png' % mutation_type,
                format='png', dpi=600)
    plt.close()

if do_positions:
    driver_stats = categorise_mutations('Driver')
    make_plot('Driver', driver_stats, 'Driver mutations in COSMIC by position')
        
    passenger_stats = categorise_mutations('Passenger')
    make_plot('Passenger', passenger_stats, 'Passenger mutations in COSMIC by position')                