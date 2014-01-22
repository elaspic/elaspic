# -*- coding: utf-8 -*-
"""
Created on Wed Jun 19 12:01:51 2013

@author: niklas
"""
import os

def getRSA(SASA, aminoAcid):
    """ see: 
        http://prowl.rockefeller.edu/aainfo/volume.htm
        http://www.biostars.org/p/10650/
        
        RSA is a dict with 'Residue Volume' as first entry and 'Surface Area'
        as second enty with key amino acid in one letter code
    """
    RSA = {'A': (88.6, 115), \
           'R': (173.4, 225), \
           'D': (111.1, 150), \
           'N': (114.1, 160), \
           'C': (108.5, 135), \
           'E': (138.4, 190), \
           'Q': (143.8, 180), \
           'G': (60.1, 75), \
           'H': (153.2, 195), \
           'I': (166.7, 175), \
           'L': (166.7, 170), \
           'K': (168.6, 200), \
           'M': (162.9, 185), \
           'F': (189.9, 210), \
           'P': (112.7, 145), \
           'S': (89.0, 115), \
           'T': (116.1, 140), \
           'W': (227.8, 255), \
           'Y': (193.6, 230), \
           'V': (140.0, 155) \
           }
    
#    if float(SASA) == 0:
#        return 'NaN'

    if len(aminoAcid) == 3:
        aa = convert_aa(aminoAcid)
        result = float(SASA) / RSA[aa][1]
        return '{:.2f}'.format(result)
    elif len(aminoAcid) == 1:
        result = float(SASA) / RSA[aminoAcid][1]
        return '{:.2f}'.format(result)
    else:
        print 'Problem with the calculation of the RSA'
        return 'NaN' 


def convert_aa(aa):
    A_DICT = {'A':'ALA', 'R':'ARG', 'N':'ASN', 'D':'ASP', 'C':'CYS', 'E':'GLU', \
              'Q':'GLN', 'G':'GLY', 'H':'HIS', 'I':'ILE', 'L':'LEU', 'K':'LYS', \
              'M':'MET', 'F':'PHE', 'P':'PRO', 'S':'SER', 'T':'THR', 'W':'TRP', \
              'Y':'TYR', 'V':'VAL', 'U':'SEC', 'O':'PYL', \
              'B':'ASX', 'Z':'GLX', 'J':'XLE', 'X':'XAA', '*':'TER'}
    
    AAA_DICT = dict([(value,key) for key,value in A_DICT.items()])
    
    if len(aa) == 3:
        try:
            return AAA_DICT[aa.upper()]
        except KeyError:
            print  'Not a valid amino acid'
            return
    if len(aa) == 1:
        try:
            return A_DICT[aa.upper()]
        except KeyError:
            print  'Not a valid amino acid'
            return
    print 'Not a valid amino acid'


def convert_mutations(mut):
    fromAA = mut[:3]
    position = mut[3:-3]
    toAA = mut[-3:]
    
    fromAA = convert_aa(fromAA)
    toAA = convert_aa(toAA)
    
    return fromAA + position + toAA


def read_mutations(path, \
                   computed_success_wt = dict(), \
                   computed_success_mut = dict(), \
                   computed_additional_information = dict(), \
                   computed_error = set(), \
                   computed_timeout = set(), \
                   computed_success = set() \
                   ):
    """ read all the calculated mutations from one folder
    """
    # store the IDs here
    computed_no_template = list()
#    computed_error = set()
#    computed_timeout = set()
#    computed_success = set()

    # from the no_processed list, check which failed or timed out
    files = [ path + '/' + item + '/not_processed.log' for item in os.listdir(path) ]
    for fname in files:
        with open(fname, 'r') as f:
            for line in f:
                ID, reason = line.split('\t')
                if reason.strip() == 'no template found':
                    computed_no_template.append(ID)
                elif reason.strip() == 'timeout':
                    computed_timeout.add(ID)
                else:
                    computed_error.add(ID)
    
    ## from the succesfull ones, read the data
    # result_wt
    files = [ path + '/' + item + '/result_wt.log' for item in os.listdir(path) ]
    for fname in files:
        with open(fname, 'r') as f:
            f.readline()
            for l in f:
                line = l.split('\t')
                ID = line[0][:-3]
                if computed_success_wt.has_key(ID):
                    if computed_success_wt[ID] != line:
                        print 'Duplicate wt result with differing prediction in', ID
                computed_success_wt[ID] = line
                computed_success.add(ID)
                
    # result_mut
    files = [ path + '/' + item + '/result_mut.log' for item in os.listdir(path) ]
    for fname in files:
        with open(fname, 'r') as f:
            f.readline()
            for l in f:
                line = l.split('\t')
                ID = line[0][:-4]
                if computed_success_mut.has_key(ID):
                    if computed_success_mut[ID] != line:
                        print 'Duplicate mut result with differing prediction in', ID
                computed_success_mut[ID] = line
    
    # additional info
    files = [ path + '/' + item + '/result_additional_information.log' for item in os.listdir(path) ]
    for fname in files:
        with open(fname, 'r') as f:
            f.readline()
            for l in f:
                line = l.split('\t')
                ID = line[0]
                if computed_additional_information.has_key(ID):
                    if computed_additional_information[ID] != line:
                        print 'Duplicate additional info result with differing prediction in', ID
                computed_additional_information[ID] = line

#    with open('failed.dat', 'a') as f:
#        for item in computed_error:
#            uniprotKB, mutation = item.split('_')
#            f.write(uniprotKB + '\t' + mutation + '\n')
#    with open('timeout.dat', 'a') as f:
#        for item in computed_timeout:
##            f.write(item[:4] + '\t' + item[5:] + '\n')
#            f.write(item + '\n')
#    with open('success.dat', 'a') as f:
#        for item in computed_success:
##            f.write(item[:4] + '\t' + item[5:] + '\n')
#            f.write(item + '\n')

    return computed_success_wt, computed_success_mut, computed_additional_information, \
            computed_error, computed_timeout, computed_success








def read_sift_predictions(ID):
    """ read all the SIFT predictions and store them in a dict
    """
    path = '/home/niklas/work-ccbr/working-folder/paper/databases/sift_predictions/'
    
    uniprotKB         = ID.split('_')[0]
    try:
        fromAA, pos, toAA = ID.split('_')[1][0], ID.split('_')[1][1:-1], ID.split('_')[1][-1]
    except:
        print ID
        raise
    
    sift_failed = ['P10762', 'P30289', 'P03040', 'P03050', \
                     'P37001', 'P03706', 'P07445', 'O14836', \
                     'P21817', 'P21817', 'A6NMZ7', 'P29400', \
                     'P53420', 'Q01955', 'Q6NW29', 'A2RTY3', \
                     'O15050', 'O15287', 'O43508', 'O75534', \
                     'O95715', 'P08572', 'P12109', 'P15018', \
                     'P25940', 'P32971', 'P39059', 'P59894', \
                     'Q07092']
    if uniprotKB in sift_failed:
        return '-1'

    columns = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K',\
               'L', 'M', 'N', 'P', 'Q', 'R', 'S', 'T', 'V', 'W',\
               'X', 'Y', 'Z', '*', '-']
    get_column = columns.index(toAA)
    try:
        with open(path + uniprotKB + '.SIFTprediction', 'r') as f:
            # skip the header
            line = f.readline()
            i = 0
            while line[0:3] != '  A':
                line = f.readline()
                if i > 100:
                    return '-1'
                i += 1
            #go to the correct amino acid
            # i.e. one line before the line that we want to read in
            for i in range(int(pos)):
                line = f.readline()
            
            # read the line and get rid of the blanks
            line = [ item.strip() for item in line.split(' ') if item.strip() != '' ]
            
            sift_score = line[get_column]
    except IOError:
        return 'SIFT file not found'

    return sift_score

def read_line(line, core_or_interface):
    
    #ID = line[0]
    normDOPE = line[1]
    intraclashesEnergy1_wt = line[2]
    intraclashesEnergy2_wt = line[3]
    interactionEnergy_wt = line[4]
    Backbone_Hbond_wt = line[5]
    Sidechain_Hbond_wt = line[6]
    Van_der_Waals_wt = line[7]
    Electrostatics_wt = line[8]
    Solvation_Polar_wt = line[9]
    Solvation_Hydrophobic_wt = line[10]
    Van_der_Waals_clashes_wt = line[11]
    entropy_sidechain_wt = line[12]
    entropy_mainchain_wt = line[13]
    sloop_entropy_wt = line[14]
    mloop_entropy_wt = line[15]
    cis_bond_wt = line[16]
    torsional_clash_wt = line[17]
    backbone_clash_wt = line[18]
    helix_dipole_wt = line[19]
    water_bridge_wt = line[20]
    disulfide_wt = line[21]
    electrostatic_kon_wt = line[22]
    partial_covalent_bonds_wt = line[23]
    energy_Ionisation_wt = line[24]
    Entropy_Complex_wt = line[25]
    Number_of_Residues_wt = line[26]
    stabilityEnergy_wt = line[27]
    stability_Backbone_Hbond = line[28]
    stability_Sidechain_Hbond = line[29]
    stability_Van_der_Waals = line[30]
    stability_Electrostatics = line[31]
    stability_Solvation_Polar = line[32]
    stability_Solvation_Hydrophobic = line[33]
    stability_Van_der_Waals_clashes = line[34]
    stability_entropy_sidechain = line[35]
    stability_entropy_mainchain = line[36]
    stability_sloop_entropy = line[37]
    stability_mloop_entropy = line[38]
    stability_cis_bond = line[39]
    stability_torsional_clash = line[40]
    stability_backbone_clash = line[41]
    stability_helix_dipole = line[42]
    stability_water_bridge = line[43]
    stability_disulfide = line[44]
    stability_electrostatic_kon = line[45]
    stability_partial_covalent_bonds = line[46]
    stability_energy_Ionisation = line[47]
    stability_Entropy_Complex = line[48]
    stability_Number_of_Residues = line[49].strip()
    
    core_or_interface = 'core'
     
    if core_or_interface == 'interface':
        return normDOPE, \
               intraclashesEnergy1_wt, intraclashesEnergy2_wt, interactionEnergy_wt, '0', \
               Backbone_Hbond_wt, Sidechain_Hbond_wt, \
               Van_der_Waals_wt, Electrostatics_wt, \
               Solvation_Polar_wt, Solvation_Hydrophobic_wt, \
               Van_der_Waals_clashes_wt, entropy_sidechain_wt, \
               entropy_mainchain_wt, sloop_entropy_wt, \
               mloop_entropy_wt, cis_bond_wt, \
               torsional_clash_wt, backbone_clash_wt, \
               helix_dipole_wt, water_bridge_wt, \
               disulfide_wt, electrostatic_kon_wt, \
               partial_covalent_bonds_wt, energy_Ionisation_wt, \
               Entropy_Complex_wt, Number_of_Residues_wt
    elif core_or_interface == 'core':
        return normDOPE, \
               '0', '0', '0', stabilityEnergy_wt, \
               stability_Backbone_Hbond, stability_Sidechain_Hbond, \
               stability_Van_der_Waals, stability_Electrostatics, \
               stability_Solvation_Polar, stability_Solvation_Hydrophobic, \
               stability_Van_der_Waals_clashes, stability_entropy_sidechain, \
               stability_entropy_mainchain, stability_sloop_entropy, \
               stability_mloop_entropy, stability_cis_bond, \
               stability_torsional_clash, stability_backbone_clash, \
               stability_helix_dipole, stability_water_bridge, \
               stability_disulfide, stability_electrostatic_kon, \
               stability_partial_covalent_bonds, stability_energy_Ionisation, \
               stability_Entropy_Complex, stability_Number_of_Residues
    
    

def make_output(line_result_wt, line_result_mut, line_additional_info, polyphen, mutations):
    """ combines the input for the mutations and returns the line
    which should be written to a file
    """
    ID                         = line_additional_info[0]
    core_or_interface          = line_additional_info[1]
#    core_or_interface          = 'interface'
    seq_id_avg                 = line_additional_info[2]
    seq_id_chain1              = line_additional_info[3]
    seq_id_chain2              = line_additional_info[4]
    matrix_score               = line_additional_info[5]
    if_hydrophobic             = line_additional_info[6]
    if_hydrophilic             = line_additional_info[7]
    if_total                   = line_additional_info[8]
    contactVector_wt_ownChain  = line_additional_info[9].split(':')
    contactVector_wt           = line_additional_info[10].split(':')
    contactVector_mut_ownChain = line_additional_info[11].split(':')
    contactVector_mut          = line_additional_info[12].split(':')
    solvent_accessibility_wt   = line_additional_info[13]
    secondary_structure_wt     = line_additional_info[14]
    solvent_accessibility_mut  = line_additional_info[15]
    secondary_structure_mut    = line_additional_info[16].strip()
    
    disease_types = mutations[ID]
    
#    print 'ID', ID
#    print 'disease_type', disease_type
    
    
#    H = alpha helix
#    B = residue in isolated beta-bridge
#    E = extended strand, participates in beta ladder
#    G = 3-helix (3/10 helix)
#    I = 5 helix (pi helix)
#    T = hydrogen bonded turn
#    S = bend 
    secondary_structure_mut_count = 0
    
    if secondary_structure_wt == '-':
        secondary_structure_wt = '0'
    elif secondary_structure_wt == 'H':
        secondary_structure_wt = '1'
    elif secondary_structure_wt == 'B':
        secondary_structure_wt = '2'
    elif secondary_structure_wt == 'E':
        secondary_structure_wt = '3'
    elif secondary_structure_wt == 'G':
        secondary_structure_wt = '4'
    elif secondary_structure_wt == 'I':
        secondary_structure_wt = '5'
    elif secondary_structure_wt == 'T':
        secondary_structure_wt = '6'
    elif secondary_structure_wt == 'S':
        secondary_structure_wt = '7'
    else:
#        print 'secondary_structure_wt unrecognized', secondary_structure_wt
#        with open('data_DSSP.dat', 'a') as f:
#            f.write('\t'.join(ID.split('_')) + '\n')
        secondary_structure_wt = '-1'
    
    if secondary_structure_mut == '-':
        secondary_structure_mut = '0'
    elif secondary_structure_mut == 'H':
        secondary_structure_mut = '1'
    elif secondary_structure_mut == 'B':
        secondary_structure_mut = '2'
    elif secondary_structure_mut == 'E':
        secondary_structure_mut = '3'
    elif secondary_structure_mut == 'G':
        secondary_structure_mut = '4'
    elif secondary_structure_mut == 'I':
        secondary_structure_mut = '5'
    elif secondary_structure_mut == 'T':
        secondary_structure_mut = '6'
    elif secondary_structure_mut == 'S':
        secondary_structure_mut = '7'
    else:
#        print 'secondary_structure_mut unrecognized', secondary_structure_mut
        secondary_structure_mut = '-1'
        secondary_structure_mut_count = 1
        
        
    
    wt_line  = read_line(line_result_wt, core_or_interface)
    mut_line = read_line(line_result_mut, core_or_interface)
            
#    additional_info[ID] = if_hydrophobic, if_hydrophilic, if_total, contactVector_own_chain_wt, contactVector_wt, contactVector_own_chain_mut, contactVector_mut, matrix_score
    
#    with open('/home/niklas/work-ccbr/playground/SKEMPI_full/SIFT/results/map_SKEMPI_full_mutation_pdb_to_sequence.dat', 'r') as f:
#        for line in f:
#            ID_map, pos_map = line.split('\t')
#            if ID == ID_map:
#                position = pos_map.strip()
#    ID_SIFT = ID[:4] + ID.split('_')[-1][1] + '_' + ID.split('_')[3][0] + position + ID[-1]
    
#    polyphen_missing = set()
#    # get the polyphen prediction
#    polphen_prediciton, polphen_dscore = polyphen(ID_SIFT)
#    if polphen_dscore == '?':
#        polphen_dscore = '-1'
#    
#    if polphen_prediciton == 'unknown':
#        polphen_prediciton = '1'
#    elif polphen_prediciton == 'benign':
#        polphen_prediciton = '2'
#    elif polphen_prediciton == 'possibly damaging':
#        polphen_prediciton = '3'
#    elif polphen_prediciton == 'probably damaging':
#        polphen_prediciton = '4'
#    elif polphen_prediciton == 'PolyPhen prediction not found':
##        print 'polyphen predictions not found for', ID
#        polphen_prediciton, polphen_dscore = '-1', '-1'
#        polyphen_missing.add(ID)
#    else:
#        polphen_prediciton = '-1'

    sift_missing = set()
    
    sift_score = read_sift_predictions(ID)

    if sift_score == '-':
        sift_score = '-1'
    elif sift_score == 'SIFT file not found':
#        print 'SIFT score not found for', ID
        sift_score = '-1'
        sift_missing.add(ID)
    
#    ddG_exp = ''
#    with open('/home/niklas/work-ccbr/playground/SKEMPI_full/SKEMPI_full_experimental_values.dat', 'r') as f:
#        f.readline()
#        for line in f:
#            ID_tmp, ddG_tmp = line.split('\t')
#            if ID == ID_tmp:
#                ddG_exp = float(ddG_tmp)
#    if ddG_exp == '':
#        print 'No exp. value found for:', ID
#        ddG_exp = 0
    ddG_exp = '0.0'
    
    
#    if core_or_interface == 'core':
    if True:
        dG_wt = float(wt_line[4])
        dG_mut = float(mut_line[4])
#        core_or_interface = '0'
        core_or_interface = '1'
    elif core_or_interface == 'interface':
        dG_wt = float(wt_line[3])
        dG_mut = float(mut_line[3])
        core_or_interface = '1'
    else:
        print 'core_or_interface not recognised', core_or_interface
    
#    ddG_calc = float(mut_line[3]) - float(wt_line[3])
    ddG_calc = dG_mut - dG_wt
    
    line = ''
    for disease_type in disease_types:
        line = line + ID + '\t' + ddG_exp + '\t' + "%.2f" % ddG_calc + '\t' + \
               "%.2f" % dG_wt + '\t' + "%.2f" % dG_mut + '\t' +\
               wt_line[1] + '\t' + wt_line[2] + '\t' + \
               '\t'.join(wt_line[5:-1]) + '\t' + \
               mut_line[1] + '\t' + mut_line[2] + '\t' + \
               '\t'.join(mut_line[5:]) + '\t' + \
               matrix_score + '\t' + \
               str(if_hydrophobic) + '\t' + str(if_hydrophilic) + '\t' + str(if_total) + '\t' + \
               '\t'.join(contactVector_wt) + '\t' + '\t'.join(contactVector_wt_ownChain) + '\t' + \
               '\t'.join(contactVector_mut) + '\t' + '\t'.join(contactVector_mut_ownChain) + '\t' + \
               secondary_structure_wt + '\t' + solvent_accessibility_wt + '\t' + \
               secondary_structure_mut + '\t' + solvent_accessibility_mut + '\t' + \
               sift_score + '\t' + "%.2f" % float(seq_id_avg) +'\t' + "%.2f" % float(seq_id_chain1) + '\t' + \
               "%.2f" % float(seq_id_chain2) + '\t' + wt_line[0] + '\t' + \
               core_or_interface + '\t' + disease_type + \
               '\n'
    
    # skip entries with normDOPE > 0
#    if float(wt_line[0]) > 0:
#        return '', sift_missing, polyphen_missing, secondary_structure_mut_count

    return line, sift_missing, polyphen_missing, secondary_structure_mut_count





mutations = dict()
#with open('/home/niklas/work-ccbr/playground/mutations_clasified_recep/mutations.txt', 'r') as f:
#    f.readline()
#    for l in f:
#        line = l.split('\t')
#        uniprotKB = line[0]
#        disease_type = line[6]
#        mutation_code = line[2] + line[1] + line[7]
#        mutations.setdefault(uniprotKB + '_' + mutation_code, set()).add(disease_type)
#
#with open('/home/niklas/work-ccbr/playground/mutations_clasified_recep/NonRedundantGWASMissense.txt', 'r') as f:
#    for l in f:
#        line = l.split('\t')
#        uniprotKB = line[1]
#        mutation_code = convert_mutations(line[2].strip().split('.')[1])
#        disease_type = 'GWAS'
#        mutations.setdefault(uniprotKB + '_' + mutation_code, set()).add(disease_type)
#
#i = 0
#with open('/home/niklas/work-ccbr/playground/data_run_with_dssp/disease/input_data/data_HAPMAP.dat', 'r') as f:
#    for l in f:
#        line = l.split('\t')
#        uniprotKB, mutation_code = line
#        disease_type = 'HAPMAP'
##        if i < 2:
##            print l
##            i += 1
#        mutations.setdefault(uniprotKB + '_' + mutation_code.strip(), set()).add(disease_type)

with open('disease_full_predictions.dat_send_to_recep', 'r') as f:
    h = f.readline()
    for l in f:
        line = l.split('\t')
        ID = line[0]
        disease_type = line[-1].strip()
        
        mutations.setdefault(ID, set()).add(disease_type)

# get a dict with the polyphen predictions
#polyphen = read_polyphen_predictions(['/home/niklas/work-ccbr/playground/mutations_clasified_recep/prepare_input/polyphen/v0/pph2-short.txt', \
#                                     '/home/niklas/work-ccbr/playground/mutations_clasified_recep/prepare_input/polyphen/v1/pph2-short.txt', \
#                                     '/home/niklas/work-ccbr/playground/mutations_clasified_recep/prepare_input/polyphen/v2/pph2-short.txt', \
#                                     '/home/niklas/work-ccbr/playground/mutations_clasified_recep/prepare_input/polyphen/v3/pph2-short.txt', \
#                                     '/home/niklas/work-ccbr/playground/mutations_clasified_recep/prepare_input/polyphen/v4/pph2-short.txt', \
#                                     '/home/niklas/work-ccbr/playground/SKEMPI_full/polyphen/v2/pph2-short.txt']
#                                    )
polyphen = dict()

# read all the mutations
#computed_success_wt, computed_success_mut, computed_additional_information, \
#computed_error, computed_timeout, computed_success = read_mutations('data/GWAS/')
#
#
#computed_success_wt, computed_success_mut, computed_additional_information, \
#computed_error, computed_timeout, computed_success = read_mutations('data/GWAS-2nd/', \
#                                                                    computed_success_wt, \
#                                                                    computed_success_mut, \
#                                                                    computed_additional_information, \
#                                                                    computed_error, \
#                                                                    computed_timeout, \
#                                                                    computed_success)
#                                                                                            
#computed_success_wt, computed_success_mut, computed_additional_information, \
#computed_error, computed_timeout, computed_success = read_mutations('data/COSMIC_DRIVER/', \
#                                                                    computed_success_wt, \
#                                                                    computed_success_mut, \
#                                                                    computed_additional_information, \
#                                                                    computed_error, \
#                                                                    computed_timeout, \
#                                                                    computed_success)
#
#computed_success_wt, computed_success_mut, computed_additional_information, \
#computed_error, computed_timeout, computed_success = read_mutations('data/COSMIC_DRIVER-2nd/', \
#                                                                    computed_success_wt, \
#                                                                    computed_success_mut, \
#                                                                    computed_additional_information, \
#                                                                    computed_error, \
#                                                                    computed_timeout, \
#                                                                    computed_success)

computed_success_wt, computed_success_mut, computed_additional_information, \
computed_error, computed_timeout, computed_success = read_mutations('data/interfaceStability/v1/')

computed_success_wt, computed_success_mut, computed_additional_information, \
computed_error, computed_timeout, computed_success = read_mutations('data/interfaceStability/v2/', \
                                                                    computed_success_wt, \
                                                                    computed_success_mut, \
                                                                    computed_additional_information, \
                                                                    computed_error, \
                                                                    computed_timeout, \
                                                                    computed_success)

computed_success_wt, computed_success_mut, computed_additional_information, \
computed_error, computed_timeout, computed_success = read_mutations('data/interfaceStability/beagle07/', \
                                                                    computed_success_wt, \
                                                                    computed_success_mut, \
                                                                    computed_additional_information, \
                                                                    computed_error, \
                                                                    computed_timeout, \
                                                                    computed_success)

computed_success_wt, computed_success_mut, computed_additional_information, \
computed_error, computed_timeout, computed_success = read_mutations('data/interfaceStability/v4/', \
                                                                    computed_success_wt, \
                                                                    computed_success_mut, \
                                                                    computed_additional_information, \
                                                                    computed_error, \
                                                                    computed_timeout, \
                                                                    computed_success)

computed_success_wt, computed_success_mut, computed_additional_information, \
computed_error, computed_timeout, computed_success = read_mutations('data/interfaceStability/v4.1/', \
                                                                    computed_success_wt, \
                                                                    computed_success_mut, \
                                                                    computed_additional_information, \
                                                                    computed_error, \
                                                                    computed_timeout, \
                                                                    computed_success)

computed_success_wt, computed_success_mut, computed_additional_information, \
computed_error, computed_timeout, computed_success = read_mutations('data/interfaceStability/v6/', \
                                                                    computed_success_wt, \
                                                                    computed_success_mut, \
                                                                    computed_additional_information, \
                                                                    computed_error, \
                                                                    computed_timeout, \
                                                                    computed_success)

computed_success_wt, computed_success_mut, computed_additional_information, \
computed_error, computed_timeout, computed_success = read_mutations('data/interfaceStability/v6.1/', \
                                                                    computed_success_wt, \
                                                                    computed_success_mut, \
                                                                    computed_additional_information, \
                                                                    computed_error, \
                                                                    computed_timeout, \
                                                                    computed_success)
                                                                    
computed_success_wt, computed_success_mut, computed_additional_information, \
computed_error, computed_timeout, computed_success = read_mutations('data/interfaceStability/v7/', \
                                                                    computed_success_wt, \
                                                                    computed_success_mut, \
                                                                    computed_additional_information, \
                                                                    computed_error, \
                                                                    computed_timeout, \
                                                                    computed_success)
                                                                    
#computed_success_wt, computed_success_mut, computed_additional_information, \
#computed_error, computed_timeout, computed_success = read_mutations('data/interfaceStability/v4/')
#
#computed_success_wt, computed_success_mut, computed_additional_information, \
#computed_error, computed_timeout, computed_success = read_mutations('data/interfaceStability-beagle07/v1/', \
#                                                                    computed_success_wt, \
#                                                                    computed_success_mut, \
#                                                                    computed_additional_information, \
#                                                                    computed_error, \
#                                                                    computed_timeout, \
#                                                                    computed_success)
#                                                                    
#computed_success_wt, computed_success_mut, computed_additional_information, \
#computed_error, computed_timeout, computed_success = read_mutations('data/interfaceStability-beagle07/v2/', \
#                                                                    computed_success_wt, \
#                                                                    computed_success_mut, \
#                                                                    computed_additional_information, \
#                                                                    computed_error, \
#                                                                    computed_timeout, \
#                                                                    computed_success)
#computed_success_wt, computed_success_mut, computed_additional_information, \
#computed_error, computed_timeout, computed_success = read_mutations('data/HAPMAP-2nd/', \
#                                                                    computed_success_wt, \
#                                                                    computed_success_mut, \
#                                                                    computed_additional_information, \
#                                                                    computed_error, \
#                                                                    computed_timeout, \
#                                                                    computed_success)

#computed_success_wt, computed_success_mut, computed_additional_information, \
#computed_error, computed_timeout, computed_success = read_mutations('data/OTHER/', \
#                                                                    computed_success_wt, \
#                                                                    computed_success_mut, \
#                                                                    computed_additional_information, \
#                                                                    computed_error, \
#                                                                    computed_timeout, \
#                                                                    computed_success)
#
#computed_success_wt, computed_success_mut, computed_additional_information, \
#computed_error, computed_timeout, computed_success = read_mutations('data/ALL-2nd/', \
#                                                                    computed_success_wt, \
#                                                                    computed_success_mut, \
#                                                                    computed_additional_information, \
#                                                                    computed_error, \
#                                                                    computed_timeout, \
#                                                                    computed_success)


## write the failed HAPMAPs
#with open('data_HAPMAP-failed.dat', 'w') as f:
#    for item in computed_error:
#        f.write('\t'.join(item.split('_')) + '\n')

res = open('disease_interfaceStability_predictions.dat', 'w')


res.write('ID' + '\t' + 'ddG_exp' + '\t' + 'ddG_calc' + '\t' + 'dG_wt\t' + 'dG_mut\t' + \
          'intraclashesEnergy1_wt\t' + 'intraclashesEnergy2_wt\t' + \
          'Backbone_Hbond_wt\t' + 'Sidechain_Hbond_wt\t' + \
          'Van_der_Waals_wt\t' + 'Electrostatics_wt\t' + \
          'Solvation_Polar_wt\t' + 'Solvation_Hydrophobic_wt\t' + \
          'Van_der_Waals_clashes_wt\t' + 'entropy_sidechain_wt\t' + \
          'entropy_mainchain_wt\t' + 'sloop_entropy_wt\t' + \
          'mloop_entropy_wt\t' + 'cis_bond_wt\t' + \
          'torsional_clash_wt\t' + 'backbone_clash_wt\t' + \
          'helix_dipole_wt\t' + 'water_bridge_wt\t' + \
          'disulfide_wt\t' + 'electrostatic_kon_wt\t' + \
          'partial_covalent_bonds_wt\t' + 'energy_Ionisation_wt\t' + \
          'Entropy_Complex_wt\t' + \
          'intraclashesEnergy1_mut\t' + 'intraclashesEnergy2_mut\t' + \
          'Backbone_Hbond_mut\t' + 'Sidechain_Hbond_mut\t' + \
          'Van_der_Waals_mut\t' + 'Electrostatics_mut\t' + \
          'Solvation_Polar_mut\t' + 'Solvation_Hydrophobic_mut\t' + \
          'Van_der_Waals_clashes_mut\t' + 'entropy_sidechain_mut\t' + \
          'entropy_mainchain_mut\t' + 'sloop_entropy_mut\t' + \
          'mloop_entropy_mut\t' + 'cis_bond_mut\t' + \
          'torsional_clash_mut\t' + 'backbone_clash_mut\t' + \
          'helix_dipole_mut\t' + 'water_bridge_mut\t' + \
          'disulfide_mut\t' + 'electrostatic_kon_mut\t' + \
          'partial_covalent_bonds_mut\t' + 'energy_Ionisation_mut\t' + \
          'Entropy_Complex_mut\t' + 'Number_of_Residues\t' + \
          'matrix_score\t' + \
          'if_hydrophobic' + '\t' + 'if_hydrophilic' + '\t' + 'if_total' + '\t' + \
          'pcv_salt_equal_wt' + '\t' + 'pcv_salt_opposite_wt' + '\t' + 'pcv_hbond_wt' + '\t' + 'pcv_vdW_wt' + '\t' + \
          'pcv_salt_equal_self_wt' + '\t' + 'pcv_salt_opposite_self_wt' + '\t' + 'pcv_hbond_self_wt' + '\t' + 'pcv_vdW_self_wt' + '\t' + \
          'pcv_salt_equal_mut' + '\t' + 'pcv_salt_opposite_mut' + '\t' + 'pcv_hbond_mut' + '\t' + 'pcv_vdW_mut' + '\t' + \
          'pcv_salt_equal_self_mut' + '\t' + 'pcv_salt_opposite_self_mut' + '\t' + 'pcv_hbond_self_mut' + '\t' + 'pcv_vdW_self_mut' + '\t' + \
          'secondary_structure_wt' + '\t' + 'solvent_accessibility_wt' + '\t' + \
          'secondary_structure_mut' + '\t' + 'solvent_accessibility_mut' + '\t' + \
          'sift_score' + '\tseq_id_avg' + '\tseq_id_chain1' + '\tseq_id_chain2' + \
          '\tnormDOPE' + '\tcore_or_interface' + '\tdisease_type' + \
          '\n')





sift_missing = set()
polyphen_missing = set()
secondary_structure_mut_count_all = 0

print len(computed_success_wt), len(set(computed_success_wt))

disease_dict = dict()

for key in computed_success_wt:
    
    line_result_wt = computed_success_wt[key]
    line_result_mut = computed_success_mut[key]
    line_additional_info = computed_additional_information[key]

    core_or_interface = line_additional_info[1]
    
    line, sift_missing_one, polyphen_missing_one, secondary_structure_mut_count = make_output(line_result_wt, line_result_mut, line_additional_info, polyphen, mutations)
    secondary_structure_mut_count_all += secondary_structure_mut_count
        
    for item in sift_missing_one:
        sift_missing.add(item.split('_')[0])
    for item in polyphen_missing_one:
        polyphen_missing.add(item)
    
    disease = line.split('\t')[-1].strip()
#    print 'disease', disease
    disease_dict.setdefault(disease, []).append(1)
    
    res.write(line)


with open('sift_missing.txt', 'w') as f:
    for line in sift_missing:
        f.write(line + '\n')
#        try:
#            uniprotKB, mutation = line.split('_')
#        except:
#            print line
#            raise
#        f.write(uniprotKB + '\t' + mutation + '\n')

#with open('polyphen_missing.txt', 'w') as f:
#    for line in polyphen_missing:
#        f.write(line + '\n')
#        uniprotKB, mutation = line.split('_')
#        fromAA, pos, toAA = mutation[0], mutation[1:-1], mutation[-1]
#        pph = uniprotKB + ' ' + pos + ' ' + fromAA + ' ' + toAA
#        f.write(pph + '\n')

with open('failed.dat', 'w') as f:
    for item in computed_error:
        uniprotKB, mutation = item.split('_')
        f.write(uniprotKB + '\t' + mutation + '\n')

computed_timeout_done = set()
#with open('timeout-1st.dat', 'r') as f:
#    for line in f:
#        if line == '' or line == '\n':
#            continue
#        uniprotKB, mutation = line.split('\t')
#        computed_timeout_done.add( (uniprotKB, mutation.strip()) )
with open('timeout.dat', 'w') as f:
    for item in computed_timeout:
        uniprotKB, mutation = item.split('_')
        if (uniprotKB, mutation) in computed_timeout_done:
            continue
        f.write(uniprotKB + '\t' + mutation + '\n')
with open('success.dat', 'w') as f:
    for item in computed_success:
        f.write(item + '\n')
#            uniprotKB, mutation = item.split('_')
#            f.write(uniprotKB + '\t' + mutation + '\n')


for key in disease_dict:
    print 'key', key, len(disease_dict[key])

print 'DSSP missing for', secondary_structure_mut_count_all


#done = ['failed.dat.running-v2', 'failed.dat.v2.01', 'failed.dat.v2.02']
success = list()
#for fname in done:
#    with open(fname, 'r') as f:
#        for l in f:
#            success.append(l)
#print 'Calculated at ccbr:', len(success)

#res = open('failed.dat', 'w')
#with open('failed.dat', 'r') as f:
#    for l in f:
#        if l in success:
#            pass
#        else:
#            res.write(l)
#res.close()
