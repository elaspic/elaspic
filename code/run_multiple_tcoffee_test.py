from os import environ
import subprocess
import argparse

def call_tcoffee(alignInFile, tcoffee_path, GAPOPEN=-0.0, GAPEXTEND=-0.0):
    """
    calls t_coffee in expresso mode (make sure BLAST is installed locally)
    
    input:  alignInFile     type: string        file containing two sequences
                                                in fasta format
            tmpPATH         type: string        used to set a unique path
                                                for the t_coffee tmp, cache,
                                                and lock directory (for 
                                                paralellization)
            GAPOPEN         type: int or str    see t_coffee manual
            GAPEXTEND       type: int or str    see t_coffee manual
            
            
    return: Biopython multiple sequence alignment object
    """

    out = tmpPath + tcoffee_path + 'sequenceAlignment.aln'
    
    # try the alignment in expresso mode (structure based with sap alignment)
    system_command, my_env = call_tcoffee_system_command(alignInFile, out, 'expresso', tcoffee_path)
    childProcess = subprocess.Popen(system_command, 
                                    stdout=subprocess.PIPE, 
                                    stderr=subprocess.PIPE, 
                                    shell=True, 
                                    env=my_env,
                                    )
    result, error = childProcess.communicate()
    rc = childProcess.returncode
    print "tcoffee rc:", rc

def call_tcoffee_system_command(alignInFile, out, mode, tcoffee_path):
    # to be able to run parallel instances of T_Coffee the environment
    # variables have to be set to unique paths for each T_Coffee call
    # also make sure to change the directory to a unique one bevore
    # calling T_Coffee
    my_env = environ.copy()
#    print "my_env----", str('00'), "end\n"
    my_env['HOME_4_TCOFFEE'] = tmpPath + tcoffee_path
#    my_env['DIR_4_TCOFFEE'] = "/home/alexey/working/bin/tcoffee/Version_10.00.r1613"
    my_env['TMP_4_TCOFFEE'] = tmpPath  + tcoffee_path + 'tmp/'
    my_env['CACHE_4_TCOFFEE'] = tmpPath  + tcoffee_path + 'cache/'
    my_env['LOCKDIR_4_TCOFFEE'] = tmpPath  + tcoffee_path + 'lck/'
    my_env['ERRORFILE_4_TCOFFEE'] = tmpPath  + tcoffee_path + 't_coffee.ErrorReport'
    my_env['BLASTDB'] = '/home/alexey/working/bin/ncbi-blast-2.2.28+/pdbaa_db'
    if mode == 'expresso':
        system_command = 'cd ' + tmpPath + tcoffee_path + ' && t_coffee' + \
                         ' -mode expresso' + \
                         ' -seq ' + alignInFile + \
                         ' -blast_server=LOCAL' + \
                         ' -method clustalw_pair slow_pair' + \
                         ' -pdb_db=pdbaa -protein_db=uniprot' + \
                         ' -outorder=input' + \
                         ' -output fasta_aln' + \
                         ' -quiet -no_warning' + \
                         ' -outfile=' + out
    if mode == 't_coffee':
        system_command = 'cd ' + tmpPath + tcoffee_path + ' && t_coffee' + \
                         ' -mode expresso' + \
                         ' -method TMalign_pair' + \
                         ' -seq ' + alignInFile + \
                         ' -blast_server=LOCAL' + \
                         ' -method clustalw_pair slow_pair' + \
                         ' -pdb_db=pdbaa -protein_db=pdbaa' + \
                         ' -outorder=input' + \
                         ' -output fasta_aln' + \
                         ' -quiet -no_warning' + \
                         ' -outfile=' + out
    return system_command, my_env
    

parser = argparse.ArgumentParser()
parser.add_argument('uniqueID', help='unique identifier for t_coffee path')
arguments = parser.parse_args()

tmpPath = '/tmp/pipeline/'

tcoffee_path = 'tcoffee' + str(arguments.uniqueID) + '/'
print tcoffee_path


system_command = 'mkdir -p ' + tmpPath + tcoffee_path + \
                    ' && cp ' + tmpPath + 'seqfiles.fasta ' + tmpPath + tcoffee_path + 'seqfiles.fasta'
childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
result, error = childProcess.communicate()
if childProcess.returncode != 0:
    print "Couldn't create t_coffee folder"
call_tcoffee('seqfiles.fasta', tcoffee_path)


system_command = 'cd ' + tmpPath
childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
result, error = childProcess.communicate()