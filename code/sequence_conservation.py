# -*- coding: utf-8 -*-

import subprocess


system_command = 'provean ' + \
    '-q examples/P04637.fasta ' + \
    '-v deleteme.var ' + \
    '-d /home/kimlab1/strokach/ncbi-blast-2.2.28+/db/nr ' + \
    '--psiblast `which psiblast` ' + \
    '--blastdbcmd `which blastdbcmd` ' + \
    '--cdhit `which cd-hit` ' + \
    '--save_supporting_set supset.fasta'

childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
result, error = childProcess.communicate()
if childProcess.returncode != 0:
    raise Exception('PROVEAN exited with an error:\n %s' % error)
                
### Results look something like this:
#[23:28:34] clustering subject sequences...
#[23:28:34] selecting clusters...
#[23:28:34] 0 subject sequences in 0 clusters were selected for supporting sequences.
#[23:28:34] use the query itself as a supporting sequence
#[23:28:34] loading subject sequences from a FASTA file...
#[23:28:34] scores were computed based on the query sequence itself.
## Number of clusters:	1
## Number of supporting sequences used:	1
#[23:28:34] computing delta alignment scores...
#[23:28:34] printing PROVEAN scores...
### PROVEAN scores ##
## VARIATION	SCORE
#M1A	-6.000

variations_started = False
for line in result:
    if 'Number of supporting sequences used:' in line:
        provean_supporting_set_length = int(line.split()[-1])
    if 'VARIATION\tSCORE' in line:
        variations_started = True
    if variations_started:
        provean_mutation, provean_score = line.split()

