# -*- coding: utf-8 -*-

import argparse

parser = argparse.ArgumentParser()
parser.add_argument('log_filename', help='filename of the log file to be split')
parser.add_argument('number_of_consumers', type=int, default=8, help='number of consumers')
arguments = parser.parse_args()


ih = open(arguments.log_filename, 'r')
ohs = [open(arguments.log_filename.split('.')[0] + '-' + str(i) + '.log', 'w') for i in range(1, arguments.number_of_consumers + 1)]

last_index = 0
for line in ih:
    log_values = line.split(' - ')
    if len(log_values) > 1 and log_values[1].split('-')[0] == 'Consumer':
        last_index = int(log_values[1].split('-')[-1]) - 1 # index starts from 0, consumers start from 1
        ohs[last_index].writelines(line)
    else:
        ohs[last_index].writelines(line)
        
ih.close
(oh.close() for oh in ohs)
