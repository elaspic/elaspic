#!/usr/bin/python

import os
import shutil
import fileinput
import subprocess

folder_list = ['/home/alexey/working/pipeline/code/results-' + str(i) + '/' for i in range(0,9)]
subfolder_list = ['alignments/', 'bestModels/', 'pdbFiles/']
output_folder = '/home/alexey/working/pipeline/code/results-all/'


# function for making folders if they don't already exist
def make_directory(directory_path):
    system_command = 'mkdir -p ' + directory_path
    childProcess = subprocess.Popen(system_command, stdout=subprocess.PIPE, stderr=subprocess.PIPE, shell=True)
    result, error = childProcess.communicate()
    if childProcess.returncode != 0:
        print 'Failed making directory', directory_path


for folder in folder_list:
    for subfolder in subfolder_list:
        make_directory(output_folder + subfolder)
        for file in os.listdir(folder + subfolder):
            shutil.copy2(folder + subfolder + file, output_folder + subfolder + file)
    with open(output_folder + 'result_additional_information.log', 'a') as fout:
        for line in fileinput.input(folder + 'result_additional_information.log'):
            fout.write(line)
    with open(output_folder + 'result_mut.log', 'a') as fout:
        for line in fileinput.input(folder + 'result_mut.log'):
            fout.write(line)
    with open(output_folder + 'result_wt.log', 'a') as fout:
        for line in fileinput.input(folder + 'result_wt.log'):
            fout.write(line)






