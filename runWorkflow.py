from CTDopts.CTDopts import _InFile, CTDModel, args_from_file
import sys
import os
import subprocess
import re

pattern = re.compile('Q\w{4}[0-9]{3}[a-zA-Z]\w')

wf_dir = sys.argv[1]
ctd_params = args_from_file(wf_dir + '/WORKFLOW-CTD')
ctd_files = args_from_file(wf_dir + '/IN-FILESTOSTAGE')

command = 'python IndividualizedProteins_1.0.py '

data_path = '%s/data/' % wf_dir
result_path = '%s/result/' % wf_dir
log_path = '%s/logs/' % wf_dir

for key in ctd_files.keys():
    if ctd_files[key] != '':
        fileName = ctd_files[key].split('/')[-1]
        ctd_files[key] = '%s%s' % (data_path, fileName)

command += '-s %s' % ctd_files['Somatic Mutations']

foundPatterns = pattern.findall(ctd_files['Somatic Mutations'])
if len(foundPatterns) > 1:
    identifier = pattern.findall(ctd_files['Somatic Mutations'])[0]
else:
    identifier = ctd_files['Somatic Mutations'].split('/')[-1].split('.')[0]

logfilename = 'personalizedProteins_%s.logs' % identifier
logfile = open(logfilename, 'w')

for param in ctd_params.keys():
        if param == 'b':
        	command += ' -%s %s' % ('b', result_path)
        elif param == 'g':
            if ctd_params[param] == 'true':
                command += ' -g %s' % ctd_files['Germline Mutations']
        elif param == 'd':
            command += ' -d %s' % identifier
        elif ctd_params[param] == 'true':
                command += ' -%s' % param
        else:
                command += ' -%s %s' % (param, ctd_params[param])

subprocess.call(command.split(),stderr=logfile, stdout=logfile)

logfile.close()
subprocess.call(["mv", logfilename, log_path])