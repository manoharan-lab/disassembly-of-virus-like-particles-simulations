#!/usr/bin/env python

import os
import string
from shutil import copy
from subprocess import Popen

batch_name = 'submission_scripts'  # Organize all simulation runs from this batch into one folder
base_dir = './' + batch_name
if not os.path.exists(base_dir):
    os.mkdir(base_dir)
if not os.path.exists('./raw'):
    os.mkdir('./raw')
if not os.path.exists('./plots'):
    os.mkdir('./plots')
if not os.path.exists('./xscale'):
    os.mkdir('./xscale')

# Batch parameters
repetitions = list(range(30))  # How many times to run each parameter set
debye_lengths = [0.05, 0.10, 0.15, 0.2, 0.25]
subunit_attractions = [1.6, 1.8, 2.0, 2.2, 2.4, 2.6]
PS_spacings = [10]

sh_tmp = string.Template(open('base.sh').read())

for rep in repetitions:
    for deb in debye_lengths:
        for sa in subunit_attractions:
            for pss in PS_spacings:
                run_name = '%s-%s_%s-%s_%s-%s_%s' % ("DEB", str(deb), "ATT", str(sa), "PSS", str(pss), str(rep))
                print(run_name)

                # Generate input and job submission script files for this run

                sh = sh_tmp.substitute(DEB=str(deb), ATT=str(sa), PSS=str(pss), REP=str(rep))
                open(base_dir + '/submit_'+ run_name + '.sh', 'w').write(sh)
                Popen('chmod +x ' + base_dir + '/submit_'+ run_name + '.sh', shell=True).wait()
                Popen('sbatch ' + base_dir + '/submit_'+ run_name + '.sh', shell=True).wait()

