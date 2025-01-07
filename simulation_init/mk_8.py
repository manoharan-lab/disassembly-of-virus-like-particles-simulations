#!/usr/bin/env python

import os
import string
from shutil import copy
from subprocess import Popen

batch_name = '8_PS_off'  # Organize all simulation runs from this batch into one folder
base_dir = './' + batch_name
if not os.path.exists(base_dir):
    os.mkdir(base_dir)

seed_file = './seed'
seed1 = 8140
seed2 = 15640
seed3 = 17290


# Simulation batch 1 (Debye 0.1 - 0.2, subunit-subunit attraction 1.6 - 2.4)

seed = seed1
open(seed_file, 'w').write(str(seed))  # Update seed file with the current seed

repetitions = list(range(30))  # How many times to run each parameter set
debye_lengths = [0.1, 0.15, 0.2]
subunit_attractions = [1.6, 1.8, 2.0, 2.2, 2.4]
PS_spacings = [10]

sh_tmp = string.Template(open('base.sh').read())
in_tmp = string.Template(open('base.txt').read())

for rep in repetitions:
    for deb in debye_lengths:
        for sa in subunit_attractions:
            for pss in PS_spacings:
                seed = seed + 1
                run_dir = '%s-%s_%s-%s_%s-%s_%s' % ("DEB", str(deb), "ATT", str(sa), "PSS", str(pss), str(rep))
                print(run_dir)
                if not os.path.exists(os.path.join(base_dir, run_dir)):
                    os.mkdir(os.path.join(base_dir, run_dir))
                job_name = batch_name + '_' + run_dir

                # Generate input and job submission script files for this run

                input_file = in_tmp.substitute({'SEED': str(seed), 'DEB': str(deb), 'EBOND': str(sa), 'PSS': str(pss)})
                open(os.path.join(base_dir, run_dir) + '/sd%s.input' % seed, 'w').write(input_file)
                sh = sh_tmp.substitute(SEED=str(seed), JNAME=job_name, NCPU='--mode=gpu',
                                       DIR=os.path.join(base_dir, run_dir))
                open(os.path.join(base_dir, run_dir) + '/submit.sh', 'w').write(sh)
                Popen('chmod +x ' + os.path.join(base_dir, run_dir) + '/submit.sh', shell=True).wait()

                # Copy other input files into the run directory (record-keeping)

                copy('./base.sh', os.path.join(base_dir, run_dir))
                copy('./wales.py', os.path.join(base_dir, run_dir))
                copy('./generate_PS.py', os.path.join(base_dir, run_dir))
                copy('./config_dict.py', os.path.join(base_dir, run_dir))
                copy('./assembled_capsid_8.gsd', os.path.join(base_dir, run_dir, 'assembled_capsid.gsd'))

                open(seed_file, 'w').write(str(seed))  # Update seed file with the current seed
                copy(seed_file, os.path.join(base_dir, run_dir))

                Popen('sbatch ' + os.path.join(base_dir, run_dir) + '/submit.sh', shell=True).wait()


#  Simulation batch 2 (Debye 0.05 and 0.25, subunit-subunit attraction 1.6 - 2.4)

seed = seed2
open(seed_file, 'w').write(str(seed))  # Update seed file with the current seed

repetitions = list(range(30))  # How many times to run each parameter set
debye_lengths = [0.05, 0.25]
subunit_attractions = [1.6, 1.8, 2.0, 2.2, 2.4]
PS_spacings = [10]

sh_tmp = string.Template(open('base.sh').read())
in_tmp = string.Template(open('base.txt').read())

for rep in repetitions:
    for deb in debye_lengths:
        for sa in subunit_attractions:
            for pss in PS_spacings:
                seed = seed + 1
                run_dir = '%s-%s_%s-%s_%s-%s_%s' % ("DEB", str(deb), "ATT", str(sa), "PSS", str(pss), str(rep))
                print(run_dir)
                if not os.path.exists(os.path.join(base_dir, run_dir)):
                    os.mkdir(os.path.join(base_dir, run_dir))
                job_name = batch_name + '_' + run_dir

                # Generate input and job submission script files for this run

                input_file = in_tmp.substitute({'SEED': str(seed), 'DEB': str(deb), 'EBOND': str(sa), 'PSS': str(pss)})
                open(os.path.join(base_dir, run_dir) + '/sd%s.input' % seed, 'w').write(input_file)
                sh = sh_tmp.substitute(SEED=str(seed), JNAME=job_name, NCPU='--mode=gpu',
                                       DIR=os.path.join(base_dir, run_dir))
                open(os.path.join(base_dir, run_dir) + '/submit.sh', 'w').write(sh)
                Popen('chmod +x ' + os.path.join(base_dir, run_dir) + '/submit.sh', shell=True).wait()

                # Copy other input files into the run directory (record-keeping)

                copy('./base.sh', os.path.join(base_dir, run_dir))
                copy('./wales.py', os.path.join(base_dir, run_dir))
                copy('./generate_PS.py', os.path.join(base_dir, run_dir))
                copy('./config_dict.py', os.path.join(base_dir, run_dir))
                copy('./assembled_capsid_8.gsd', os.path.join(base_dir, run_dir, 'assembled_capsid.gsd'))

                open(seed_file, 'w').write(str(seed))  # Update seed file with the current seed
                copy(seed_file, os.path.join(base_dir, run_dir))

                Popen('sbatch ' + os.path.join(base_dir, run_dir) + '/submit.sh', shell=True).wait()


# Simulation batch 3 (Debye 0.05 - 0.25, subunit-subunit attraction 2.6)
seed = seed3
open(seed_file, 'w').write(str(seed))  # Update seed file with the current seed

repetitions = list(range(30))  # How many times to run each parameter set
debye_lengths = [0.05, 0.10, 0.15, 0.20, 0.25]
subunit_attractions = [2.6]
PS_spacings = [10]

sh_tmp = string.Template(open('base.sh').read())
in_tmp = string.Template(open('base.txt').read())

for rep in repetitions:
    for deb in debye_lengths:
        for sa in subunit_attractions:
            for pss in PS_spacings:
                seed = seed + 1
                run_dir = '%s-%s_%s-%s_%s-%s_%s' % ("DEB", str(deb), "ATT", str(sa), "PSS", str(pss), str(rep))
                print(run_dir)
                if not os.path.exists(os.path.join(base_dir, run_dir)):
                    os.mkdir(os.path.join(base_dir, run_dir))
                job_name = batch_name + '_' + run_dir

                # Generate input and job submission script files for this run

                input_file = in_tmp.substitute({'SEED': str(seed), 'DEB': str(deb), 'EBOND': str(sa), 'PSS': str(pss)})
                open(os.path.join(base_dir, run_dir) + '/sd%s.input' % seed, 'w').write(input_file)
                sh = sh_tmp.substitute(SEED=str(seed), JNAME=job_name, NCPU='--mode=gpu',
                                       DIR=os.path.join(base_dir, run_dir))
                open(os.path.join(base_dir, run_dir) + '/submit.sh', 'w').write(sh)
                Popen('chmod +x ' + os.path.join(base_dir, run_dir) + '/submit.sh', shell=True).wait()

                # Copy other input files into the run directory (record-keeping)

                copy('./base.sh', os.path.join(base_dir, run_dir))
                copy('./wales.py', os.path.join(base_dir, run_dir))
                copy('./generate_PS.py', os.path.join(base_dir, run_dir))
                copy('./config_dict.py', os.path.join(base_dir, run_dir))
                copy('./assembled_capsid_8.gsd', os.path.join(base_dir, run_dir, 'assembled_capsid.gsd'))

                open(seed_file, 'w').write(str(seed))  # Update seed file with the current seed
                copy(seed_file, os.path.join(base_dir, run_dir))

                Popen('sbatch ' + os.path.join(base_dir, run_dir) + '/submit.sh', shell=True).wait()

