#!/bin/bash
#SBATCH -J analysis
#SBATCH -p shared
#SBATCH -n 1
#SBATCH -o analysis.out
#SBATCH -e analysis.err
#SBATCH -t 0-01:00
#SBATCH --mem=4000

module load python/3.7.7-fasrc01
source activate hoomd_analysis
python ./intercapsomer_bonds.py $DEB $ATT $PSS $REP
python ./cluster_analysis.py $DEB $ATT $PSS $REP
python ./polymer_path.py $DEB $ATT $PSS $REP
python ./polymer_analysis2.py $DEB $ATT $PSS $REP
exit

