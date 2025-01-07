#!/bin/bash

#SBATCH -J $JNAME
#SBATCH -p gpu_requeue
#SBATCH -o $DIR/$JNAME.txt
#SBATCH -e $DIR/$JNAME.err
#SBATCH --gres=gpu:1
#SBATCH --constraint=cc6.0
#SBATCH --mem=8000
#SBATCH -t 1-00:00
#SBATCH --open-mode=append

module load gcc/8.2.0-fasrc01 openmpi/4.0.1-fasrc01 cuda/10.2.89-fasrc01 python/3.7.7-fasrc01

cd $DIR

source activate hoomd_analysis

python ./generate_PS.py sd$SEED.input assembled_capsid.gsd

source activate hoomd-blue

echo "Start $SEED @ `date`"
python ./wales.py sd$SEED.input $SEED $NCPU
echo "End $SEED @ `date`"

exit
