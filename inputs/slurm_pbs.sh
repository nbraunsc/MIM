#!/bin/bash

#SBATCH -p normal_q
#SBATCH -N 1  # this requests 1 node, 1 core. 
#SBATCH --mem=20G
#SBATCH -t 01:00:00
#SBATCH --account=nmayhall_group
#SBATCH --mail-user=nbraunsc@vt.edu
#SBATCH --mail-type=FAIL

sleep 10
hostname

module reset
module load site/tinkercliffs-rome/easybuild/setup  
module load site/tinkercliffs/easybuild/setup  
module load Anaconda3/2020.07
module load gcc/8.2.0

source activate pyconda

cd $SLURM_SUBMIT_DIR

echo $LEVEL
echo $BATCH

python run.py $LEVEL $BATCH $FOLDER

exit;
