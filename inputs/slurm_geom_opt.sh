#!/bin/bash

#SBATCH -p normal_q
#SBATCH -N 1  # this requests 1 node, 1 core. 
#SBATCH --mem=10G
#SBATCH -t 01:00:00
#SBATCH --account=nmayhall_group
#SBATCH --mail-user=nbraunsc@vt.edu
#SBATCH --mail-type=FAIL

## SBATCH --exclusive # this requests exclusive access to node for interactive jobs

sleep 20
hostname

module reset
module load site/tinkercliffs-rome/easybuild/setup  #only for infer
module load site/tinkercliffs/easybuild/setup  #only for infer
module load Anaconda3/2020.07
module load gcc/8.2.0

source activate pyconda

cd $SLURM_SUBMIT_DIR

finished=0

while [ "$finished" != "1" ]
do
    finished=$(python checker.py $FOLDER)
    sleep 10
done

echo "Calculations are done!, end status:"
echo $finished

python eg_reap.py $FOLDER
echo "Reap is done!"

exit;
