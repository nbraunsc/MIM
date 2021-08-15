#!/bin/bash

##SBATCH -p preemptable_q
#SBATCH -p normal_q
#SBATCH -N 1  # this requests 1 node, 1 core. 
#SBATCH --mem=2GB
#SBATCH -t 04:00:00
#SBATCH --account=nmayhall_group
#SBATCH --mail-user=nbraunsc@vt.edu
#SBATCH --mail-type=FAIL,END
# SBATCH --get-user-env=30

## SBATCH --exclusive # this requests exclusive access to node for interactive jobs

# . /etc/bashrc


if [ -z ${HOME+x} ];
then
    export HOME=$(echo ~)
    source /etc/profile
    source /etc/bashrc
    source $HOME/.bashrc
fi

hostname
# sleep 10

module reset
module load site/tinkercliffs-rome/easybuild/setup  #only for infer
module load site/tinkercliffs/easybuild/setup  #only for infer
module load Anaconda3/2020.07
module load gcc/8.2.0

source activate mim_env

cd $SLURM_SUBMIT_DIR

finished=0

while [ "$finished" != "1" ]
do
    finished=$(python checker.py $FOLDER)
    sleep 10
done

echo "Calculations are done!, end status:"
echo $finished

python eg_reap.py $FOLDER >> eg_reap.out
echo "Reap is done!"

exit;
