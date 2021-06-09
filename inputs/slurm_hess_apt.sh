#!/bin/bash
#SBATCH -p normal_q
#SBATCH -N 1  # this requests 1 node, 1 core. 
#SBATCH --mem=20G
#SBATCH -t 05:00:00
### SBATCH --cpus-per-task=22
#SBATCH --account=nmayhall_group
#SBATCH --mail-user=nbraunsc@vt.edu
#SBATCH --mail-type=FAIL
#SBATCH --export=ALL
# SBATCH --get-user-env=30

#Remake my environment if HOME variable not set
if [ -z ${HOME+x} ];
then
    export HOME=$(echo ~)
    source /etc/profile
    source /etc/bashrc
    source $HOME/.bashrc
fi

# sleep 10
hostname

module reset
module load site/tinkercliffs-rome/easybuild/setup  #only for infer
module load site/tinkercliffs/easybuild/setup  #only for infer
module load Anaconda3/2020.07
module load gcc/8.2.0

source activate mim_env

cd $SLURM_SUBMIT_DIR

echo $LEVEL
echo $BATCH

#for i in $BATCH_LIST
#do
#    python run_opt.py $LEVEL $BATCH $FOLDER >& run_opt{$i}.out &
#done

python run_opt.py $LEVEL $BATCH $FOLDER

exit;
