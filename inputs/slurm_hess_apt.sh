#!/bin/bash
#SBATCH -p normal_q

###SBATCH -p preemptable_q
#SBATCH -N 1  # this requests 1 node
##SBATCH --mem=10GB
#SBATCH -t 6:00:00
#SBATCH --account=nmayhall_group
#SBATCH --exclusive # this requests exclusive access to node for interactive jobs

##SBATCH --mail-user=nbraunsc@vt.edu
##SBATCH --mail-type=FAIL
##SBATCH --export=ALL

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

export ERROR=${OUTFILE%%.*}.error
echo $ERROR

export TEMP=$LEVEL/"$OUTFILE.reap"
mkdir $TEMP
echo $TEMP

#python run_opt.py $LEVEL $BATCH $FOLDER

#change batch_list into bash array
export BATCH_LIST=$1
IFS=', ' read -r -a array <<< "$BATCH_LIST"
echo "${array[@]}"

export len=${#array[@]}
export len1=$(( $len + 1 ))
# echo $len

for i in "${array[@]}"
do
    echo ${i}
    export FILE=$LEVEL/"$i.out"
    python run_opt.py $LEVEL $i $FOLDER $TEMP &>> $FILE &
done

#check if #.status files = #jobs in batch
echo "Finished submitting python run.py. Now changing into" $TEMP

#make empty status file (testing to see if this fixed fails)
cd $TEMP
touch test.status
ls
echo $len
while [ $(ls -lR ./*.status | wc -l) != $len1 ]
do
    echo "Jobs Not Done Yet"
    if grep -q Killed $LEVEL/$ERROR; then
        sendmail nbraunsc@vt.edu < $LEVEL/$OUTFILE
        echo "Job Killed, email sent"
        kill $$
        fi
    sleep 5
done

echo "Finished!"

#mkdir $LEVEL/"$OUTFILE.scr"
#cp * $LEVEL/"$OUTFILE.scr"
rm -r *

exit;
