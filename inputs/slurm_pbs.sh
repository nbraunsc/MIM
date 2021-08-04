#!/bin/bash

###SBATCH -p normal_q
#SBATCH -p preemptable_q
#SBATCH -N 1  # this requests 1 node. 
#SBATCH --mem=2GB
##SBATCH --mem-per-cpu=10GB #memory requested for each core (or CPU)
#SBATCH -t 04:00:00
#SBATCH --account=nmayhall_group
##SBATCH --mail-user=nbraunsc@vt.edu
##SBATCH --mail-type=FAIL,ARRAY_TASKS

#Remake my environment if HOME variable not set
if [ -z ${HOME+x} ];
then
    export HOME=$(echo ~)
    source /etc/profile
    source /etc/bashrc
    source $HOME/.bashrc
fi

hostname

module reset
module load site/tinkercliffs-rome/easybuild/setup  
module load site/tinkercliffs/easybuild/setup  
module load Anaconda3/2020.07
module load gcc/8.2.0

source activate mim_env

cd $SLURM_SUBMIT_DIR

echo $LEVEL
echo $OUTFILE
export ERROR=${OUTFILE%%.*}.error
echo $ERROR

export TEMP=$LEVEL/"$OUTFILE.scr"
mkdir $TEMP
echo $TEMP

#Start an rsync command which runs in the background and keeps a local version of the output file up to date
#touch $LEVEL/$OUTFILE
#while true; do rsync -av $OUTFILE $LEVEL/"$OUTFILE"; sleep 60; done &

#change batch_list into bash array
export BATCH_LIST=$1
IFS=', ' read -r -a array <<< "$BATCH_LIST"
echo "${array[@]}"

export len=${#array[@]}
export len1=$(( $len + 1 ))

for i in "${array[@]}"
do
    echo ${i}
    export FILE=$LEVEL/"$i.out"
    python run.py $LEVEL $i $FOLDER $TEMP &>> $FILE &
done

#tail -1 $LEVEL/"$i.out"

###make empty status file (testing to see if this fixed fails)
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
