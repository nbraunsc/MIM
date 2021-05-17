#PBS -l walltime=00:24:00:00
#PBS -l nodes=1:ppn=4
#PBS -l mem=50GB
#PBS -q nmayhall_lab
#PBS -A qcvt_doe
#PBS -W group_list=nmayhall_lab
#PBS -M nbraunsc@vt.edu
#PBS -m a

module purge
module load gcc/5.2.0
module load Anaconda/5.2.0

#need num of threads for python jobs, keep one until parallizing on single node
MKL_NUM_THREADS=1

source activate pyconda

cd $PBS_O_WORKDIR

echo $LEVEL
echo $BATCH

python run.py $LEVEL $BATCH $FOLDER

exit;


