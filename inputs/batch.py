import os
import sys
import glob
import mim
from mim import *
#from mim import runpie, Molecule, fragmentation, Fragment, Pyscf

batch_size = int(sys.argv[1])
folder = sys.argv[2]
script = sys.argv[3]
queue = sys.argv[4]

def batch(iterable, n=1):
    l = len(iterable)
    for ndx in range(0, l, n):
        yield iterable[ndx:min(ndx + n, l)]

os.path.abspath(os.curdir)
os.chdir(folder)
command_list = []

dirlist = []
for filename in os.listdir():
    if os.path.isdir(os.path.join(filename)):
        dirlist.append(filename)


for i in dirlist:
    os.chdir(i)
    files_dill = glob.glob('*.dill')
    if batch_size == None:
        cmd = 'python %s %s %s %s >> run_error.txt'%(script, path, string_num, folder)
        os.chdir('../../')
        os.system(cmd)
    else:
        for x in batch(files_dill, batch_size):
            path = os.getcwd()
            string_num = str()
            for frag in x:
                y = frag.replace("fragment", "").replace(".dill", "")
                string_num+=y+"_"
            submit_name = i + "_" + string_num
            os.chdir('../../')
            if queue == 'local':
                cmd = 'python %s %s %s %s >> run_error.txt'%(script, path, string_num, folder)
            if queue =='slurm':
                cmd = 'sbatch -e %s -J %s -o "%s" --export=LEVEL="%s",BATCH="%s",FOLDER="%s"  %s'%(path+"/"+submit_name+".error", submit_name, path+"/"+submit_name+".log", path, string_num, folder, script)     ##For TinkerCliffs/Huckleberry
            if queue == 'pbs':
                cmd = 'qsub -e %s -N %s -o %s -v LEVEL="%s",BATCH="%s",FOLDER="%s" %s'%(path, submit_name, path, path, string_num, folder, script)     ##For Newriver
            
            command_list.append(cmd)
            os.chdir(folder)
            os.chdir(i)
    os.chdir('../')
os.chdir('../')

for command in command_list:
    print("\n Submitting job:", command, "\n")
    os.system(command)

