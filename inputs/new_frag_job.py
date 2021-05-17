import dill
import numpy as np
import os
import sys

fragment_name = sys.argv[1]
level = sys.argv[2]
batch = sys.argv[3]
folder = sys.argv[4]

infile = open(fragment_name, 'rb')
frag_class = dill.load(infile)

#make changes as needed to frag_class

# example:
# frag_class.qc_backend.spin = 2 


cmd = 'sbatch -J %s -o "%s" --export=LEVEL="%s",BATCH="%s",FOLDER="%s"  slurm_pbs.sh'%(fragment_name+"rerun", path+"/"+submit_name+".out", path, string_num, folder)     ##For TinkerCliffs/Huckleberry

os.system(cmd)


