import pickle
import dill
import numpy as np
import os
import sys

directory = sys.argv[1] #should be the directory of mim level
batch = sys.argv[2].split("_")
folder = sys.argv[3]
batch_list = []
for i in batch:
    if not i:
        continue
    else:
        entry = "fragment" + i + ".dill"
        batch_list.append(entry)

print("Batch list:", batch_list)
os.chdir(folder)
os.chdir(directory)
for i in batch_list:
    #undill and run e, g, hess, apt etc
    infile = open(i, 'rb')
    new_class = dill.load(infile)
    e, g, h = new_class.qc_backend()
    new_class.hess_apt(h)  #this should only be done at optimized geometry

    print("##################################this is from run_opt", e, g.shape)
    infile.close()
    
    ##redill with updated fragment e, g, hess, apt, etc
    outfile = open(i, "wb")
    dill.dump(new_class, outfile)
    outfile.close()
