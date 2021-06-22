import os
import glob
import dill
import numpy as np
from numpy import linalg as LA
import sys

folder = sys.argv[1]
os.chdir(folder)
level_list = os.listdir()
levels = []
for level in level_list:
    if os.path.isdir(os.path.join(level)):
        levels.append(level)

e = 0
g = 0

for level in levels:
    os.chdir(level)
    frags = glob.glob('*.dill')
    for i in frags:
        #undill and run e, g, hess, apt etc
        infile = open(i, 'rb')
        new_class = dill.load(infile)
        #print("Fragment ID:", level, i)
        e += new_class.energy
        g += new_class.grad
        #print(new_class.molecule.atomtable)
        outfile = open(i, "wb")
        dill.dump(new_class, outfile)
        infile.close()
        outfile.close()
    print("############## Energy for level(", level, ") is :", e, "################")
    os.chdir('../')

#os.chdir('../')
print(os.getcwd())
print("MIM Energy (Hartree):\n", e)
print("Norm of Gradient:\n", LA.norm(g))
#print("about to save energy and gradient")
np.save('energy.npy', e)
np.save('gradient.npy', g)
