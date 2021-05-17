import mim
from mim import runpie, Molecule, fragmentation, Fragment, Pyscf, Psi4

import numpy as np
import os
import sys
import pickle
import dill
import shutil

input_file = sys.argv[1]
folder = sys.argv[3]

from input_file import *
coords = sys.argv[2]
path = os.getcwd() + "/" + coords

if software == 'Pyscf':
    software = Pyscf.Pyscf

if software == 'Psi4':
    software = Psi4.Psi4

#if software == "Molcas" or "OpenMolcas":
#    software = Molcas.Molcas

if mim_levels == None:
    raise Exception('MIM level not defined')
    
if not frag_type:
    raise Exception('Fragmentation type not defined')
    
if frag_deg == None:
    raise Exception('Fragmentation distance cutoff not defined')
    
if frag_deg_large == None:
    raise Exception('Second fragmentation distance cutoff not defined')

if not basis_set:
    raise Exception('Basis set not defined')
    
if mim_levels == 2 and not low_theory: 
    raise Exception('Low level of theory not defined')

if not high_theory:
    raise Exception('High level of theory not defined')
    
if not software:
    raise Exception('Software for calculations not defined')
    
if stepsize == None:
    raise Exception('Step size for finite difference not defined')


if opt == False:
    input_molecule = Molecule.Molecule(path, special_prim=atom, special_charge=charge, special_spin=spin)
    input_molecule.initalize_molecule()
    obj_list = []
    
    if mim_levels == 1:
        frag1 = fragmentation.Fragmentation(input_molecule)
        frag1.do_fragmentation(fragtype=frag_type, value=frag_deg)
        frag1.initalize_Frag_objects(theory=high_theory, basis=basis_set, qc_backend=software, xc=xc, step_size=stepsize, local_coeff=1)
        os.path.abspath(os.curdir)
        os.chdir(folder)
        level_list = os.listdir()
        for level in level_list:
            if os.path.isdir(os.path.join(level)):
                shutil.rmtree(level)
        os.mkdir('frag1')
        os.chdir('frag1')
        for i in range(0, len(frag1.frags)):
            filename = "fragment" + str(i) + ".dill"
            outfile = open(filename, "wb")
            print("writing object:", frag1.frags[i])
            dill.dump(frag1.frags[i], outfile)
            #status_name = filename + '.status'
            #status = open(status_name, "wb")
            #dill.dump(int(0), status)
            #status.close()
            outfile.close()
        obj_list.append(frag1)
        os.chdir('../')
        
    if mim_levels == 2:
        #""" MIM high theory, small fragments"""
        frag1 = fragmentation.Fragmentation(input_molecule)
        frag1.do_fragmentation(fragtype=frag_type, value=frag_deg)
        frag1.initalize_Frag_objects(theory=high_theory, basis=basis_set, qc_backend=software, xc=xc, step_size=stepsize, local_coeff=1)
        os.path.abspath(os.curdir)
        os.chdir(folder)
        level_list = os.listdir()
        for level in level_list:
            if os.path.isdir(os.path.join(level)):
                shutil.rmtree(level)
        os.mkdir('frag1')
        os.chdir('frag1')
        for i in range(0, len(frag1.frags)):
            filename = "fragment" + str(i) + ".dill"
            outfile = open(filename, "wb")
            dill.dump(frag1.frags[i], outfile)
            #status_name = filename.replace('.dill', '.status')
            #status = open(status_name, "wb")
            #dill.dump(int(0), status)
            #status.close()
            outfile.close()
        obj_list.append(frag1)
        os.chdir('../')
        
        #""" MIM low theory, small fragments"""
        frag2 = fragmentation.Fragmentation(input_molecule)
        frag2.do_fragmentation(fragtype=frag_type, value=frag_deg)
        frag2.initalize_Frag_objects(theory=low_theory, basis=basis_set, qc_backend=software, xc=xc, step_size=stepsize, local_coeff=-1)
        print("done with initalize low theory small fragments")
        os.mkdir('frag2')
        os.chdir('frag2')
        for i in range(0, len(frag2.frags)):
            filename = "fragment" + str(i) + ".dill"
            outfile = open(filename, "wb")
            dill.dump(frag2.frags[i], outfile)
            #status_name = filename.replace('.dill', '.status')
            #status = open(status_name, "wb")
            #dill.dump(int(0), status)
            #status.close()
            outfile.close()
        obj_list.append(frag2)
        os.chdir('../')
        
        #""" MIM low theory, large fragments (iniffloate system)"""
        frag3 = fragmentation.Fragmentation(input_molecule)
        print("starting large fragments")
        frag3.do_fragmentation(fragtype=frag_type, value=frag_deg_large)
        print(frag3.derivs)
        frag3.initalize_Frag_objects(theory=low_theory, basis=basis_set, qc_backend=software, xc=xc, step_size=stepsize, local_coeff=1)
        print("done with initalize")
        os.mkdir('frag3')
        os.chdir('frag3')
        for i in range(0, len(frag3.frags)):
            filename = "fragment" + str(i) + ".dill"
            outfile = open(filename, "wb")
            dill.dump(frag3.frags[i], outfile)
            outfile.close()
        obj_list.append(frag3)
        os.chdir('../')
    
    os.chdir('../')
    if queue == 'pbs':
        cmd = "python batch.py %s %s pbs.sh %s"%(str(batch_size), folder, queue)       ##For Newriver
    if queue == 'slurm':
        cmd = 'python batch.py %s %s slurm_pbs.sh %s'%(str(batch_size), folder, queue)         ##For TinkerCliffs/Huckleberry/Infer
    
    if queue == 'local':
        cmd = 'python batch.py %s %s run.py %s'%(str(batch_size), folder, queue)
    
    os.system(cmd)
    print(cmd)
    
if opt == True:
    cmd = 'python opt.py %s %s %s'%(path, folder, queue)
    print(cmd)
    os.system(cmd)
