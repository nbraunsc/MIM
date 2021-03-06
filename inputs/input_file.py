#Input for MIM
import sys

"""MIM Parameters"""
mim_levels = 1
frag_type = 'distance'          #distance or graphical
frag_deg = 2.6                  #smaller fragmentation level
frag_deg_large = 2.6            #larger fragmentation level

"""Quantum Chemsitry Parameters"""
software = 'Pyscf'  
high_theory = 'DFT'
low_theory = None
xc_high = 'LDA'
xc_low = None
basis_set_high = '631g*'
basis_set_low = None
opt = True
stepsize = 0.001                #for second derivative by finite difference

"""Special Parameters for Defects"""
atom = 'N'
charge = -1
spin = 3                        #2S+1

"""Run Parameters"""
batch_size = 64
queue = 'slurm'
