#Input for MIM
import sys

#number of mim levels
mim_levels = 1

#can be 'distance' or 'graphical'
#frag_type = 'graphical'
frag_type = 'distance' 

#smaller fragmentation level
frag_deg = 10

#larger fragmentation level
frag_deg_large = 2.6

#Basis set for quantum calculation
#basis_set_high = 'ccpvdz'
basis_set_high = 'sto-3g'


#Always need to define high_theory
high_theory = 'RHF'
#high_theory = 'MP2'
#high_theory = 'DFT'

#Only define low_theory if mim_levels = 2
low_theory = 'RHF'
#low_theory = 'DFT'

#exchange-correlation functional for DFT
xc = 'LDA'

#could be Pyscf or Psi4, and eventually Qchem, or Molcas
#software = 'Psi4'  
software = 'Pyscf'  

#for second derivative by finite difference
stepsize = 0.001        

#batch_size for running calculations
batch_size = 8

#geometry optimization set to True or False
opt = False
#opt = True

#special params for primitive
atom = 'Si'
charge = -1
spin = 2 #2S+1

#pbs/slurm/Local?
queue = 'slurm'
#queue = 'local'
