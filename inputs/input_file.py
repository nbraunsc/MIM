#Input for MIM
import sys

#number of mim levels
mim_levels = 1

#can be 'distance' or 'graphical'
frag_type = 'distance' 

#smaller fragmentation level
frag_deg = 1.6

#larger fragmentation level
frag_deg_large = 4.4

#Basis set for quantum calculation
basis_set = '6-31g'

#exchange-correlation functional for DFT
xc = 'pbe'

#Always need to define high_theory
high_theory = 'DFT'
#high_theory = 'CCSD'

#Only define low_theory if mim_levels = 2
low_theory = 'RHF'

#could be Pyscf or Psi4, and eventually Qchem, or Molcas
#software = 'Psi4'  
software = 'Pyscf'  

#for second derivative by finite difference
stepsize = 0.001        

#batch_size for running calculations
batch_size = 10

#geometry optimization set to True or False
opt = False
#opt = True

#special params for primitive
atom = 'N'
charge = -1
spin = 3 #2S+1

#pbs/slurm/Local?
#queue = 'slurm'
queue = 'local'
