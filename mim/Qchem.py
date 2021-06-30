import numpy as np
import os
import sys
from pyqchem import QchemInput
from pyqchem import get_output_from_qchem
from pyqchem.parsers.basic import basic_parser_qchem

class Qchem():
    """ QChem 5 quantum chemistry backend class.

    This will write an input file and submit it to run Qchem.  Then it will process the 
    output from QChem into numpy arrays and give back to the Fragment() class.
    """

    def __init__(self, theory=None, basis=None, tol=None, active_space=None, nelec=None,  nelec_alpha=None, nelec_beta=None, max_memory=None, xc=None, charge=0, spin=0):
        self.theory = theory
        self.basis = basis
        self.spin = spin
        self.tol = tol
        self.active_space = active_space    #number of orbitals in active space
        self.nelec = nelec  #number of electrons in the active space
        self.nelec_alpha = nelec_alpha
        self.nelec_beta = nelec_beta
        self.max_memory = max_memory
        self.xc = xc
        self.charge = charge

    def energy_gradient(self, new_coords):
        #Defining molecule
        molecule = Structure(coordinates=[[0.0, 0.0, 0.0],
                                          [0.0, 0.0, 0.9],
                                          [0.0, 0.0, -0.9]],
                             symbols=['O', 'H', 'H'],
                             charge=0,
                             multiplicity=1)
        
        #Preparing QCHEM input
        qc_input1 = QchemInput(molecule,
                              jobtype='force',
                              exchange='hf',
                              basis='6-31G')
        

    def energy_gradient(self, new_coords):
        self.writexyz(new_coords)
        lines ="""$molecule
READ '{filename}'
$end

$rem
JOBTYPE FORCE
METHOD '{theory}'
UNRESTRICTED FALSE
BASIS '{basis}'
$end

@@@

$molecule
READ '{filename}'
$end

$rem
JOBTYPE     freq
METHOD      '{theory}'
BASIS       '{basis}'
$end
"""
        context = {
        "filename":path/to/molecule.xyz,
        "charge":self.charge
        "spin":self.spin
        "theory":self.theory
        "basis":self.basis
        }

        with open(filename, "w") as myfile:
            myfile.write(lines.format(**context))
    
        myfile.close()

        cmd = 'command to submit qchem job'
        os.system(cmd)
        #read output file to get energy and gradient
        with open(filename) as f:
            if line.startswith("Total energy"):
                print(line)
                #e = pull out energy
                print("energy = ", e)
            
            if 'Gradient of scf energy' in f.read():
                print(f.read())
                #do something to save gradient as np.array
            if line.startswith(" Hessian of the SCF Energy"):
                #h = pull out hessian
                print("hessian =", h)
        outfile.close()
        return e, g

    def apply_field(self, E, input_xyz, com, origin, direction):
        """ Apply electric field in a cartesian direction.
        
        Parameters
        ----------
        E : np array
            This is a 1D array of an x, y, z.  Put magintude of wanted E field in the position of the.  Right now it is only set up for RHF.

            direction wanted.
        input_xyz : list
            Coordinates for molecule
        Returns
        -------
        e : float
            Energy
        g : ndarray
            gradient
        dipole : ndarray
            dipole moment
        """
        pass
        return e, g, dipole, nuc_dip, g_nuc, g_elec
    
    def get_dipole(self, coords_new):
        nuc_dip = []
        with open(filename) as f:
            if '    Dipole Moment (Debye)' in f.read():
                print(f.read())
                #dipole moment in line under this found line
                #do something to save dipole moment as np.array
        return dip, nuc_dip

        



        
        
