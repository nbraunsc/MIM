import numpy as np
import psi4

#psi4.core.set_output_file('output.dat', False)

class Psi4():
    """ Psi4 Numpy quantum chemistry backend class

    An instance of this class is passed into the Fragment class
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
    
    def energy_gradient(self, input_xyz):
        string = "noreorient \n" + str(input_xyz).replace("[", " ").replace(",", "").replace("]", "\n").replace("'", "")
        mol = psi4.geometry(string)
        options = {'BASIS': self.basis}
        psi4.set_options(options)
        mol.set_molecular_charge(self.charge)
        mol.set_multiplicity(self.spin)
        #psi4.set_memory('100 GB')
        e = 0
        g = 0
        h = 0

        if self.theory == 'RHF':
            e, wfn = psi4.energy('scf', return_wfn=True)
            g = np.array(psi4.gradient('scf'))
            h = np.array(psi4.hessian('scf'))
        
        if self.theory == 'UHF':
            psi4.set_options({'reference': 'uhf', 'basis': self.basis})
            e, wfn = psi4.energy('scf', return_wfn=True)
            g = np.array(psi4.gradient('scf'))
            h = 0
        
        if self.theory == 'ROHF':
            psi4.set_options({'reference': 'rohf', 'basis': self.basis})
            e, wfn = psi4.energy('scf', return_wfn=True)
            g = np.array(psi4.gradient('scf'))
            h = 0
        
        if self.theory == 'MP2':
            e, wfn = psi4.energy('MP2', return_wfn=True)
            g = np.array(psi4.gradient('MP2'))
            h = 0
            #h = np.array(psi4.hessian('MP2'))
        
        if self.theory == 'CCSD':
            e = psi4.energy('CCSD')
            g = np.array(psi4.gradient('CCSD'))
            h = np.array(psi4.hessian('CCSD'))
        
        if self.theory == 'CCSD(T)':
            e = psi4.energy('CCSD(T)')
            g = np.array(psi4.gradient('CCSD(T)'))
            h = np.array(psi4.hessian('CCSD(T)'))
        
        if self.theory == 'CISD':
            e = psi4.energy('CISD')
            g = np.array(psi4.gradient('CISD'))
            h = np.array(psi4.hessian('CISD'))
        
        if self.theory == 'CISDT':
            e = psi4.energy('CISDT')
            g = np.array(psi4.gradient('CISDT'))
            h = np.array(psi4.hessian('CISDT'))
        
        return e, g, h
    
    def apply_field(self, E, input_xyz, com, origin, direction):
        string = "symmetry C1\n noreorient\n no_com \n" + str(input_xyz).replace("[", " ").replace(",", "").replace("]", "\n").replace("'", "")
        mol = psi4.geometry(string)
        options = {
                'basis' : self.basis,
                'PERTURB_H' : True,
                'PERTURB_WITH' : 'DIPOLE',
                'PERTURB_DIPOLE' : E,
                'save_jk' : True
        }
        psi4.set_options(options)
        mol.set_molecular_charge(self.charge)
        
        #setting everything to zero
        e = 0
        g = 0
        dip = 0
        nuc_dip = 0
        g_nuc = 0
        g_elec = 0
        
        if self.theory == 'RHF':
            e, scf_wfn = psi4.energy('scf', return_wfn=True)
            densitya = psi4.core.Matrix.to_array(scf_wfn.Da())
            densityb = psi4.core.Matrix.to_array(scf_wfn.Db())
            focka = psi4.core.Matrix.to_array(scf_wfn.Fa())
            fockb = psi4.core.Matrix.to_array(scf_wfn.Fb())
            hcore = psi4.core.Matrix.to_array(scf_wfn.H())   #core hamiltonian
            #print(hcore)

            mints = psi4.core.MintsHelper(scf_wfn.basisset())
            T = np.asarray(mints.ao_kinetic())
            V = np.asarray(mints.ao_potential())
            S = np.asarray(mints.ao_overlap())
            #print(focka + fockb)

            g = np.array(psi4.gradient('scf'))*1.88973    #H/B -> H/A
            #print(e)
            #print(g)
       
        if self.theory == 'UHF':
            psi4.set_options({'reference': 'uhf', 'basis': self.basis})
            e, scf_wfn = psi4.energy('scf', return_wfn=True)
            g = np.array(psi4.gradient('scf'))*1.88973    #H/B -> H/A
        
        if self.theory == 'ROHF':
            psi4.set_options({'reference': 'rohf', 'basis': self.basis})
            e, wfn = psi4.energy('scf', return_wfn=True)
            g = np.array(psi4.gradient('scf'))*1.88973
            
        if self.theory == 'MP2':
            e = psi4.energy('MP2')
            g = np.array(psi4.gradient('MP2'))*1.88973    #H/B -> H/A
        
        if self.theory == 'CCSD':
            e = psi4.energy('CCSD')
            g = np.array(psi4.gradient('CCSD'))*1.88973    #H/B -> H/A
        
        if self.theory == 'CCSD(T)':
            e = psi4.energy('CCSD(T)')
            g = np.array(psi4.gradient('CCSD(T)'))*1.88973    #H/B -> H/A
        
        if self.theory == 'CISD':
            e = psi4.energy('CISD')
            g = np.array(psi4.gradient('CISD'))*1.88973    #H/B -> H/A
        
        if self.theory == 'CISDT':
            e = psi4.energy('CISDT')
            g = np.array(psi4.gradient('CISDT'))*1.88973    #H/B -> H/A
        

        #nuc_dip = np.zeros((3))
        nuc_dip = np.array(mol.nuclear_dipole())
        
        #compute nuclear gradient
        Gradient = {}
        Gradient["N"] = psi4.core.Matrix.to_array(mol.nuclear_repulsion_energy_deriv1([0, 0, 0]))
        N_grad = psi4.core.Matrix.from_array(Gradient["N"])
        N_grad.name = "NUCLEAR GRADIENT"
        N_grad.print_out()
        g_nuc = np.array(N_grad)*1.88973 #H/B
        
        g_elec = 0
        print("nuclear repo energy:", mol.nuclear_repulsion_energy())
        return e, g, dip, nuc_dip, g_nuc, g_elec

    def get_dipole(self, coords_new):
        """ This only for RHF and is used when building the APT's from numerical diff w.r.t 
        atomic coordinates.

        There are no CCSD(T) dipole moments implemented in Psi4.

        """
        string = "noreorient \n" + str(coords_new).replace("[", " ").replace(",", "").replace("]", "\n").replace("'", "")
        mol = psi4.geometry(string)
        options = {'BASIS': self.basis}
        psi4.set_options(options)
        mol.set_molecular_charge(self.charge)
        #psi4.set_memory('100 GB')
        
        #calculate nuclear dip moment
        nuc_dip = np.array(mol.nuclear_dipole())
        dipole = np.zeros((3))  #Debye
        method = self.theory
        
        if self.theory == 'RHF':
            method = 'SCF'
            psi4.prop('scf', properties=["DIPOLE"])

        if self.theory == 'UHF':        #unsure if this one works
            method = 'SCF'
            psi4.set_options({'reference': 'uhf', 'basis': self.basis})
            psi4.prop('scf', properties=["DIPOLE"])
        
        if self.theory == 'ROHF':
            method = 'SCF'
            psi4.set_options({'reference': 'rohf', 'basis': self.basis})
            psi4.prop('scf', properties=["DIPOLE"])
        
        if self.theory == 'MP2':
            psi4.prop('mp2', properties=["DIPOLE"])
        
        if self.theory == 'CCSD':
            method = 'CC'
            e, scf_wfn = psi4.energy('scf', return_wfn=True)
            psi4.prop('ccsd', properties=["DIPOLE"], ref_wfn=scf_wfn)
            #psi4.prop('ccsd', properties=["DIPOLE"], return_wfn=False)
        
        if self.theory == 'CISD':
            method = 'CI'
            psi4.prop('cisd', properties=["DIPOLE"])
        
        if self.theory == 'CISDT':
            method = 'CI'
            psi4.prop('cisdt', properties=["DIPOLE"])
        
        #if self.theory == 'CC2':
        #    psi4.prop('cc2', properties=["DIPOLE"])

        dipole[0] = psi4.get_variable(method + " DIPOLE X")
        dipole[1] = psi4.get_variable(method + " DIPOLE Y")
        dipole[2] = psi4.get_variable(method + " DIPOLE Z")
        
        return dipole, nuc_dip
