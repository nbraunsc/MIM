import numpy as np
import pyscf
from pyscf import gto, scf, ci, cc, mp, hessian, lib, grad, mcscf, ao2mo, tools, dft
from pyscf.prop.freq import rhf
from pyscf.prop.polarizability import rhf
#from pyscf.grad.rhf import GradientsBasics
from pyscf.geomopt.berny_solver import optimize
from berny import Berny, geomlib
#from pyscf.data import nist    #for unit conversions

class Pyscf():
    """ Pyscf quantum chemistry backend class
    
    An instance of this class is passed into the Fragment class
    """

    def __init__(self, theory=None, basis=None, tol=None, active_space=None, nelec=None,  nelec_alpha=None, nelec_beta=None, max_memory=None, xc=None, charge=0, spin=0):
        self.theory = theory
        self.basis = basis
        self.spin = spin      #pyscf wants 2S not 2S+1, S = nalpha-nbeta
        #self.spin = spin-1      #pyscf wants 2S not 2S+1, S = nalpha-nbeta
        self.charge = charge
        self.tol = tol
        self.active_space = active_space    #number of orbitals in active space
        self.nelec = nelec  #number of electrons in the active space
        self.nelec_alpha = nelec_alpha
        self.nelec_beta = nelec_beta
        self.max_memory = max_memory
        self.xc = xc
        #self.input_xyz = []

    def energy_gradient(self, input_xyz):
        mol = gto.Mole()
        mol.atom = input_xyz
        mol.basis = self.basis
        mol.charge = self.charge
        if self.spin != 0:
            mol.spin = self.spin-1
        else:
            mol.spin = self.spin
        print("spin:", mol.spin)
        print("charge:", mol.charge)
        #mol.unit = 'Bohr'
        mol.unit = 'Angstrom'
        #mol.symmetry = True
        mol.build()
        #mf.conv_tol = 1e-14
        
        if self.theory == 'DFT':
            if self.spin != 0:
                mf = dft.UKS(mol)
                mf.xc = self.xc
                e = mf.kernel()
                g = mf.nuc_grad_method().kernel()
                hess = mf.Hessian().kernel()
                h = hess.transpose(0,2,1,3).reshape(len(input_xyz)*3, len(input_xyz)*3, order='C')

            else:
                mf = dft.RKS(mol)
                mf.xc = self.xc
                e = mf.kernel()
                g = mf.nuc_grad_method().kernel()
                hess = mf.Hessian().kernel()
                h = hess.transpose(0,2,1,3).reshape(len(input_xyz)*3, len(input_xyz)*3, order='C')
                return e, g, h


        if self.theory == 'RHF': #Restricted HF calc
            mf = scf.RHF(mol).run()
            e = mf.kernel()
            mf.dump_scf_summary()
            g = mf.nuc_grad_method().kernel()
            hess = mf.Hessian().kernel()
            h = hess.transpose(0,2,1,3).reshape(len(input_xyz)*3, len(input_xyz)*3, order='C')
            return e, g, h
        
        if self.theory == 'UHF': #Restricted HF calc
            mf = scf.UHF(mol).run()
            e = mf.kernel()
            mf.dump_scf_summary()
            g = mf.nuc_grad_method().kernel()
            hess = mf.Hessian().kernel()
            h = hess.transpose(0,2,1,3).reshape(len(input_xyz)*3, len(input_xyz)*3, order='C')
            return e, g, h
        
        if self.theory == 'ROHF': #Restricted Open shell HF calc
            mf = scf.ROHF(mol).run()
            e = mf.kernel()
            mf.dump_scf_summary()
            g = mf.nuc_grad_method().kernel()
            h = 0
            #hess = mf.Hessian().kernel()
            #h = hess.transpose(0,2,1,3).reshape(len(input_xyz)*3, len(input_xyz)*3, order='C')
            return e, g, h
    
        if self.theory == 'MP2': #Perturbation second order calc
            mf = scf.RHF(mol).run()
            postmf = mp.MP2(mf).run()
            e = mf.kernel() + postmf.kernel()[0]
            g2 = postmf.nuc_grad_method()
            g = g2.kernel()
            h = 0
            return e, g, h
    
        if self.theory == 'CISD':    #CI for singles and double excitations
            mf = scf.RHF(mol).run()
            postmf = ci.CISD(mf).run()
            e = postmf.kernel()
            g2 = postmf.nuc_grad_method()
            g = g2.kernel()
            h = 0
            return e, g, h
    
        if self.theory == 'CCSD':    #Couple Cluster for singles and doubles
            mf = scf.RHF(mol).run()
            postmf = cc.CCSD(mf).run()
            e = mf.kernel() + postmf.kernel()[0]
            g2 = postmf.nuc_grad_method()
            g = g2.kernel()
            h = 0
            return e, g, h
        
        if self.theory == 'CCSD(T)':    #Couple Cluster for singles and doubles, pertub triples
            mf = scf.RHF(mol).run()
            mycc = cc.CCSD(mf).run()
            et = mycc.ccsd_t()
            #postmf = cc.ccsd_t(mf).run()
            #e = mf.kernel() + postmf.kernel()[0]
            e = mycc.e_tot + et
            g=0
            #g2 = postmf.nuc_grad_method()
            #g = g2.kernel()
            h = 0
            return e, g, h

        if self.theory == 'CASSCF':
            mf = scf.RHF(mol).run()
            postmf = mcscf.CASSCF(mf, self.active_space, self.nelec).run()
            e = postmf.kernel()
            g2 = postmf.nuc_grad_method()
            g = g2.kernel()
            h = 0
            return e, g, h

        if self.theory == 'CASCI':
            mf = scf.RHF(mol).run()
            h = 0
            return e, g, h
    
    def apply_field(self, E, input_xyz, com, origin, direction):
        """ This will apply an e lectric field in a specific direction for pyscf. gives E vector
        to make a new hcore. Used when building the derivative of gradient w.r.t electric field to
        compute the IR intensities.
    
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
        e=0
        g=0
        dipole1=0
        dm0 = [None]
        mol1 = gto.Mole()
        mol1.atom = input_xyz
        #mol1.atom = mol1.format_atom(input_xyz, origin=(com), unit='Ang')
        mol1.unit = 'Angstrom'
        mol1.charge = self.charge
        mol1.symmetry = False
        #charge_center =  np.array([0,  0, -0.0001248])
        #new_axes = np.identity(3)
        #mol1.groupname = 'C1'
        #mol1.unit = 'Bohr'
        #mol1.atom = mol1.format_atom(input_xyz, origin=(charge_center), axes=new_axes)
        mol1.basis = self.basis

        mol1.build(dump_input=True)
        print(input_xyz)
        print(mol1.atom_coords())
         
        #print("coords before field:\n", mol1.atom_coords()/1.88973) 

        dipole = 0
        g_nuc = 0
        g_elec = 0

        #Compute nuclear dipole moment
        charges = mol1.atom_charges()
        r = mol1.atom_coords()
        nucl_dip = np.einsum('i,ix->x', charges, r) 
        nuc_charge_center = nucl_dip/charges.sum()
        print(nuc_charge_center)
        #mol1.set_common_orig(nuc_charge_center)  # The gauge origin for dipole integral (Needs to be Bohr)
        #mol1.set_common_orig(com*1.88973)  # The gauge origin for dipole integral
        #mol1.set_rinv_orig(com*1.88973)  # The gauge origin for dipole integral
        #mol1.set_rinv_origin([0,0,0])
        mol1.set_common_orig([0, 0, 0])  # The gauge origin for dipole integral
        
        #h1 =(mol1.intor('int1e_kin') + mol1.intor('int1e_nuc') + np.einsum('x,xij->ij', E, mol1.intor('int1e_r', comp=3))) #mol.intor builds AO integrals
        
        h =(mol1.intor('cint1e_kin_sph') + mol1.intor('cint1e_nuc_sph') + np.einsum('x,xij->ij', E, mol1.intor('cint1e_r_sph', comp=3))) #mol.intor builds AO integrals
        
        #Compute nuclear dipole moment
        charges = mol1.atom_charges()
        r = mol1.atom_coords()
        nuc_dip = np.einsum('i,ix->x', charges, r) 
        if self.theory == 'RHF':
            mfx = scf.RHF(mol1)
            mfx.conv_tol = 1e-14
            
            #update hcore
            mfx.get_hcore = lambda *args: h
            
            #get density matrix
            mfx.scf(dm0[0])
            dm0[0] = mfx.make_rdm1()
            
            #print("coords after field:\n", mol1.atom_coords()/1.88973) 
            #vhfo = mfx.get_veff(mol1, dm0[0])
            #j, k = scf.hf.get_jk(mol1, dm0[0], hermi=0)
            #vhf_field = np.einsum('x,xij->ij', E, mfx.get_veff(mol1, dm0[0]))
            
            mol1.incore_anyway = True    #needed for post HF calculations to make sure custom H is used
            #e = mfx.kernel()
            e = mfx.kernel(dm=dm0[0]) 

            density = mfx.make_rdm1(mo_coeff=mfx.mo_coeff, mo_occ= mfx.mo_occ)
            print(density, density.shape)
            kinetic = mol1.intor('int1e_kin')
            potential = mol1.intor('int1e_nuc')
            overlap = mol1.intor('int1e_ovlp')
            print(kinetic, kinetic.shape)

            print(np.dot(density, kinetic))
            print("\potential:\n", potential)
            print(np.dot(density, potential))
            print("\overlap:\n", overlap)
            print(np.dot(density, overlap))
            vj, vk = mfx.get_jk(mol1, density)
            vhf1 = mfx.get_veff(mfx.mol, density)
            print("J:\n", vj)
            print('exchange K:\n', vk)
            fock = mfx.get_fock(h1e=h, s1e=overlap, vhf=vhf1, dm=density)
            print("fock\n:", fock)
            print(E)

            
            #print("total:",e)
            #print("nuclear", mol1.energy_nuc())
            g = mfx.nuc_grad_method().kernel()
            
            #gradient = mfx.Gradients()
            #g_elec = gradient.grad_elec()
            #g_nuc = gradient.grad_nuc()
            dipole = mfx.dip_moment(mol1, dm=dm0[0], unit='DEBYE')

            #compute nuclear energy
            print("nuclear repo energy:", mol1.energy_nuc())
            mfx.dump_scf_summary()
        
        if self.theory == 'MP2':
            mf = scf.RHF(mol1)
            mf.get_hcore = lambda *args: h
            mf.scf(dm0[0])
            dm0[0] = mf.make_rdm1()
            mol1.incore_anyway = True    #needed for post HF calculations to make sure custom H is used
            postmf = mp.MP2(mf)
            e = mf.kernel(dm=dm0[0]) + postmf.kernel()[0]
            g = postmf.nuc_grad_method().kernel()
            g_elec = None
            #gradient = mf.Gradients()
            #g_nuc = gradient.grad_nuc()
            nuc_dip = None
        
        return e, g, dipole, nuc_dip, g_nuc, g_elec
    
    def get_dipole(self, coords_new):
        """ This is used when building the APT's from numerical diff w.r.t 
        atomic coordinates.

        Pyscf only has dipole moments for SCF and DFT.
        """

        mol2 = gto.Mole()
        mol2.atom = coords_new
        mol2.basis = self.basis
        mol2.charge = self.charge
        if self.spin != 0:
            mol2.spin = self.spin-1
        else:
            mol2.spin = self.spin
        mol2.unit = 'Angstrom'
        mol2.build()

        if self.theory == 'RHF' or 'MP2':
            mfx = scf.RHF(mol2).run()
            dipole = mfx.dip_moment(mol2, unit='DEBYE')
        
        if self.theory == 'DFT':
            if mol2.spin != 0:
                mf = dft.UKS(mol2)
                mf.xc = self.xc
                mf.kernel()
                mf.analyze()
                dipole = mf.dip_moment(mol2, unit='DEBYE')

            else:
                mf = dft.RKS(mol2)
                mf.xc = self.xc
                mf.kernel()
                mf.analyze()
                dipole = mf.dip_moment(mol2, unit='DEBYE')
       
        if self.theory == 'UHF':
            mfx = scf.UHF(mol2).run()
            dipole = mfx.dip_moment(mol2, unit='DEBYE')
        
        if self.theory == 'ROHF': #Restricted Open shell HF calc
            mfx = scf.ROHF(mol2).run()
            dipole = mfx.dip_moment(mol2, unit='DEBYE')

        #Compute nuclear dipole moment
        charges = mol2.atom_charges()
        r = mol2.atom_coords()
        nucl_dip = np.einsum('i,ix->x', charges, r) 
        nuc_charge_center = nucl_dip/charges.sum()

        return dipole, nucl_dip
    




    


