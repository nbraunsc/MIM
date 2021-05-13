import string
import time
import numpy as np
#from .Pyscf import *
#from ase import Atoms
#from ase.calculators.vasp import Vasp
#from ase.vibrations import Infrared  
from numpy import linalg as LA 
from mendeleev import element

#import nicolefragment
#from nicolefragment import runpie, Molecule, fragmentation, Fragment, Pyscf

class Fragment():
    """
    Class to store a list of primitives corresponding to a molecular fragment
    
    Parameters
    ----------
    theory : str
        Level of theory for calculation
    basis : str
        Basis set name for calculations
    prims : list
        List of fragments from Fragmentation class with atom indexes in list
    attached : list
        List of attached pairs with the atom that is in the fragment and its corresponding atom pair that was cut
    coeff : int
        Coefficent of fragment.  This will either be 1 or -1.
    
    """
    
    def __init__(self, qc_class, molecule, prims, attached=[], coeff=1, step_size=0.001, local_coeff=1):
        self.prims = prims 
        self.molecule = molecule 
        self.coeff = coeff 
        self.attached = attached    #[(supporting, host), (supporting, host), ...]
        self.inputxyz = []
        self.apt = []
        self.aptgrad = np.array([])
        self.step = step_size
        self.energy = 0
        self.grad = 0
        self.hessian = 0
        self.hess = []
        self.notes = []     # [index of link atom, factor, supporting atom, host atom]
        self.jacobian_grad = [] #array for gradient link atom projections
        self.jacobian_hess = []  #ndarray shape of full system*3 x fragment(with LA)*3
        self.qc_class = qc_class
        self.step_size = step_size
        self.local_coeff = local_coeff
        self.M = [] #this is the mass matrix for massweighting shape: (3N, 3N)
        self.center = []
        self.gradlist = []
        self.origin_vec = []
        self.nuc_deriv = []

    def add_linkatoms(self, atom1, attached_atom, molecule):
        """ Adds H as a link atom
        
        This link atoms adds at a distance ratio between the supporting and host atom to each fragment where a previous atom was cut
        
        Parameters
        ----------
        atom1 : int
            This is the integer corresponding to the supporting atom (real atom)
        attached_atom : int
            This is the integer corresponiding to the host atom (ghost atom)
        molecule : <class> instance
            This is the molecule class instance

        Returns
        -------
        new_xyz : list
            This is the list of the new link atom with atom label and xyz coords
        factor : float
            The factor between the supporting and host atom. Used in building Jacobians for link atom projections.

        """

        atom1_element = molecule.atomtable[atom1][0]
        attached_atom_element = molecule.atomtable[attached_atom][0]
        cov_atom1 = molecule.covrad[atom1_element][0]
        cov_attached_atom = molecule.covrad[attached_atom_element][0]
        self.atom_xyz = np.array(molecule.atomtable[atom1][1:])
        attached_atom_xyz = np.array(molecule.atomtable[attached_atom][1:])
        vector = attached_atom_xyz - self.atom_xyz
        dist = np.linalg.norm(vector)
        h = 0.32
        factor = (h + cov_atom1)/(cov_atom1 + cov_attached_atom)
        new_xyz = list(factor*vector+self.atom_xyz)
        coord = []
        coord.append('H')
        coord.append(new_xyz)
        return coord, factor
    
    def build_xyz(self):
        """ Builds the xyz input with the atom labels, xyz coords, and link atoms as a string or list 
        
        Parameters
        ----------
        none
        
        Returns
        -------
        inputxyz : str
            String with atom label then corresonding xyz coordinates.  This input includes the link atoms.
        input_list : list of lists
            ie [[['H', [0, 0 ,0]], ['O', [x, y, z]], ... ] 
        self.notes: list of lists
            List of lists that is created with len = number of link atoms. Each sub list corresponds to one link atom.
            (i.e. [index of link atom, factor, supporting atom number, host atom number])
        
        """
        self.notes = []
        input_list = []
        coord_matrix = np.empty([len(self.prims)+len(self.attached), 3])
        for atom in self.prims:
            input_list.append([self.molecule.atomtable[atom][0]])
            input_list[-1].append(list(self.molecule.atomtable[atom][1:]))
            x = np.array(self.molecule.atomtable[atom][1:])
        for pair in range(0, len(self.attached)):
            la_input, factor = self.add_linkatoms(self.attached[pair][0], self.attached[pair][1], self.molecule)
            input_list.append(la_input)
            position = len(self.prims)+pair
            self.notes.append([position])
            self.notes[-1].append(factor)
            self.notes[-1].append(self.attached[pair][0])
            self.notes[-1].append(self.attached[pair][1])
        #self.input_list = input_list
        return input_list
    
    def build_jacobian_Grad(self):
        """Builds Jacobian matrix for gradient link atom projections
        
        Parameters
        ----------
        none
        
        Returns
        -------
        self.jacobian_grad : ndarray
            Array where entries are floats on the diagonal with the corresponding factor. 
            Array has size (# of atoms in full molecule + all link atoms, # of atoms in primiative)
        
        """
        self.jacobian_grad = 0
        array = np.zeros((self.molecule.natoms, len(self.prims)))
        linkarray = np.zeros((self.molecule.natoms, len(self.notes)))
        for i in range(0, len(self.prims)):
            array[self.prims[i]][i] = 1
        for j in range(0, len(self.notes)):
            factor = 1 - self.notes[j][1]
            linkarray[self.notes[j][2]][j] = factor
            linkarray[self.notes[j][3]][j] = self.notes[j][1]
        self.jacobian_grad = np.concatenate((array, linkarray), axis=1)
        jacob = self.jacobian_grad
        return jacob
    
    def build_jacobian_Hess(self):
        """ Builds Jacobian matrix for hessian link atom projections.

        Parameters
        ----------

        Returns
        -------
        self.jacobian_hess : ndarray (tensor)
            Array where the entries are matrices corresponding factor.
            
        """
        
        zero_list = []
        full_array = np.zeros((self.molecule.natoms, len(self.prims)+len(self.notes), 3, 3))

        for i in range(0, len(self.prims)):
            full_array[self.prims[i], i] = np.identity(3)
        for j in range(0, len(self.notes)):
            factor_s = 1-self.notes[j][1]
            factor_h = self.notes[j][1]
            x = np.zeros((3,3))
            np.fill_diagonal(x, factor_s)
            position = len(self.prims) + j
            full_array[self.notes[j][2]][position] = x
            np.fill_diagonal(x, factor_h)
            full_array[self.notes[j][3]][position] = x
        self.jacobian_hess = full_array
        return self.jacobian_hess

    def qc_backend(self):
        """ Runs the quantum chemistry backend.
        
        Returns
        -------
        self.energy : float
            This is the energy for the fragment*its coeff
        self.gradient : ndarray
            This is the gradient for the fragment*its coeff
        self.hessian : ndarray (4D tensor)
            This is the hessian for the fragement*its coeff
        """
        np.set_printoptions(suppress=True, precision=9, linewidth=200)
        self.energy = 0
        hess_py = 0
        self.grad = 0
        self.inputxyz = self.build_xyz()
        
        #sets origin of coords to center of mass
        self.center = self.com()   
        #finds inertia vector, R and T modes (only for 3 atom molecules currently)
        #self.inertia()
        
        energy, grad, hess_py = self.qc_class.energy_gradient(self.inputxyz)
        self.energy = self.local_coeff*self.coeff*energy
        jacob = self.build_jacobian_Grad()
        self.grad = self.local_coeff*self.coeff*jacob.dot(grad)
        self.M = self.mass_matrix()
        print("Done!")
        return self.energy, self.grad, hess_py  #, self.hessian#, self.apt

    def hess_apt(self, hess_py):
        
        #If not analytical hess, do numerical below
        if type(hess_py) is int:
            print("Numerical hessian needed, Theory=", self.qc_class.theory)
            hess_flat = np.zeros(((len(self.inputxyz))*3, (len(self.inputxyz))*3))
            i = -1
            for atom in range(0, len(self.inputxyz)):
                for xyz in range(0, 3):
                    i = i+1
                    self.inputxyz[atom][1][xyz] = self.inputxyz[atom][1][xyz]+self.step_size
                    grad1 = self.qc_class.energy_gradient(self.inputxyz)[1].flatten()
                    self.inputxyz[atom][1][xyz] = self.inputxyz[atom][1][xyz]-2*self.step_size
                    grad2 = self.qc_class.energy_gradient(self.inputxyz)[1].flatten()
                    self.inputxyz[atom][1][xyz] = self.inputxyz[atom][1][xyz]+self.step_size
                    vec = (grad1 - grad2)/(4*self.step_size)
                    hess_flat[i] = vec
                    hess_flat[:,i] = vec
        
        #Analytical hess from qc_backend gets reshaped and flatten to 3Nx3N matrix 
        else:
            hess_flat = hess_py
        
        #start building jacobian and reshaping
        self.jacobian_hess = self.build_jacobian_Hess() #shape: (Full, Sub, 3, 3)
        j_reshape = self.jacobian_hess.transpose(0,2,1,3)
        j_flat = j_reshape.reshape(self.molecule.natoms*3, len(self.inputxyz)*3, order='C')    #shape: (Full*3, Sub*3)
        j_flat_tran = j_flat.T  #shape: (Sub*3, Full*3)
        
        first = np.dot(j_flat, hess_flat)   # (Full*3, Sub*3) x (Sub*3, Sub*3) -> (Full*3, Sub*3)
        second = np.dot(first, j_flat_tran)     # (Full*3, Sub*3) x (Sub*3, Full*3) -> (Full*3, Full*3)
        self.hessian = second*self.coeff*self.local_coeff
        #start building the APT's
        self.apt = self.build_apt()    
        #self.aptgrad = self.apt_grad()     #one i am trying to get to work
        return self.hessian, self.apt
    
    def inertia(self):
        """ Finds principal axes and moments of inertia in amu*Bohr^2
        (I did this in a very non-optimized way!)
        """
        xx = 0
        yy = 0
        zz = 0
        xy = 0
        xz = 0
        yz = 0
        for i in range(0, len(self.inputxyz)):
            x = element(self.inputxyz[i][0])
            mass = x.atomic_weight
            xx += (self.inputxyz[i][1][1]**2 + self.inputxyz[i][1][2]**2)*mass
            yy += (self.inputxyz[i][1][0]**2 + self.inputxyz[i][1][2]**2)*mass
            zz += (self.inputxyz[i][1][0]**2 + self.inputxyz[i][1][1]**2)*mass
            xy += self.inputxyz[i][1][0]*self.inputxyz[i][1][1]*mass
            xz += self.inputxyz[i][1][0]*self.inputxyz[i][1][2]*mass
            yz += self.inputxyz[i][1][1]*self.inputxyz[i][1][2]*mass
        print("moment of interia for xx:", xx)
        print("moment of interia for yy:", yy)
        print("moment of interia for zz:", zz)
        print("moment of interia for xy:", xy)
        print("moment of interia for xz:", xz)
        print("moment of interia for yz:", yz)
        tensor = np.zeros((3,3))
        tensor[0][0] = xx
        tensor[0][1] = tensor[1][0] = xy
        tensor[1][1] = yy
        tensor[0][2] = tensor[2][0] = xz
        tensor[2][2] = zz
        tensor[1][2] = tensor[2][1] = yz
        print("Inertia tensor:\n", tensor)
        evalues, vec = LA.eig(tensor)  ###only for origin in pyscf calc
        #evalues, vec = LA.eigh(tensor)
        print(evalues)
        print(" Principal axes and moments of inertia in amu*Bohr^2:")
        print("Eigenvalues: \n", evalues*1.88973*1.88973)
        #vec[:, [2, 0]] = vec[:, [0, 2]]
        xyz = np.array(["X", "Y", "Z"])
        print(xyz[0], vec[0])
        print(xyz[1], vec[1])
        print(xyz[2], vec[2])
        
        #compute rotational constants
        conv = (6.626755E-34/(8*np.pi**2))/1.6605402E-27  #kg -> amu, cancel out all masses
        conv_final = (conv*1E20)/2.99792458E10 #B^2 -> A^2 -> m^2, cancel out all lengths, speed of light cm/s
        self.origin_vec = np.sqrt(conv/evalues)  #units of Bohr
        print("Pyscf origin vector:", self.origin_vec)
        rotate_const = conv_final/evalues
        print("Rotational constants (units: cm-1)\n", rotate_const)
        
        #generating internal coordinates to sep out R and T modes
        #self.int_coords(vec)


    def com(self):
        """ This is translating the origin of fragment to the center of mass.
        This will also update the coordinates for self.inputxyz to be in the center of mass basis.
        """

        first = 0
        second = 0
        for i in range(0, len(self.inputxyz)):
            x = element(self.inputxyz[i][0])
            mass = x.atomic_weight
            first += np.array(self.inputxyz[i][1])*mass
            second += mass
        self.center = (first/second)

        #update coordinates to COM in Bohr
        #for j in range(0, len(self.inputxyz)):
        #    self.inputxyz[j][1] = np.array(self.inputxyz[j][1]) - self.center
        return self.center
    
#    def int_coords(self, X):
#        """" Generate coordinates in teh rotating and translating frame.
#
#           This was trying to match Gaussian's way of computing the frequencies, taking out 
#            the rotational and translational modes, and IR intensities. 
#        """
#        R = np.zeros((len(self.inputxyz), 3))   #Coords in COM
#        M = np.zeros((len(self.inputxyz), 3))   #Mass 3x3 matrix with m^1/2
#        T = np.zeros((len(self.inputxyz), 3))   #Translation matrix 3x3
#        D = np.zeros((len(self.inputxyz)*3, 6))
#        D1 = np.array([1, 0, 0, 1, 0, 0, 1, 0, 0]).reshape((3,3))
#        D2 = np.array([0, 1, 0, 0, 1, 0, 0, 1, 0]).reshape((3,3))
#        D3 = np.array([0, 0, 1, 0, 0, 1, 0, 0, 1]).reshape((3,3))
#        
#        for i in range(0, R.shape[0]):
#            x = element(self.inputxyz[i][0])
#            mass = np.sqrt(x.atomic_weight)
#            M[i][i] = mass
#            D1[i] = D1[i]*mass
#            D2[i] = D2[i]*mass
#            D3[i] = D3[i]*mass
#            R[i] = np.array(self.inputxyz[i][1])
#        P = np.dot(R, X.T)
#        D1 = D1.flatten()
#        D2 = D2.flatten()
#        D3 = D3.flatten()
#        D4 = np.dot(np.outer(P[:,1], X[2]) - np.outer(P[:,2], X[1]), M).flatten()
#        print("D4:\n", np.dot(np.outer(P[:,1], X[2]) - np.outer(P[:,2], X[1]), M))
#        print("D5\n", np.dot(np.outer(P[:,2], X[0]) - np.outer(P[:,0], X[2]), M))
#        D5 = np.dot(np.outer(P[:,2], X[0]) - np.outer(P[:,0], X[2]), M).flatten()
#        print("D6\n", np.dot(np.outer(P[:,0], X[1]) - np.outer(P[:,1], X[0]), M))
#        D6 = np.dot(np.outer(P[:,0], X[1]) - np.outer(P[:,1], X[0]), M).flatten()
#        #print("D1\n", D1)
#        #print("D2\n", D2)
#        #print("D3\n", D3)
#        #print("D4\n", D4)
#        #print("D5\n", D5)
#        #print("D6\n", D6)
#        #print(D[:,0].shape)
#        #print(D1.shape)
#        D[:,0] = D1
#        D[:,1] = D2
#        D[:,2] = D3
#        D[:,3] = D4
#        D[:,4] = D5
#        D[:,5] = D6
#        #print(D, D.shape)
#        
#        #normalize D tensor
#        for j in range(0, D.shape[1]):
#            norm = 0
#            scalar = np.dot(D[:,j].T, D[:,j])
#            print(scalar)
#            if scalar < 1E-8:
#                continue
#            else:
#                norm = 1/np.sqrt(scalar)
#            D[:,j] = D[:,j]*norm
#
#        q, r = np.linalg.qr(D)
#        print(q, q.shape)
#        #exit()

    def apt_grad(self):
        """ Working on implementing this.
        Function to create the apts by applying an electric field in a certain direciton to 
        molecule then finite difference of gradient w.r.t the applied E field.

        Returns
        -------
        apt_grad : ndarray (3N, 3)
            The deriv of gradient w.r.t applied field after LA projections are done. 

        """
        extra_dip = self.qc_class.get_dipole(self.inputxyz)[0]
        
        #e_field = 1.889725E-4   #Got this number from Qchem
        e_field = 0.001
        E = [0, 0, 0]
        energy_vec = np.zeros((3))
        apt = np.zeros((3, ((len(self.prims)+len(self.notes))*3)))
        nucapt = np.zeros((3, ((len(self.prims)+len(self.notes))*3)))
        nuc3 = np.zeros((3)) 

        for i in range(0, 3):
            #no field 
            e1, g1, dip, n, g_nuc, g_elec = self.qc_class.apply_field(E, self.inputxyz, self.center, self.origin_vec, i)   #no field
            
            print("\n############ Field applied in the ", i, "direction ###############\n")
            #positive direction field
            E[i] = e_field
            e2, g2, dipole2, nuc2, g_nuc2, g_elec2 = self.qc_class.apply_field(E, self.inputxyz, self.center, self.origin_vec, i) #positive direction
            
            #negative direction field
            E[i] = -1*e_field
            e3, g3, dipole3, nuc, g_nuc3, g_elec3 = self.qc_class.apply_field(E, self.inputxyz, self.center, self.origin_vec, i)   #neg direction
            
            #setting field back to zero
            E[i] = 0
            
            print(g1)
            print(g2)
            print(g3)
            
            #central finite diff of gradient, a.u. -> Debye
            #print("positive grad:\n", g3, "\n Negative grad:\n", g2, "\n")
            gradient1 = ((g3-g2)/(2*e_field))/0.3934303
            print("$$$$$$$$$$$\n", gradient1)
            #add nuclear gradient to electronic
            gradient = g_nuc/0.393430 - gradient1    #for pyscf
            
            #checking finite diff of E w.r.t field (should be dipole moment)
            energy2 = (e2-e3)/(2*e_field)
            energy_vec[i] = energy2/0.393430 #a.u.(E_field) -> Debye, may need a neg sign
            
            #Subtracting elec dip from nuclear dip moment
            newvec = energy_vec  #for psi4
            #newvec = nuc3 - energy_vec #for pyscf 
            
            print("\nElectronic energy vec (Debye):", energy_vec, np.linalg.norm(energy_vec))
            print("\nNuclear dipole moment energy vec (Debye):", nuc3, np.linalg.norm(nuc3))
            print("\nDipole moment energy vec (Debye):", newvec, np.linalg.norm(newvec))
            print("\nDipole moment from no field (Debye):\n", extra_dip, np.linalg.norm(extra_dip))
            print("\ngradient no field", g1, "\n")
            print("\ngradient elec after finite diff:\n", gradient1)
            print("\ngradient nuc after finite diff:\n", g_nuc/0.393430)
            print("\ng_nuc - g_elec:\n", gradient)
            apt[i] = gradient1.flatten()
            #apt[i] = gradient.flatten()    #nuclear and electronic grad
       
       #mass weight APT        
        mass_apt = apt.T
        
        #Do link atom projection, multiply by local and principle inclusion/exculsion  coefficients
        reshape_mass_hess = self.jacobian_hess.transpose(0, 2, 1, 3)
        jac_apt = reshape_mass_hess.reshape(reshape_mass_hess.shape[0]*reshape_mass_hess.shape[1],reshape_mass_hess.shape[2]*reshape_mass_hess.shape[3])
        apt_grad = np.dot(self.M, self.local_coeff*self.coeff*np.dot(jac_apt, mass_apt))
        return apt_grad
            
    def build_apt(self):
        """
            Builds the atomic polar tensor with numerical derivative of dipole moment w.r.t atomic Cartesian
            coordinates. Function builds xyz input with link atoms in ndarray format, not string type or list like previous functions.

        Units of APT: Debye / (Angstrom np.sqrt(amu))

        Returns
        -------
        oldapt: ndarray (3N, 3)
            This is the mass weighted APT for current fragment after LA projections are done.

        """
        apt = []
        for atom in range(0, len(self.prims)+len(self.notes)):  #atom interation
            storing_vec = np.zeros((3,3))
            y = element(self.inputxyz[atom][0])
            value = 1/(np.sqrt(y.atomic_weight))
            for comp in range(0,3):   #xyz interation
                self.inputxyz[atom][1][comp] = self.inputxyz[atom][1][comp]+self.step_size
                dip1, nuc1 = self.qc_class.get_dipole(self.inputxyz)
                self.inputxyz[atom][1][comp] = self.inputxyz[atom][1][comp]-2*self.step_size
                dip2, nuc2 = self.qc_class.get_dipole(self.inputxyz)
                vec = (dip1 - dip2)/(2*self.step_size)
                storing_vec[comp] = vec
                self.inputxyz[atom][1][comp] = self.inputxyz[atom][1][comp]+self.step_size
            apt.append(storing_vec)
        px = np.vstack(apt)
        
        reshape_mass_hess = self.jacobian_hess.transpose(0, 2, 1, 3)
        jac_apt = reshape_mass_hess.reshape(reshape_mass_hess.shape[0]*reshape_mass_hess.shape[1],reshape_mass_hess.shape[2]*reshape_mass_hess.shape[3])
        
        oldapt = np.dot(self.M, self.local_coeff*self.coeff*np.dot(jac_apt, px)) #mass weight here and LA projection
        return oldapt

    def mass_matrix(self):
        M = np.zeros((self.molecule.natoms*3, self.molecule.natoms*3))
        counter = np.array([0, 1, 2])
        for i in range(0, self.molecule.natoms):
            x = element(self.molecule.atomtable[i][0])
            value = 1/(np.sqrt(x.atomic_weight))
            for j in counter:
                M[j][j] = value
            counter = counter + 3
        self.M = M
        return self.M

    def mw_hessian(self, full_hessian):
        """
        Will compute the mass-weighted hessian, frequencies, and 
        normal modes for the full system.
        
        Parameters
        ----------
        full_hessian : ndarray
            This is the full hessian for the full molecule.

        Returns
        -------
        freq : ndarray
            1D np array holding the frequencies
        modes : ndarray
            2D ndarray holding normal modes in the columns
        """
        np.set_printoptions(suppress=True)
        first = np.dot(full_hessian, self.M) #shape (3N,3N) x (3N, 3N)
        second = np.dot(self.M, first)   #shape (3N,3N) x (3N, 3N)
        
        e_values, modes = LA.eigh(second)
        print("\nEvalues of hessian [H/Bohr^2]):\n", e_values)

        #unit conversion of freq from H/B**2 amu -> 1/s**2
        #factor = (4.3597482*10**-18)/(1.6603145*10**-27)/(1.0*10**-20) # Hartreee->J, amu->kg, Angstrom->m  
        factor = (1.8897259886**2)*(4.3597482*10**-18)/(1.6603145*10**-27)/(1.0*10**-10)**2 #Bohr->Angstrom, Hartreee->J, amu->kg, Angstrom->m
        freq = (np.sqrt(e_values*factor))/(2*np.pi*2.9979*10**10) #1/s^2 -> cm-1
        return freq, modes, self.M, e_values

