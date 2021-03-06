from collections import Counter
import cProfile
import time
import pandas as pd
import os
import numpy as np
from numpy import linalg as LA 
import sys
from sys import argv
import xml.etree.ElementTree as ET
from itertools import cycle
from mendeleev import element
from atom_count import *
import time

#import mim
#from mim import runpie, Molecule, fragmentation, Fragment, Pyscf, oldpie, newpie, newnewpie
from mim import Molecule, fragmentation, Fragment, Pyscf, newnewpie

class Fragmentation():
    """
    Fragmentation is a child class of Molecule().
    
    Used to build the molecular fragments and their derivatives. An instance of this class 
    represents a mapping between a large molecule and a list of independent fragments with 
    the appropriate coefficients needed to reassemble expectation values of observables. 
    """
    
    def __init__(self, molecule, folder): 
        self.fragment = []
        self.molecule = molecule
        self.folder = folder
        self.unique_frag = []
        self.derivs = []
        self.coefflist = []
        self.atomlist = []
        self.frags = []     #list of Fragment class instances
        self.etot = 0
        self.gradient = []
        self.hessian = []
        self.fullgrad = {}  #dictonary for full molecule gradient
        self.fullhess = {}
        self.moleculexyz = []   #full molecule xyz's
        self.etot_opt = 0
        self.grad_opt = []
        self.name = hex(id(self))

    def build_frags(self, frag_type=None, value=None):    #deg is the degree of monomers wanted
        """ Performs the initalize fragmmentation 
        
        This will fragment the full molecule at the level of deg specified.
        Also runs self.compress_frags() to delete repeating fragments.
        
        Parameters
        ----------
        frag_type : str
            This is the type of fragmentation wanted either distance or graphical.

        value : int
            The degree of fragmentation wanted. Larger deg results in larger fragments. For distance type it will
            be the largest radius cuttoff. For covalent_net it will be the number of bonds away from a center.
        
        Returns
        -------
        self.fragment : list
            List of fragments containing index of primiatives that are within the value
        
        """

        self.fragment = []
        
        if frag_type == 'graphical':
            self.molecule.build_primchart()

            #do atom centric fragmentation
            if (value).is_integer():
                self.molecule.build_molmatrix()
                for prim in range(0, self.molecule.molchart.shape[0]):
                    x = list(np.where(self.molecule.molchart[prim] <= value)[0])
                    self.fragment.append(x)

            #do bond centric fragmention
            else:
                frag_list = []
                for prim in range(0, self.molecule.primchart.shape[0]):
                    x = list(np.where(self.molecule.primchart[prim] == 1)[0])
                    for i in x:
                        y = [prim, i]
                        frag_list.append(y)
                
                if value == 0.5:
                    self.fragment = frag_list
                
                if value > 0.5:
                    self.molecule.build_molmatrix()
                    bonds = value-0.5
                    extend_fraglist = []
                    #value = 0.5 -> only one bond, 1.5 -> one bond away from each prim in original frags
                    for frag in frag_list:
                        new_frag = frag
                        primlist = list(frag)
                        for prim in primlist:
                            x = list(np.where(self.molecule.molchart[prim] <= bonds)[0])
                            new_frag = new_frag + x
                        extend_fraglist.append(set(new_frag))
                    self.fragment = [list(i) for i in extend_fraglist]


        if frag_type == 'distance':
            self.molecule.prim_dist = self.molecule.build_prim_dist()
            for a in range(0, len(self.molecule.prims)):
                y = list(np.where(self.molecule.prim_dist[a] <= value)[0])
                self.fragment.append(y)
        
        new_frag_list = []
        for frag in self.fragment:
            #get all connected prims
            connected = []
            for prim in frag:
                connected.extend(list(np.where(self.molecule.primchart[prim] > 0)[0]))
            
            #find all duplicates in connected list
            duplicates = set([x for x in connected if connected.count(x) > 1])
            
            #if there are duplicates, add in primitives that would have LA added in same position
            if len(duplicates) == 0:
                new_frag_list.append(frag)
                continue
            else:
                new_frag = self.fix_superposition(frag, duplicates)
                new_frag_list.append(new_frag)
        
        # Now get list of unique frags, running compress_frags function below
        #self.compress_frags(self.fragment)
        self.compress_frags(new_frag_list)
        self.write_xyzfiles(self.unique_frag)
        
        #exit()

        #code for me to check stuff
        #fragment = []
        #self.molecule.build_molmatrix()
        #for prim in range(0, self.molecule.molchart.shape[0]):
        #    x = list(np.where(self.molecule.molchart[prim] <= 3)[0])
        #    fragment.append(x)
        
        #for i in range(0, len(self.molecule.prims)):
        #    print("distance:\n", self.fragment[i])
        #    for prim in self.fragment[0]:
        #        atoms = self.molecule.prims[prim].atoms
        #        for atom in atoms:
        #            print(str(self.molecule.atomtable[atom]).replace("[", "").replace("]", "").replace(",", "").replace("'", ""))
        #    
        #    print("bond order:\n", fragment[i])
        #    for prim in fragment[0]:
        #        atoms = self.molecule.prims[prim].atoms
        #        for atom in atoms:
        #            print(str(self.molecule.atomtable[atom]).replace("[", "").replace("]", "").replace(",", "").replace("'", ""))
        #    exit()
            
            #if self.fragment[i] != fragment[i]:
            #    print("Fragment #", i, "doesnt match")
            #    #raise ValueError
            #else:
            #    print("matching!")
            #    continue
        #print("All match")

    def write_xyzfiles(self, fraglist):
        """ Write the xyz files for all fragments before derivatives are formed
        """
        path = os.getcwd()
        loc = path + "/" + self.folder + "traj.xyz"
        filename = open(loc, "w")
        #filename.write(str(self.molecule.natoms)+"\n")
        for frag in fraglist:
            filename.write(str(self.molecule.natoms)+"\n")
            filename.write("Fragment" + str(frag))
            total_str = "" + "\n"
            all_atoms = []
            list_2 = list(range(0, self.molecule.natoms))
            for prim in frag:
                atoms = self.molecule.prims[prim].atoms
                all_atoms.extend(atoms)
                #build xyz for fragment
                for atom in atoms:
                    coords = self.molecule.atomtable[atom]
                    string = str(coords).replace("'", "").replace("[", "").replace("]", "").replace(",", "") + "\n"
                    total_str += string
            
            #build xyz for dummie atoms (label = Si)
            main_list = list(set(list_2) - set(all_atoms))
            for atomi in main_list:
                dummie_coords = self.molecule.atomtable[atomi].copy()
                dummie_coords[0] = 'X'
                stringi = str(dummie_coords).replace("'", "").replace("[", "").replace("]", "").replace(",", "") + "\n"
                total_str += stringi

            filename.write(total_str)
        filename.close()


    def fix_superposition(self, frag, duplicates):
        """ Removes occurances where link atoms would be added in the same location

        Parameters
        ----------
        frag : list
            List of primitives in current fragment

        Returns
        -------
        new_frag : list
            New fragment with a list of primitives (longer list)

        """
        new_frag = []
        new_frag = frag
        for i in list(duplicates):
            if i not in frag:
                new_frag = new_frag + [i]
            else:
                continue 
        
        if len(new_frag) == 0:
            return frag

        new_connected = []
        for prim in new_frag:
            new_connected.extend(list(np.where(self.molecule.primchart[prim] > 0)[0]))
        new_duplicates = set([x for x in new_connected if new_connected.count(x) > 1])
        
        #check if new fragment also has redundancy with link atoms
        if duplicates != new_duplicates:
            new_frag = self.fix_superposition(new_frag, new_duplicates)
        
        return new_frag

    def compress_frags(self, frag_list): #takes full list of frags, compresses to unique frags
        """ Takes the full list of fragments and compresses them to only the unique fragments.
        
        This functions compresses self.fragment into self.uniquefrag. Deletes fragments that
        are completely represented by another fragment.

        Parameters
        ----------
        none
        
        Returns
        -------
        self.unique_frag : list
            List of sets containing only the unqie fragments 
        
        """
        frag_list = [set(i) for i in frag_list]
        sign = 1
        uniquefrags = []
        for i in frag_list:
            add = True
            for j in frag_list:
                if i.issubset(j) and i != j:
                    add = False
           
            if add == True:
                if i not in uniquefrags:
                    uniquefrags.append(i)   
        self.unique_frag = uniquefrags
        print(self.unique_frag)
        print(frag_list)
        print(len(frag_list), len(self.unique_frag))
        print("Diff:", len(frag_list)-len(self.unique_frag))

        for i in self.unique_frag[0]:
            atoms = self.molecule.prims[i].atoms
            for atom in atoms:
                print(str(self.molecule.atomtable[atom]).replace("]", "").replace("[", "").replace(",", "").replace("'", ""))
        
        #print(self.unique_frag)
        #print("len of unique:", len(self.unique_frag))
        #print("len of orig:", len(self.fragment))
        
        #for frag in self.unique_frag:
        #    print(frag, len(frag))
        #total = []
        #for frag in self.unique_frag:
        #    atoms = []
        #    for prim in frag:
        #        atoms.extend(self.molecule.prims[prim].atoms)
        #    for atom in atoms:
        #        print(str(self.molecule.atomtable[atom]).replace("]", "").replace("[", "").replace(",", "").replace("'", ""))
        #    atoms = self.molecule.prims[i].atoms
        #    total.extend(atoms)
        #print(total)
        #print(len(total))

    def find_attached(self):    #finding all the attached atoms that were cut during fragmentation
        """ Finds the atoms that were cut during the fragmentation.
        
        Parameters
        ----------
        none
        
        Returns
        -------
        self.attached : list of lists
            List of pairs of atoms (a, b) len of # of fragments, where a is in fragment and atom b was the 
            atom that was cut.
        
        """
        
        self.attached = []
        for frag in range(0, len(self.atomlist)):
            fragi = []
            for atom in self.atomlist[frag]:
                row = self.molecule.A[atom]
                for i in range(0, len(row)):
                    if row[i] == 1.0:
                        attached = i
                        if attached not in self.atomlist[frag]:
                            fragi.append([atom, i])
                        else:
                            continue
            self.attached.append(fragi) 


#    def remove_repeatingfrags(self, oldcoeff, derivs):
#        """ Removes the repeating derivatives. Function gets called in self.do_fragmentation().
#        
#        The repeating derivatives are the exact same derivatives that would be added then subtracted
#        again during the principle of inclusion-exculsion process of MIM.  This also updates the self.derivs
#        and self.coefflist that are used with Fragment() class instances.
#        
#        Parameters
#        ----------
#        oldcoeff : list
#            List of coefficients that is the output of runpie() function that is called in self.do_fragmentation().
#        
#        Returns
#        -------
#        self.coefflist : list
#            List of coefficents where index of coeff correlates to that fragment
#            
#        """
#        start = time.time()
#        derv = []
#        coeff =[]
#        
#        #x = Counter(map(tuple, derivs))
#        #keylist = list(x.keys())
#        #print(x)
#        #print("keylist:", len(keylist))
#        #print("keylist:", keylist, len(keylist))
#
#
#        #testing sorted thing
#        sortlist = []
#        for frag in derivs:
#            #sorted makes into a list
#            newfrag = sorted(frag)
#            sortlist.append(newfrag)
#        
#        y = Counter(map(tuple, sortlist))
#        #keylist is a list of tuples
#        keylist = list(y.keys())
#        
#        for i in keylist:
#            count = y[i]
#            #count = x[i]
#            if count == 1:
#                derv.append(set(i))
#                position = sortlist.index(list(i))
#                #position = derivs.index(set(i))
#                coeff.append(oldcoeff[position])
#            else:
#                #indices = [index for index, element in enumerate(derivs) if element == set(i)]
#                indices = [index for index, element in enumerate(sortlist) if element == list(i)]
#                sum_value = 0
#                for num in indices:
#                    sum_value += oldcoeff[num]
#                
#                #print("Fragment", i, "has sum value of :", sum_value)
#
#                if sum_value > 0:
#                    for k in range(0, abs(sum_value)):
#                        derv.append(set(i))
#                        coeff.append(1)
#                    #print("now deriv:", derv)
#                    #print("now coeff:", coeff)
#
#                if sum_value < 0:
#                    for j in range(0, abs(sum_value)):
#                        derv.append(set(i))
#                        coeff.append(-1)
#                    #print("now deriv:", derv)
#                    #print("now coeff:", coeff)
#
#        #print(derv)
#        end = time.time()
#        print("\nTotal time of remove repeating:", end-start)
#        self.coefflist = coeff
#        print("\nOriginal length of dervs from pie:", len(derivs))
#        print("Final length after remove repeating:", len(derv))
#


        
        ##unique = [list(tupl) for tupl in {tuple(item) for item in derivs }]
        #unique = [list(x) for x in set(tuple(x) for x in derivs)]
        #
        ##print("Unique fragments:\n", unique)
        #print("There are ", len(unique), "unique frags out of ", len(derivs))
        #exit()
        #derv = []
        #coeff =[]

        #for x in range(0, len(unique)):
        #    index = []
        #    value = 1
        #    count = derivs.count(set(unique[x]))
        #    
        #    #if derv only appears once add it to final list
        #    if count == 1:
        #        derv.append(unique[x])
        #        index = derivs.index(set(unique[x]))
        #        coeff.append(oldcoeff[index])

        #    #if derv appears more than once, add/subtract coeff to determine if needed
        #    if count > 1:
        #        indices = [index for index, element in enumerate(derivs) if element == set(unique[x])]
        #        temp_coef = []
        #        for num in indices: 
        #            temp_coef.append(oldcoeff[num])
        #        sum_coeff = float(sum(temp_coef))
        #        
        #        #if overall neg, add derv that amount of times with -1 as coeff
        #        if sum_coeff < 0:
        #            for i in range(0, int(abs(sum_coeff))):
        #                derv.append(unique[x])
        #                coeff.append(-1)

        #        #if overall positive, add derv that amount of times with +1 as coeff
        #        if sum_coeff > 0:
        #            for i in range(0, int(sum_coeff)):
        #                derv.append(unique[x])
        #                coeff.append(1)
        #        
        #        #if sum_coeff == 0:
        #            #print("indices >0, but sum of coeff ==0, dont add")
        #            #print(unique[x])
        #    
                
        #self.coefflist = coeff

        #derivatives = [list(tupl) for tupl in {tuple(item) for item in derivs }]

        #for i in derivatives:
        #    key = str(sorted(set(i)))
        #    self.coefflist.append(derv_dict[key])
        #print(derivatives)
        #print(self.coefflist)
        #return derv

    def initalize_Frag_objects(self, theory=None, basis=None, qc_backend=None, spin=None, tol=None, active_space=None, nelec=None, nelec_alpha=None, nelec_beta=None, max_memory=None, xc=None, charge=0, step_size=0.001, local_coeff=1):
        """ Initalizes the Fragment() objects
        
        Fragment is another child class of Fragmentaion() where link atoms get added, 
        fragment energy, gradient and hessian is calcuated.
        
        Parameters
        ----------
        theory : str
            Level of theory each fragment should be run at
        basis : str
            Basis set wanted for calculation
        
        Returns
        -------
        none
        
        """
        
        self.frags = []

        for fi in range(0, len(self.atomlist)):
            coeffi = self.coefflist[fi]
            attachedlist = self.attached[fi]
            frag_charge = 0
            frag_spin = 0
            for prim in self.derivs[fi]:
                frag_charge += self.molecule.prims[prim].charge
                frag_spin += self.molecule.prims[prim].spin
            
            qc_fi = qc_backend(theory=theory, basis=basis, tol=tol, active_space=active_space, nelec_alpha=nelec_alpha, nelec_beta=nelec_beta, max_memory=max_memory, xc=xc, charge=frag_charge, spin=frag_spin)
        
            self.frags.append(Fragment.Fragment(qc_fi, self.molecule, self.atomlist[fi], attachedlist, coeff=coeffi, step_size=step_size, local_coeff=local_coeff))
        print("done with initlaze frags")
        #for frag in self.frags:
        #    inputlist = frag.build_xyz()
        #    print("\nNew frag:\n")
        #    for atom in inputlist:
        #        print(str(atom).replace("]", "").replace("[", "").replace(",", "").replace("'", ""))
        #exit()

    def write_xyz(self, name):
        """ Writes an xyz file with the atom labels, number of atoms, and respective Cartesian coords for geom_opt().

        Parameters
        ----------
        name : str
            Name of the molecule, must be same name as input.cml file
        
        Returns
        -------
        none
        
        """
        molecule = np.array(self.molecule.atomtable)
        atomlabels = []
        for j in range(0, len(molecule)):
            atomlabels.append(molecule[j][0])
        coords = molecule[:, [1,2,3]]
        self.moleculexyz = []
        for i in coords:
            #x = [i[1], i[2], i[3]]
            y = np.array(i)
            z = y.astype(float)
            self.moleculexyz.append(z)
        self.moleculexyz = np.array(self.moleculexyz)   #formatting full molecule coords
        
        f = open(name, "w+")
        title = ""
        f.write("%d\n%s\n" % (self.moleculexyz.size / 3, title))
        for x, atomtype in zip(self.moleculexyz.reshape(-1, 3), cycle(atomlabels)): 
            f.write("%s %.18g %.18g %.18g\n" % (atomtype, x[0], x[1], x[2]))
        f.close()
    
    def attached_frags(self, fraglist):
        prim_dict = {}
        for prim in range(0, len(self.molecule.prims)):
            tmp = []
            for frag in range(0, len(fraglist)):
                if prim in fraglist[frag]:
                    tmp.append(frag)
            prim_dict[prim]=tmp
        return prim_dict


    def do_fragmentation(self, fragtype=None, value=None):
        """ Main executeable for Fragmentation () class

        This function fragments the molecule, runs principle of inclusion-exculsion,
        removes repeating fragments, creates attached pairs list, initalizes instances of Fragment() class.
        
        Parameters
        ----------
        frag_type : str
            Either 'graphical' to do covalent network or 'radius' to do spacial fragmentation
        value : int
            Degree or radius of fragmentation wanted, gets used in self.build_frags()
        
        Returns
        -------
        none
        
        """
        ### Builds inital fragments and pulls out frags that are fully represented by another fragment
        self.build_frags(frag_type=fragtype, value=value)
        #print(self.unique_frag, len(self.unique_frag))

        ### Builds dict to hold which primitives are in which fragments
        att_dict = self.attached_frags(self.unique_frag)

        ### Start of Priniciple inclusion/exculsion
        #derivs, oldcoeff = runpie.start_pie(self.unique_frag, att_dict)
        #derivs, oldcoeff = oldpie.runpie(self.unique_frag)
        #derivs, oldcoeff = newpie.start_pie(self.unique_frag, att_dict)

        ## Function profile
        #string = 'mim.newnewpie.start_pie(' + str(self.unique_frag) + ', ' + str(att_dict) + ')'
        #cProfile.run(string, 'restats')
        #exit()


        start = time.time()
        derv_dict = newnewpie.start_pie(self.unique_frag, att_dict)
        end = time.time()
        
        #pulling out all non zero coeff fragments
        #new_dervs = []
        #new_coeff = []
        for frag in list(derv_dict.keys()):
            value = derv_dict[frag]
            if value == 0:
                continue
            else:
                self.derivs.append(list(frag))
                self.coefflist.append(value)

        #self.coefflist = new_coeff 
        #self.derivs = new_dervs

        ### turning derivs into a list of atoms instead of a list of primitives
        self.atomlist = []
        for j in self.derivs:
            temp = []
            for prim in j:
                temp.extend(list(self.molecule.prims[prim].atoms))
            self.atomlist.append(temp)


        #for derv in self.atomlist:
        #    print("\nNew frag:", derv, "\n")
        #    for atom in derv:
        #        print(str(self.molecule.atomtable[atom]).replace("]", "").replace("[", "").replace(",", "").replace("'", ""))
        
        ### Counting # atoms in largest fragment
        count = []
        for i in self.atomlist:
            count.append(len(i))
        large = np.argmax(count)
        print("\nLargest deriv is frag #", large)
        
        print("With", len(self.atomlist[large]), "atoms before adding link atoms\n")
        
        ### testing to make sure all atoms are only counted once
        vec = test_atoms(self.atomlist, self.coefflist, self.molecule.natoms)
        print("Vec should be all 1's:\n", vec)
        print("\nTime of pie:", end-start)
        print("Total # of dervs:", len(self.atomlist))
        summ = sum(vec)
        if summ != len(vec):
            raise ValueError("Not all atoms are counted only once!")

        ### Removes derivs that would otherwise get added then subtracted
        #self.derivs = self.remove_repeatingfrags(oldcoeff, derivs)
        #print(self.derivs)
        
        ### turning derivs into a list of atoms instead of a list of primitives
        #self.atomlist = []
        #for j in self.derivs:
        #    temp = []
        #    for prim in j:
        #        temp.extend(list(self.molecule.prims[prim].atoms))
        #    self.atomlist.append(temp)
        
        ### Counting # atoms in largest fragment
        #count = []
        #for i in self.atomlist:
        #    count.append(len(i))
        #large = np.argmax(count)
        #print("\nLargest deriv is frag #", large)
        
        #print("With", len(self.atomlist[large]), "atoms\n")

        ### testing to make sure all atoms are only counted once
        #vec = test_atoms(self.atomlist, self.coefflist, self.molecule.natoms)
        #print("Vec should be all 1's:\n", vec)
        #print("\nTime of pie:", end-start)
        #summ = sum(vec)
        #if summ != len(vec):
        #    raise ValueError("Not all atoms are counted only once!")
        #exit()
        
        #for value in vec:
        #    if value != 1:
        #        raise ValueError("Not all atoms are counted only once!")
        #        break
        
        #exit()

        ### Finding the attached atoms to help with adding link atoms in Fragment() class
        self.find_attached()
    
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
        mass_xyz = np.zeros((self.molecule.natoms, 3))
        mass_array = np.zeros((self.molecule.natoms))
        for row in range(0, self.molecule.natoms):
            #adding 1/np.sqrt(amu) units to xyz coords
            x = element(self.molecule.atomtable[row][0])
            value = np.sqrt(x.atomic_weight)
            mass_array[row] = value
            mass_xyz[row] = np.array(self.molecule.atomtable[row][1:])*(1/value)   #mass weighted coordinates
            for column in range(0, len(self.molecule.atomtable)):    
                y = element(self.molecule.atomtable[column][0])
                z = x.atomic_weight*y.atomic_weight
                value = np.sqrt(z)**-1
                term = full_hessian[row][column]*value
                full_hessian[row][column] = term
                
        reshape_mass_hess = full_hessian.transpose(0, 2, 1, 3)
        x = reshape_mass_hess.reshape(reshape_mass_hess.shape[0]*reshape_mass_hess.shape[1],reshape_mass_hess.shape[2]*reshape_mass_hess.shape[3])
        e_values, modes = LA.eigh(x)

        #unit conversion of freq from H/B**2 amu -> 1/s**2
        #factor = (4.3597482*10**-18)/(1.6603145*10**-27)/(1.0*10**-20)  #Angstrom to m
        factor = 1.8897259886**2*(4.3597482*10**-18)/(1.6603145*10**-27)/(1.0*10**-20) #Bohr to Angstrom
        freq = (np.sqrt(e_values*factor))/(2*np.pi*2.9979*10**10) #1/s^2 -> cm-1
        return freq, modes
    
    #def qc_params(self, qc_backend=None, theory=None, basis=None, spin=None, tol=None, active_space=None, nelec=None, nelec_alpha=None, nelec_beta=None, max_memory=None, xc=None, charge=0):
    ##def qc_params(self, frag_index=[], qc_backend=None, theory=None, basis=None, spin=0, tol=0, active_space=0, nelec=0, nelec_alpha=0, nelec_beta=0, max_memory=0):
    #    """ Funciton that is optional
    #    
    #    Use if a certain fragment or list of fragments need to be run with a different qc backend, or theory, or needs specific params.

    #    Parameters
    #    ----------
    #    frag_index : list
    #        The location of the specific fragment within the fragment list
    #    qc_backend : str
    #        Quantum chemistry backend wanted for this fragment
    #    theory : str
    #        Theory wanted
    #    basis : str
    #        Basis set wanted
    #    spin : int
    #        Spin state for fragment
    #    tol : int
    #        Tolerance wanted
    #    active_space : int
    #        Size of active space for CASSCF etc.
    #    nelec_alpha : int
    #        Number of alpha electrons
    #    nelec_beta : int
    #        Number of beta electrons
    #    max_memory : int
    #        Max memory for this calculation
    #    """
    #    #for fi in frag_index:
    #    qc_fi = qc_backend(theory=theory, basis=basis, spin=spin, tol=tol, active_space=active_space, nelec_alpha=nelec_alpha, nelec_beta=nelec_beta, max_memory=max_memory, xc=xc, charge=charge)
    #    return qc_fi
            
    
    #def energy_gradient(self):
    ##def energy_gradient(self, newcoords):
    #    """ Function returns total energy and gradient of global molecule.
    #    
    #    This function holds virtual functions for different chemical software.  A software
    #    must be implemented inorder to run this function.  
    #    
    #    Parameters
    #    ----------
    #    newcoords : npdarray
    #        Contains xyz coords for the full molecule. These get updated after each geometry optimization cycle.
    #    
    #    Returns
    #    -------
    #    self.etot : float
    #        Energy of the full molecule
    #    self.gradient : ndarray
    #        Gradient of the full molecule
    #    
    #    """

    #    #for atom in range(0, len(self.molecule.atomtable)): #makes newcoords = self.molecule.atomtable
    #    #    x = list(newcoords[atom])
    #    #    self.molecule.atomtable[atom][1:] = x


    #    self.gradient = np.zeros((self.molecule.natoms,3)) #setting them back to zero
    #    self.etot = 0
    #    self.hessian = np.zeros((self.molecule.natoms, self.molecule.natoms, 3, 3)) 
    #    apt=0
    #    for i in self.frags:
    #        #i.molecule.atomtable = self.molecule.atomtable  #setting newcoords
    #        e, grad, i_hess, i_apt = i.qc_backend()
    #        self.etot += i.energy
    #        self.gradient += i.grad
    #        self.hessian += i.hessian
    #        apt += i_apt
    #    return self.etot, self.gradient, self.hessian, apt
        
    #def do_geomopt(self, name, theory, basis):
    #    """ Completes the geometry optimization using pyberny from pyscf.
    #    
    #    Parameters
    #    ----------
    #    name : str
    #        Name of Molecule() object to help with file path
    #    theory : str
    #        Level of theory for calculation
    #    basis : str
    #        Basis set name for calculation
    #    
    #    Returns
    #    -------
    #    self.etot_opt : float
    #        Optimized energy for full molecule
    #    self.grad_opt : ndarray
    #        Optimized gradient for full molecule
    #    
    #    """
    #    self.write_xyz(name)
    #    os.path.abspath(os.curdir)
    #    os.chdir('../inputs/')
    #    optimizer = Berny(geomlib.readfile(os.path.abspath(os.curdir) + '/' + name + '.xyz'), debug=True)
    #    for geom in optimizer:
    #        solver = self.energy_gradient(geom.coords)
    #        optimizer.send(solver)
    #        self.etot_opt = solver[0]
    #        #self.grad_opt = solver[1]
    #    relaxed = geom
    #    print("\n", "##########################", '\n', "#       Converged!       #", '\n', "##########################") 
    #    print('\n', "Energy = ", self.etot_opt)
    #    self.molecule.optxyz = relaxed.coords
    #    os.chdir('../')
    #    return self.etot_opt#, self.grad_opt
       

    #def global_apt(self):
    #    global_apt = np.zeros((self.molecule.natoms*3, 3))
    #    for frag in self.frags:
    #        #frag.build_apt()
    #        global_apt += frag.apt
    #    return global_apt


if __name__ == "__main__":
    #hf = Molecule()
    #hf.initalize_molecule('hf')
    #frag = Fragmentation(hf)
    #frag.do_fragmentation(frag_type='graphical', value=2)
    #frag.initalize_Frag_objects(theory='RHF', basis='sto-3g', qc_backend=Pyscf)
    #frag.energy_gradient(frag.moleculexyz)
    np.set_printoptions(suppress=True, precision=6)
    
    hf = Molecule.Molecule('/Users/nicole/Documents/research/Fragments/inputs/example_molecules/hf.cml')
    hf.initalize_molecule()
    frag = fragmentation.Fragmentation(hf)
    frag.do_fragmentation(fragtype='distance', value=22)
    frag.initalize_Frag_objects(theory='RHF', basis='sto-3g', qc_backend=Pyscf.Pyscf, step_size=0.001, local_coeff=1)
    frag.energy_gradient()
   
    #import ray
    #start_time = time.time()
    #ray.init()

    #frags_id = ray.put(frag)    #future for Fragmentation instance, putting in object store

    #@ray.remote
    #def get_frag_stuff(f,_frags):
    #    f_current = _frags.frags[f]
    #    return f_current.qc_backend()
    #result_ids = [get_frag_stuff.remote(fi, frags_id) for fi in range(len(frag.frags)) ]
    #out = ray.get(result_ids)
    #etot = 0
    #gtot = 0
    #htot = 0
    #for o in out:
    #    etot += o[0]
    #    gtot += o[1]
    #    htot += o[2]
    #total_time = time.time() - start_time 
    #ray.shutdown()
    #e, v = frag.mw_hessian(htot)
    #print("e", e)
    #print("modes", v)

    #start_fragtime = time.time()
    #frag.write_xyz('carbonylavo')
    #frag.energy_gradient(frag.moleculexyz)
    #end_time = time.time() - start_fragtime
    #print("Final converged energy = ", etot)
    #print("Final gradient = ", '\n', gtot)
    #print("Final hessian = ", '\n', htot)
    #print("Hessian shape = ", htot.shape)
    #print(htot[0][1])
    #print(htot[1][0])


    ### run qc_params only if you want some fragments with different params
    #frag.qc_params(frag_index=[], qc_backend, theory, basis, spin=0, tol=0, active_space=0, nelec_alpha=0 nelec_beta=0, max_memory=0)
    
