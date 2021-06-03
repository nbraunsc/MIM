import os
import numpy as np
import sys
from sys import argv
import xml.etree.ElementTree as ET
#import pandas as pd
#from pandas.io.json import json_normalize
#import json

#import mim
from mim import cov_rad, Primitive

class Molecule():
    """
    Molecule class is the parent class for the MIM code.
    
    Responsible for obtaining input file of molecule, changing into
    appropriate file format
    
    Holds molecule information that will not change througout fragmentation
    code.
    """
    
    def __init__(self, coord_path, special_prim=None, special_charge=None, special_spin=None):
        self.coord_path = coord_path
    #def __init__(self, mol_class=str()):
        #number of atoms
        self.natoms = int() 
        #atom coords with index and element, type:list
        self.atomtable = []
        #table with labels of integers for atoms, starting at 0, type: np array
        self.bond_lables = []
        #table with just bond orders, type: np array
        self.bond_order = []
        #table with labels and bond order, type:numpy.array
        self.bond_table = []
        #atom connectivity matrix
        self.A = []
        #list of primiatives
        self.prims = []
        #prim connect chart
        self.primchart = []
        #higher order prim conn chart
        self.molchart = []
        self.covrad = cov_rad.form_covalent_radii()
        self.optxyz = []
        self.prim_dist = []
        self.special_prim = special_prim
        self.special_charge = special_charge
        self.special_spin = special_spin
        
    def initalize_molecule(self):
        self.parse_cml(self.coord_path)

    def parse_cml(self, filename):
        """ Finds file location, runs parse_cml(), runs build_matrix()
        
        This the only funciton that needs to be called within the Molecule() class
        
        Parameters
        ----------
        file_name : str
            This is the molecule name as it appears in its .cml file
        
        """

        tree = ET.parse(filename)
        root = tree.getroot()
        #following is to find atomarray and bondarray for different formatted cml files
        rows = []
        for node in root:
            rows.append(str(node.tag))
        for i in range(0, len(rows)):
            if rows[i].endswith('atomArray'):
                atomArray = root[i]
            if rows[i].endswith('bondArray'):
                bondArray = root[i]
        
        self.natoms = len(atomArray)
        self.A = np.zeros((self.natoms,self.natoms)) 
        self.atomtable = [[0 for x in range(4)] for y in range(self.natoms)]
        for atomi in range(0, self.natoms):
            self.atomtable[atomi][0] = atomArray[atomi].attrib['elementType']
            self.atomtable[atomi][1] = float(atomArray[atomi].attrib['x3'])
            self.atomtable[atomi][2] = float(atomArray[atomi].attrib['y3'])
            self.atomtable[atomi][3] = float(atomArray[atomi].attrib['z3'])
        #start extracting bond order
        for bondi in bondArray:
            a12 = bondi.attrib['atomRefs2'].split()
            #assert is making all first part of indexes = a
            assert(a12[0][0] == "a")
            assert(a12[1][0] == "a")
            #takes a1 and makes it a0 so index of element is same in bondorder
            a1 = np.array([int(a12[0][1:])-1])
            a2 = np.array([int(a12[1][1:])-1])
            self.bond_labels = np.append(a1,a2)
            self.bond_order = np.array([int(bondi.attrib["order"])])
            self.bond_table = np.append(self.bond_labels, self.bond_order)
            #these need to be within for loop, don't know if need self.
            x = self.bond_table[0]
            y = self.bond_table[1]
            z = self.bond_table[2]
            self.A[x][y] = z
            self.A[y][x] = z
        self.build_prims()
        
        #pandas testing
        #path = "../inputs/" + 'aspirin.json'
        #dataframe = pd.read_json(path, orient='index')
        #print(dataframe)
        #arr = dataframe.to_numpy()
        #df_cols = ['atomArray', 'bondArray']
    
    def build_prims(self):
        """ Builds the primiative fragments (smallest possible fragments)
        
        Primitives are initalized with heavy atoms (all other than hydrogen), no double or 
        triple bonds are cut, no hydrogen bonds are cut.
        
        Parameters
        ----------
        none
        
        Returns
        ------
        self.prims : set
            This is a set of sets containing the atom indexes that are within each primitive
        
        """
        prims = []
        for i in range(0, len(self.atomtable)):
            if self.atomtable[i][0] != "H":
                #initalizing prim at every heavy atom
                prims.append([i])
                for j in range(0, len(self.atomtable)):
                    #append new atom if H-bond or higher order bond
                    if self.A[i][j] == 1 and self.atomtable[j][0] == "H" or self.A[i][j] > 1:
                        prims[-1].append(j)
                    for k in range(0, len(self.atomtable)):
                        #if atom j is bonded to an H or higher order bond, append atom k also
                        if self.A[i][j] > 1 and self.A[j][k] == 1 and self.atomtable[k][0] == "H":
                                prims[-1].append(k)

        #deleteing prims that are fully contained within a larger prim
        prim = []
        primlist = [set(i) for i in prims]
        
        for i in range(0, len(primlist)):
            add = True
            for j in range(0, len(primlist)):
                if primlist[i].issubset(primlist[j]) and i != j and primlist[i] != primlist[j]:# or primlist[i] == primlist[j] and i != j:
                    add = False
            if add == True:
                prim.append(list(primlist[i]))   
        
        #remove repeating but not subset primitives
        prim = list(set(tuple(sorted(sub)) for sub in prim))
        
        self.prims = []
        for z in range(0, len(prim)):
            tempatom = []
            for atom in prim[z]:
                tempatom.append(self.atomtable[atom][0])
            if self.special_prim != None and self.special_prim in tempatom:
                self.prims.append(Primitive.Primitive(z, prim[z], charge=self.special_charge, spin=self.special_spin))
                print(self.prims[-1].charge, self.prims[-1].spin)
                print("this in molecule.py")

            else:
                self.prims.append(Primitive.Primitive(z, prim[z], charge=0, spin=0))
        return self.prims

    def build_atom_dist(self):
        """ Builds a matrix with atom to atom distances

        Returns
        -------
        atom_dist = np.ndarray
            Ndarray of shape(# atoms, # atoms) with the Euclidean distance as the elements.
        """

        atom_dist = np.zeros((self.natoms, self.natoms))
        for i in range(0, len(self.atomtable)):
            for j in range(0, len(self.atomtable)):
                if i < j:
                    x = np.array(self.atomtable[i][1:])
                    y = np.array(self.atomtable[j][1:])
                    dist = np.linalg.norm(x - y)
                    atom_dist[i][j] = dist
                    atom_dist[j][i] = dist
        return atom_dist

    def build_prim_dist(self):
        """ Builds a matrix with prim to prims distances.

        Distance is determined by the minimum interatomic distances between pairs of atoms in their respective prims.
        """
        atom_dist = self.build_atom_dist()
        self.prim_dist = np.zeros((len(self.prims), len(self.prims)))
        
        for primi in range(0, len(self.prims)):
            for primj in range(primi+1, len(self.prims)):
                dist_list = []
                for atomi in self.prims[primi].atoms:
                    arr = np.take(atom_dist[atomi], list(self.prims[primj].atoms))
                    dist_list.append(np.min(arr))
                x = np.min(dist_list)
                self.prim_dist[primi][primj] = np.min(dist_list)
                self.prim_dist[primj][primi] = np.min(dist_list)
        return self.prim_dist
    
    def build_primchart(self):
        """ Builds a connectivity ndarray between primitaves
        
        Ndarray contains only single-linkage connectivity.  Thus if one primative is 
        connected to another primitave there will be a 1 in the row/column index of the array.
        This is used in the covalent network fragmentation scheme and not spheric scheme.
       
       Parameters
        ----------
        none
        
        Returns
        -------
        self.primchart : ndarray
            Ndarray of shape (# of primitaves, # of primitaves)
        
        """
        
        self.primsleng = len(self.prims)
        self.primchart = np.zeros( (self.primsleng,self.primsleng))
        for prim1 in range(0, len(self.prims)):
            for atomi in self.prims[prim1].atoms:
                for prim2 in range(0, len(self.prims)):
                    for atomj in self.prims[prim2].atoms:
                        if prim1 == prim2:
                            continue
                        if self.A[atomi][atomj] != 0:
                            self.primchart[prim1][prim2] = 1
                            self.primchart[prim2][prim1] = 1
        return self.primchart

    def build_molmatrix(self, i):
        """ Builds the full connectivity array
        
        Ndarray contains all connectivity between primitives. Thus row index and column index entry
        correspond to number of primiatives away column prim is from row prim. 
        This is used in the covalent network fragmentation scheme and not spheric scheme.
        
        Built through a series of matrix manipulations using self.primchart and self.molchart
        
        Parameters
        ----------
        i : int
            This must be the value 2 to start the recursive funciton for the matrix manipulations. This
            is set within the initalize_molecule() funciton
        
        Returns
        -------
        self.molchart : ndarray
            Ndarray of shape (# of prims, # of prims)
        
        """
        
        self.build_primchart()
        eta = self.natoms
        if i == 2:
            self.molchart = self.primchart.dot(self.primchart)
            np.fill_diagonal(self.molchart, 0)
            for x in range(0, len(self.prims)):
                for y in range(0, len(self.prims)):
                    if self.molchart[x][y] != 0:
                        self.molchart[x][y] = i
            
            for x in range(0, len(self.prims)):      #checkpoint
                for y in range(0, len(self.prims)):
                    if self.primchart[x][y] != 0:
                        self.molchart[x][y] = 0
            self.molchart = np.add(self.molchart, self.primchart)
        self.molchartnew = []   
        if i > 2:
            self.molchartnew = self.molchart.dot(self.primchart)
            np.fill_diagonal(self.molchartnew, 0)
            for x in range(0, len(self.prims)):
                for y in range(0, len(self.prims)):
                    if self.molchartnew[x][y] != 0:
                        self.molchartnew[x][y] = i
            
            for x in range(0, len(self.prims)):     #checkpoint 
                for y in range(0, len(self.prims)):
                    if self.molchart[x][y] != 0:
                        self.molchartnew[x][y] = 0
            self.molchart = np.add(self.molchart, self.molchartnew)

        if i < eta:     #recursive part of function
            i = i+1
            self.build_molmatrix(i)
        return self.molchart

if __name__ == "__main__":
    carbonylavo = Molecule()
    carbonylavo.initalize_molecule('carbonylavo')
    carbonylavo.build_prim_dist()
