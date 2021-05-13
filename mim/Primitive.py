import numpy as np

class Primitive():
    """
    Class to store individual primitives each with their own charge and spin multiplicity, which atoms are in primitive and label of primitive 
"""

    def __init__(self, prim_label, atoms, charge=0, spin=1):
        self.prim_label = prim_label
        self.atoms = atoms
        self.charge = charge
        self.spin = spin



