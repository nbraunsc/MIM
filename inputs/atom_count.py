import os
import numpy as np


def test_atoms(atom_list, coeff_list, num_atoms):
    vec = np.zeros((num_atoms))
    for frag in range(0, len(atom_list)):
        print(atom_list[frag])
        print(coeff_list[frag])
        for atom in atom_list[frag]:
            print(atom)
            print(coeff_list[frag])
            vec[atom] += coeff_list[frag]
    return vec
